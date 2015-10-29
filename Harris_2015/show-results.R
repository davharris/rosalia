library(dplyr)
library(mgcv)
library(ggplot2)
library(tidyr)


x = read.csv("estimates.csv", stringsAsFactors = FALSE)
x$dynamic_simulation = grepl("pop", x$rep_name)
x$simulation_type = gsub("[0-9]", "", x$rep_name)


######## Pairs ########

pairs_txt = readLines("fakedata/matrices/Pairs.txt")
library(stringr)

# Find areas of the data file that correspond
# to species pairs' results
beginnings = grep("Sp1", pairs_txt) + 1
ends = c(
  grep("^[^ ]", pairs_txt)[-1],
  length(pairs_txt)
) - 1

file_list = readLines("manuscript-materials/fakedata_list.txt")[-1]
partial_names = sapply(strsplit(grep("^>", pairs_txt, value = TRUE), " +"), function(x) x[[3]])
filename_lines = grep("^>", pairs_txt)

# File_list agrees with partial_names to the extent possible
all(sapply(1:450, function(i)grepl(partial_names[i], file_list[i])))


split_names = strsplit(file_list, "-")

pairs_results = lapply(
  1:length(file_list),
  function(i){
    
    message(i)
    
    n_sites = as.integer(split_names[[i]][[2]])
    rep_name = split_names[[i]][[3]]
    
    # Find the line where the current data set is mentioned in
    # pairs.txt
    filename_line = filename_lines[i]
    
    # Which chunk of the data file corresponds to this file?
    chunk = min(which(beginnings > filename_line))
    
    # Split the chunk on whitespace.
    splitted = strsplit(pairs_txt[beginnings[chunk]:ends[chunk]], " +")
    
    # Pull out the corresponding chunk of the "x" data frame, based on n_sites and rep_name
    x_subset = x[x$n_sites == n_sites & x$rep_name == rep_name & x$method == "correlation", ]
    
    # Pull out the species numbers and their Z-scores, then join to x_subset
    pairs_results = lapply(
      splitted,
      function(x){
        # in the x data frame, species 1 is always a lower number than species 2
        spp = sort(as.integer(x[3:4]))
        data.frame(
          sp1 = paste0("V", spp[1]), 
          sp2 = paste0("V", spp[2]), 
          z = as.numeric(x[14]),
          stringsAsFactors = FALSE
        )
      }
    ) %>% 
      bind_rows %>%
      mutate(spp = paste(sp1, sp2, sep = "-"))
    
    n_spp = length(unique(c(x_subset$sp1, x_subset$sp2)))
    
    # Re-order the pairs_results to match the other methods
    m = matrix(NA, n_spp, n_spp)
    new_order = match(
      paste0("V", row(m)[upper.tri(m)], "-V", col(m)[upper.tri(m)]),
      pairs_results$spp
    )
    ordered_pairs_results = pairs_results[new_order, ]
    
    x_subset$estimate = ordered_pairs_results$z
    x_subset$method = "null"
    
    x_subset
  }
) %>% bind_rows()

pairs_results

x = rbind(x, pairs_results)


#####



x = x[!is.na(x$estimate), ]

resids = function(data){
  resid(lm(truth ~ estimate + 0, data = data))
}

result_summary = x %>% 
  group_by(method, simulation_type) %>%
  do(data.frame(., resids = resids(.))) %>%
  ungroup %>%
  group_by(method, simulation_type, n_sites) %>%
  summarise(r2 = 1 - sum(resids^2) / sum(truth^2))

result_summary$method = reorder(result_summary$method, -result_summary$r2)

result_summary$simulation_type = reorder(result_summary$simulation_type, -result_summary$r2)


result_summary %>% 
  filter(simulation_type == "no_env") %>%
  group_by(method) %>% 
  summarise(round(mean(r2), 3))

result_summary %>% 
  filter(simulation_type == "env") %>%
  group_by(method) %>% 
  summarise(round(mean(r2), 3))

result_summary %>% 
  filter(simulation_type == "pop") %>%
  group_by(method) %>% 
  summarise(round(mean(r2), 3))



pdf("manuscript-materials/figures/performance.pdf", width = 7.5, height = 2.5)
ggplot(result_summary, aes(x = n_sites, y = r2, col = method)) + 
  facet_grid(~simulation_type) + 
  geom_line(size = .5) + 
  geom_point(size = 2, fill = "white") + 
  geom_hline(yintercept = 0, size = 1/2) +
  geom_vline(xintercept = 0, size = 1) + 
  #scale_shape_manual(values = c(21, 25, 15, 18)) + 
  coord_cartesian(ylim = c(-.01, 0.76)) + 
  ylab(expression(R^2)) + 
  xlab("Number of sites (log scale)") + 
  scale_x_log10(breaks = unique(x$n_sites), limits = range(x$n_sites)) +
  theme_bw(base_size = 11) + 
  theme(panel.margin = grid::unit(1.25, "lines")) + 
  theme(panel.border = element_blank(), axis.line = element_blank()) + 
  theme(
    panel.grid.minor = element_blank(), 
    panel.grid.major.y = element_line(color = "lightgray", size = 1/4), 
    panel.grid.major.x = element_blank()
  ) + 
  theme(strip.background = element_blank(), legend.key = element_blank()) + 
  theme(plot.margin = grid::unit(c(.01, .01, .75, .1), "lines")) + 
  theme(
    axis.title.x = element_text(vjust = -0.2, size = 14), 
    axis.title.y = element_text(angle = 0, hjust = -.1, size = 14)
  ) + 
  scale_color_brewer(palette = "Dark2")
dev.off()


####
library(lme4)
library(multcomp)

landscape_estimates = x %>% 
  group_by(method, simulation_type) %>%
  do(data.frame(., resids = resids(.))) %>%
  ungroup %>%
  group_by(method, simulation_type, rep_name, n_sites) %>%
  summarise(r2 = 1 - sum(resids^2) / sum(truth^2))


summary(lmer(r2 ~ method + n_sites + simulation_type + (1|rep_name), data = landscape_estimates))

####


xx = x[x$method == "Markov network" & grepl("^env[0-9]*$", x$rep_name), ]

p_positive = predict(gam(I(lower > 0) ~ s(truth), family = binomial, data = xx), type = "response")
p_negative = predict(gam(I(upper < 0) ~ s(truth), family = binomial, data = xx), type = "response")
p_ns = predict(gam(I(lower < 0 & upper > 0) ~ s(truth), family = binomial, data = xx), type = "response")

xx2 = x[x$method == "null" & grepl("^env[0-9]*$", x$rep_name), ]

p_positive2 = predict(gam(I(estimate < pnorm(.025)) ~ s(truth), family = binomial, data = xx2), type = "response")
p_negative2 = predict(gam(I(estimate > pnorm(.975)) ~ s(truth), family = binomial, data = xx2), type = "response")
p_ns2 = predict(gam(I(abs(estimate) < pnorm(.975)) ~ s(truth), family = binomial, data = xx2), type = "response")



pdf("manuscript-materials/figures/error_rates.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(p_ns[order(xx$truth)] ~ sort(xx$truth), type = "l", lwd = 3, ylim = c(0, 1), yaxs = "i", bty = "l", ylab = "Proportion", las = 1, xlab = "\"True\" coefficient", xaxs = "i")
lines(p_positive[order(xx$truth)] ~ sort(xx$truth), lwd = 3, col = 4)
lines(p_negative[order(xx$truth)] ~ sort(xx$truth), lwd = 3, col = 2)
text(0, .95, "Not significant\nat the 0.05 level")
text(5, .8, "Significantly\npositive", col = 4)
text(-7, .35, "Significantly\nnegative", col = 2)

plot(p_ns2[order(xx2$truth)] ~ sort(xx2$truth), type = "l", lwd = 3, ylim = c(0, 1), yaxs = "i", bty = "l", ylab = "Proportion", las = 1, xlab = "\"True\" coefficient", xaxs = "i")
lines(p_positive2[order(xx2$truth)] ~ sort(xx2$truth), lwd = 3, col = 4)
lines(p_negative2[order(xx2$truth)] ~ sort(xx2$truth), lwd = 3, col = 2)
text(0, .95, "Not significant\nat the 0.05 level")
text(5, .8, "Significantly\npositive", col = 4)
text(-7, .35, "Significantly\nnegative", col = 2)
dev.off()


plot((p_positive / (p_positive + p_negative))[order(xx$truth)] ~ sort(xx$truth), type = "l")



z = x %>%
  dplyr::select(-lower, -upper, -X) %>%
  spread(method, estimate)

ggplot(z, aes(x = correlation, y = null, col = factor(n_sites))) + 
  geom_point()

ggplot(z, aes(x = correlation, y = BayesComm, col = factor(n_sites))) + 
  geom_point()


summary(lm(null ~ correlation * factor(n_sites), data = z))

summary(lm(BayesComm ~  correlation * factor(n_sites), data = z))

summary(lm(`Markov network` ~  correlation * factor(n_sites), data = z))
