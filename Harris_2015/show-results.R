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


####

summary(lmer(r2 ~ method + n_sites + simulation_type + (1|rep_name), data = landscape_estimates))

pairs_summary = x[x$method == "null" & !grepl("pop", x$rep_name), ]
markov_summary = x[x$method == "Markov network" & !grepl("pop", x$rep_name), ]


# Pairs's Z score is based on C-scores, which are positive when species are 
# disaggregated.  So significantly positive interactions == negative estimates
pairs_summary$sig_pos = pairs_summary$estimate < qnorm(.025)
pairs_summary$sig_neg = pairs_summary$estimate > qnorm(.975)

# Markov network is significant when lower bound is above zero or lower bound is
# below zero
markov_summary$sig_pos = markov_summary$lower > 0
markov_summary$sig_neg = markov_summary$upper < 0


my_geom_smooth = function(data, color, ...){
  geom_smooth(
    data = data,
    aes(
      x = truth, 
      y = as.integer((sig_neg & truth > 0) | (sig_pos & truth < 0)),
      color = color
    ),
    method = gam, 
    family = binomial, 
    formula = y ~ s(x),
    se = FALSE,
    n = 1024,
    ...
  )
}


truth_seq = seq(min(data$truth), max(data$truth), length = 1000)

plotfun = function(..., add){
  if(add){
    lines(...)
  }else{
    plot(...)
  }
}


error_smoother = function(data){
  predict(
    gam(
      I((sig_neg & truth > 0) | (sig_pos & truth < 0)) ~ s(truth),
      data = data,
      family = binomial
    ),
    data.frame(truth = truth_seq),
    type = "response"
  )
}

y_pairs = error_smoother(pairs_summary)
y_markov = error_smoother(markov_summary)


null_cor = x %>%
  dplyr::select(-lower, -upper, -X) %>%
  spread(method, estimate) %>%
  na.omit()

# R-squared
round(summary(lm(null ~ I(correlation*sqrt(n_sites)), data = null_cor))$r.squared, 2)


pdf("manuscript-materials/figures/error_rates.pdf", height = 8, width = 4)
par(mfrow = c(2, 1))
plotfun(
  truth_seq,
  y_pairs,
  type = "l",
  xlab = "\"True\" interaction strength",
  ylab = "P(confidently predict wrong sign)",
  bty = "l",
  yaxs = "i",
  col = 2,
  add = FALSE,
  ylim = c(0, .2),
  lwd = 2
)
mtext("A. Error rate vs. interaction strength", side = 3, adj = 0, font = 2, line = 1.2) 
plotfun(truth_seq, y_markov, add = TRUE, lwd = 2)
legend("topleft", lwd = 2, legend = c("Null model", "Markov network"), col = c(2, 1), bty = "n", cex = .75)


with(
  null_cor, 
  plot(
    correlation * sqrt(n_sites), 
    null, 
    pch = ".", 
    col = "#00000020",
    ylab = "Z-score",
    xlab = expression("correlation" %*% sqrt(number~~of~~sites)),
    bty = "l"
  )
)
mtext("B. Null model estimates vs.\nscaled correlation coefficient", adj = 0, side = 3, font = 2, line = 1.2) 
abline(lm(null ~ I(correlation*sqrt(n_sites)), data = null_cor))
text(10, 12, expression(R^2==.95))

dev.off()


# P(Pairs confidently wrong)
with(pairs_summary,  mean((sig_neg & truth > 0) | (sig_pos & truth < 0)))
# P(Markov network confidently wrong)
with(markov_summary, mean((sig_neg & truth > 0) | (sig_pos & truth < 0)))


# P(reject null)
with(pairs_summary,  mean(sig_neg | sig_pos))
with(markov_summary,  mean(sig_neg | sig_pos))
