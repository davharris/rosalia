library(dplyr)
library(mgcv)
library(ggplot2)



x = read.csv("estimates.csv", stringsAsFactors = FALSE)
x$dynamic_simulation = grepl("pop", x$rep_name)
x$simulation_type = gsub("[0-9]", "", x$rep_name)

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


pdf("manuscript-materials/figures/error_rates.pdf", width = 5, height = 5)
plot(p_ns[order(xx$truth)] ~ sort(xx$truth), type = "l", lwd = 3, ylim = c(0, 1), yaxs = "i", bty = "l", ylab = "Proportion", las = 1, xlab = "\"True\" coefficient", xaxs = "i")
lines(p_positive[order(xx$truth)] ~ sort(xx$truth), lwd = 3, col = 4)
lines(p_negative[order(xx$truth)] ~ sort(xx$truth), lwd = 3, col = 2)
text(0, .95, "Not significant\nat the 0.05 level")
text(5, .8, "Significantly\npositive", col = 4)
text(-7, .35, "Significantly\nnegative", col = 2)
dev.off()


library(tidyr)
z = x %>%
  filter(method %in% c("correlation", "BayesComm") & simulation_type == "no_env") %>%
  dplyr::select(-lower, -upper, -X) %>%
  spread(method, estimate)

ggplot(z, aes(x = correlation, y = BayesComm, col = factor(n_sites))) + 
  geom_point()



