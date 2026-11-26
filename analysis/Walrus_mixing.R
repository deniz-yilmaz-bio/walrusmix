library("MixSIAR")


#No grouping factors or covariates in the dataset, therefore the NULL arguements

consumer <- load_mix_data(filename = "walrus_consumer.csv",
                          iso_names = c("d13C", "d15N"),
                          factors = NULL,
                          fac_random = NULL,
                          fac_nested = NULL,
                          cont_effects = NULL)

#No diverse sources. Therefore conc_dep = FALSE. Data reported as means/SD/n

source <- load_source_data(filename = "walrus_sources.csv",
                           source_factors = NULL,
                           conc_dep = FALSE,
                           data_type = "means",
                           consumer)

tdf <- load_discr_data(filename = "walrus_discrimination_main.csv", consumer)

#Plotting the samples and TDF-shifted sources in the isospcae to see if they nest
#within source distributions.Initial simple plot for overview.

plot_data(filename = "isospace_plot",plot_save_pdf = FALSE, plot_save_png = TRUE,
          consumer, source, tdf)

calc_area(source = source, mix = consumer, discr = tdf)

#Using ggplot to replicate the MixSIAR isospace plot with defined polygon edges
#to see which sources actually fall inside.

library(ggplot2)

src  <- read.csv("walrus_sources.csv")
discr  <- read.csv("walrus_discrimination_main.csv")
cons <- read.csv("walrus_consumer.csv")

#Making sure column name for the sources are the same between src and discr

names(src)[1] <- "source"
names(discr)[1] <- "source"

#Extra security check for column order match.

discr <- discr[match(src$source, discr$source), ]
stopifnot(all(src$source == discr$source))

#Shifting the source means and SDs with TDFs. 

shifted <- data.frame(
  source = src$source,
  d13C   = src$Meand13C + discr$Meand13C,
  d15N   = src$Meand15N + discr$Meand15N,
  sd13C  = sqrt(src$SDd13C^2 + discr$SDd13C^2),
  sd15N  = sqrt(src$SDd15N^2 + discr$SDd15N^2)
)

#The chull function in this case returns the source organisms that would form the
#corners of the polygon. By using its output to define a hull variable, we can draw
#the polygon without using the sources that fall inside it.

hull_corners <- chull(shifted$d13C, shifted$d15N)
hull <- shifted[hull_corners, ]
hull$source   

#Plotting of the sources, errorbars, consumers, and the isospace polygon.

iso_x <- expression(delta^13*C~"(‰)")
iso_y <- expression(delta^15*N~"(‰)")

polygon.plot <- ggplot() +
  geom_polygon(data = hull, aes(d13C, d15N),
               fill = "steelblue", alpha = 0.10,
               colour = "steelblue", linewidth = 0.8) +
  geom_errorbar(data = shifted, aes(x = d13C, ymin = d15N - sd15N, ymax = d15N + sd15N),
                width = 0, colour = "grey30") +
  geom_errorbar(data = shifted, aes(y = d15N, xmin = d13C - sd13C, xmax = d13C + sd13C),
                width = 0, colour = "grey30", orientation = "y") +
  geom_point(data = shifted, aes(d13C, d15N), shape = 18, size = 4) +
  geom_text(data = shifted, aes(d13C, d15N, label = source),
            nudge_x = 0.12, nudge_y = 0.25, hjust = 0, fontface = "bold", size = 3.5) +
  geom_point(data = cons, aes(d13C, d15N), colour = "firebrick", size = 2) +
  coord_equal() +    # same scale on both axes, so distances look right
  labs(x = iso_x, y = iso_y,
       title = "Mixing polygon (TDF-shifted source means)",
       subtitle = "Diamonds: sources ± 1 SD   Red: walrus mean") +
  theme_bw()

polygon.plot

#Plotting of source percentage distribution with uninformative prirors. qBeta call
#to get minima and maxima within the confidence interval, as well as the median
#in the case of a Dirichlet prior with 4 sources.

plot_prior(alpha.prior = 1, source)

qbeta(c(0.025,0.5,0.975),1,3)

#Model writing and defining different run and output options.

model.name <- "Walrus_mixingModel.txt"
write_JAGS_model(model.name, resid_err = FALSE, process_err = TRUE, consumer, source)

jags.uninf.test <- run_model(run = "test", consumer, source, tdf, model.name, alpha.prior = 1,
                             process_err = TRUE, resid_err = FALSE)

jags.uninf.normal <- run_model(run = "normal", consumer, source, tdf, model.name, alpha.prior = 1,
                               process_err = TRUE, resid_err = FALSE)

jags.uninf.long <- run_model(run = "long", consumer, source, tdf, model.name, alpha.prior = 1,
                               process_err = TRUE, resid_err = FALSE)

output_options_normal <- list(
  summary_save = TRUE,  summary_name = "summary_statistics_normal",
  sup_post = FALSE,     plot_post_save_pdf = FALSE, plot_post_name = "posterior_density_normal",
  sup_pairs = FALSE,    plot_pairs_save_pdf = FALSE, plot_pairs_name = "pairs_plot_normal",
  sup_xy = TRUE,        plot_xy_save_pdf = FALSE,    plot_xy_name = "xy_plot_normal",
  gelman = TRUE, heidel = FALSE, geweke = TRUE,
  diag_save = TRUE,     diag_name = "diagnostics_normal",
  indiv_effect = FALSE,
  plot_post_save_png = TRUE, plot_pairs_save_png = TRUE, plot_xy_save_png = TRUE,
  diag_save_ggmcmc = FALSE, return_obj = TRUE
)

output_options_long <- list(
  summary_save = TRUE,  summary_name = "summary_statistics_long",
  sup_post = FALSE,     plot_post_save_pdf = FALSE, plot_post_name = "posterior_density_long",
  sup_pairs = FALSE,    plot_pairs_save_pdf = FALSE, plot_pairs_name = "pairs_plot_long",
  sup_xy = TRUE,        plot_xy_save_pdf = FALSE,    plot_xy_name = "xy_plot_long",
  gelman = TRUE, heidel = FALSE, geweke = TRUE,
  diag_save = TRUE,     diag_name = "diagnostics_long",
  indiv_effect = FALSE,
  plot_post_save_png = TRUE, plot_pairs_save_png = TRUE, plot_xy_save_png = TRUE,
  diag_save_ggmcmc = FALSE, return_obj = TRUE
)

output_JAGS(jags.uninf.normal, consumer, source, output_options_normal)

output_JAGS(jags.uninf.long, consumer, source, output_options_long)


#The previous process repeated by using synthetic individual walrus data points
#generated from the reported summary mean and standard deviation.

consumersynth <- load_mix_data(filename = "walrus_consumer_n4_synthetic.csv",
                          iso_names = c("d13C", "d15N"),
                          factors = NULL,
                          fac_random = NULL,
                          fac_nested = NULL,
                          cont_effects = NULL)

plot_data(filename = "isospace_plot_synth",plot_save_pdf = FALSE, plot_save_png = TRUE,
          consumersynth, source, tdf)

conssynth <- read.csv("walrus_consumer_n4_synthetic.csv")

polygon.plot.synth <- ggplot() +
  geom_polygon(data = hull, aes(d13C, d15N),
               fill = "steelblue", alpha = 0.10,
               colour = "steelblue", linewidth = 0.8) +
  geom_errorbar(data = shifted, aes(x = d13C, ymin = d15N - sd15N, ymax = d15N + sd15N),
                width = 0, colour = "grey30") +
  geom_errorbar(data = shifted, aes(y = d15N, xmin = d13C - sd13C, xmax = d13C + sd13C),
                width = 0, colour = "grey30", orientation = "y") +
  geom_point(data = shifted, aes(d13C, d15N), shape = 18, size = 4) +
  geom_text(data = shifted, aes(d13C, d15N, label = source),
            nudge_x = 0.12, nudge_y = 0.25, hjust = 0, fontface = "bold", size = 3.5) +
  geom_point(data = conssynth, aes(d13C, d15N), colour = "firebrick", size = 2) +
  coord_equal() +    # same scale on both axes, so distances look right
  labs(x = iso_x, y = iso_y,
       title = "Mixing polygon (TDF-shifted source means)",
       subtitle = "Diamonds: sources ± 1 SD   Red: Synthetic Walrus Points") +
  theme_bw()

polygon.plot.synth

model.name.synth <- "Walrus_mixingModel_synth.txt"
write_JAGS_model(model.name.synth, resid_err = TRUE, process_err = TRUE, consumersynth, source)

jags.uninf.synth.test <- run_model(run = "test", consumersynth, source, tdf, model.name.synth, alpha.prior = 1,
                             process_err = TRUE, resid_err = TRUE)

jags.uninf.synth.normal <- run_model(run = "normal", consumersynth, source, tdf, model.name.synth, alpha.prior = 1,
                               process_err = TRUE, resid_err = TRUE)

jags.uninf.synth.long <- run_model(run = "long", consumersynth, source, tdf, model.name.synth, alpha.prior = 1,
                                     process_err = TRUE, resid_err = TRUE)


output_options_normal_synth <- list(
  summary_save = TRUE,  summary_name = "summary_statistics_normal_synth",
  sup_post = FALSE,     plot_post_save_pdf = FALSE, plot_post_name = "posterior_density_normal_synth",
  sup_pairs = FALSE,    plot_pairs_save_pdf = FALSE, plot_pairs_name = "pairs_plot_normal_synth",
  sup_xy = TRUE,        plot_xy_save_pdf = FALSE,    plot_xy_name = "xy_plot_normal_synth",
  gelman = TRUE, heidel = FALSE, geweke = TRUE,
  diag_save = TRUE,     diag_name = "diagnostics_normal_synth",
  indiv_effect = FALSE,
  plot_post_save_png = TRUE, plot_pairs_save_png = TRUE, plot_xy_save_png = TRUE,
  diag_save_ggmcmc = FALSE, return_obj = TRUE
)

output_options_long_synth <- list(
  summary_save = TRUE,  summary_name = "summary_statistics_long_synth",
  sup_post = FALSE,     plot_post_save_pdf = FALSE, plot_post_name = "posterior_density_long_synth",
  sup_pairs = FALSE,    plot_pairs_save_pdf = FALSE, plot_pairs_name = "pairs_plot_long_synth",
  sup_xy = TRUE,        plot_xy_save_pdf = FALSE,    plot_xy_name = "xy_plot_long_synth",
  gelman = TRUE, heidel = FALSE, geweke = TRUE,
  diag_save = TRUE,     diag_name = "diagnostics_long_synth",
  indiv_effect = FALSE,
  plot_post_save_png = TRUE, plot_pairs_save_png = TRUE, plot_xy_save_png = TRUE,
  diag_save_ggmcmc = FALSE, return_obj = TRUE
)

output_JAGS(jags.uninf.synth.normal, consumersynth, source, output_options_normal_synth)

output_JAGS(jags.uninf.synth.long, consumersynth, source, output_options_long_synth)

#The previous steps repeated for a model where the common d15N value of 3.4 is 
#used instead of the Hobson value.

discr_post <- data.frame(
  Sources = discr$source,
  Meand13C   = discr$Meand13C,
  SDd13C  = discr$SDd13C,
  Meand15N   = 3.4,
  SDd15N  = 1
  )

write.csv(discr_post,file = "walrus_discrimination_post.csv", row.names = FALSE)

tdf_post <- load_discr_data("walrus_discrimination_post.csv", consumer)

plot_data(filename = "isospace_plot_tdfpost",plot_save_pdf = FALSE, plot_save_png = TRUE,
          consumer, source, tdf_post)

jags.uninf.post.test <- run_model(run = "test", consumer, source, tdf_post, 
                                  model.name, alpha.prior = 1, 
                                  process_err = TRUE, resid_err = FALSE)

jags.uninf.post.normal <- run_model(run = "normal", consumer, source, tdf_post, 
                                    model.name, alpha.prior = 1, 
                                    process_err = TRUE, resid_err = FALSE)

jags.uninf.post.long <- run_model(run = "long", consumer, source, tdf_post, 
                                  model.name, alpha.prior = 1, 
                                  process_err = TRUE, resid_err = FALSE)

output_options_normal_post <- list(
  summary_save = TRUE,  summary_name = "summary_statistics_normal_post",
  sup_post = FALSE,     plot_post_save_pdf = FALSE, plot_post_name = "posterior_density_normal_post",
  sup_pairs = FALSE,    plot_pairs_save_pdf = FALSE, plot_pairs_name = "pairs_plot_normal_post",
  sup_xy = TRUE,        plot_xy_save_pdf = FALSE,    plot_xy_name = "xy_plot_normal_post",
  gelman = TRUE, heidel = FALSE, geweke = TRUE,
  diag_save = TRUE,     diag_name = "diagnostics_normal_post",
  indiv_effect = FALSE,
  plot_post_save_png = TRUE, plot_pairs_save_png = TRUE, plot_xy_save_png = TRUE,
  diag_save_ggmcmc = FALSE, return_obj = TRUE
)

output_options_long_post <- list(
  summary_save = TRUE,  summary_name = "summary_statistics_long_post",
  sup_post = FALSE,     plot_post_save_pdf = FALSE, plot_post_name = "posterior_density_long_post",
  sup_pairs = FALSE,    plot_pairs_save_pdf = FALSE, plot_pairs_name = "pairs_plot_long_post",
  sup_xy = TRUE,        plot_xy_save_pdf = FALSE,    plot_xy_name = "xy_plot_long_post",
  gelman = TRUE, heidel = FALSE, geweke = TRUE,
  diag_save = TRUE,     diag_name = "diagnostics_long_post",
  indiv_effect = FALSE,
  plot_post_save_png = TRUE, plot_pairs_save_png = TRUE, plot_xy_save_png = TRUE,
  diag_save_ggmcmc = FALSE, return_obj = TRUE
)

output_JAGS(jags.uninf.post.normal, consumer, source, output_options_normal_post)

output_JAGS(jags.uninf.post.long, consumer, source, output_options_long_post)

#Nereis and Hiatella combined for both the initial model and the TDF d15N = 3.4
#model. Probability checks for Hiatella dominance and interval plots


groups_HN <- list(Hiatella_Nereis = c("Hiatella", "Nereis"),
                  Molpadia = "Molpadia",
                  Buccinum = "Buccinum")


combnormal   <- combine_sources(jags.uninf.long, consumer, source, alpha.prior = 1, groups = groups_HN)
combpost <- combine_sources(jags.uninf.post.long, consumer, source, alpha.prior = 1, groups = groups_HN)

combnormal$prior.new
qbeta(c(0.025,0.5,0.975),2,2)


summary_stat(combnormal,   toprint = "global", meanSD = FALSE,
             quantiles = c(.025, .5, .975), savetxt = TRUE, filename = "summary_combined_Hobson")
summary_stat(combpost, toprint = "global", meanSD = FALSE,
             quantiles = c(.025, .5, .975), savetxt = TRUE, filename = "summary_combined_Post")


plot_intervals(combnormal, toplot = "global", prob = 0.5)
ggplot2::ggsave("intervals_combined_Hobson.png", plot = ggplot2::last_plot(),
                width = 6, height = 4, dpi = 300)

plot_intervals(combpost, toplot = "global", prob = 0.5)
ggplot2::ggsave("intervals_combined_post.png", plot = ggplot2::last_plot(),
                width = 6, height = 4, dpi = 300)


prob_checks <- function(jags_obj, model_name) {
  p <- jags_obj$BUGSoutput$sims.list$p.global
  colnames(p) <- source$source_names
  H <- p[, "Hiatella"]; N <- p[, "Nereis"]
  data.frame(
    model                 = model_name,
    `P(H > N)`            = round(mean(H > N), 3),
    `P(H > 0.5)`          = round(mean(H > 0.5), 3),
    `P(H largest source)` = round(mean(apply(p, 1, which.max) ==
                                         which(colnames(p) == "Hiatella")), 3),
    check.names = FALSE
  )
}

prob_table <- rbind(prob_checks(jags.uninf.long, "(D15N 2.4)"),
                    prob_checks(jags.uninf.post.long, "(D15N 3.4)"))
prob_table
write.csv(prob_table, "hiatella_probability_checks.csv", row.names = FALSE)

# P(H > N):            probability walrus eat more Hiatella than Nereis
# P(H > 0.5):          probability Hiatella alone is over half the diet
# P(H largest source): probability Hiatella is the single biggest source
