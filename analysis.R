#' Analysis conducted using naomi v2.3.14 customised to also
#' output linear predictors for new observations.
#'

devtools::load_all("naomi", export_all = FALSE)

library(tidyverse)
library(scales)
library(sf)
library(gridExtra)
library(cowplot)

#' ## Input data files

areas <- read_area_merged(system.file("extdata/demo_areas.geojson", package = "naomi"))
pop <- read_population(system.file("extdata/demo_population_agesex.csv", package = "naomi"))
survey <- read_survey_indicators(system.file("extdata/survey/demo_survey_hiv_indicators_alt.csv", package = "naomi"))
art <- read_art_number(system.file("extdata/demo_art_number.csv", package = "naomi"))
anc <- read_anc_testing(system.file("extdata/demo_anc_testing.csv", package = "naomi"))
spec <- extract_pjnz_naomi(system.file("extdata/demo_mwi2019.PJNZ", package = "naomi"))



#' ## Select model areas and time points

scope <- "MWI"
level <- 4
calendar_quarter_t1 <- "CY2016Q1"
calendar_quarter_t2 <- "CY2018Q3"
calendar_quarter_t3 <- "CY2019Q3"

#' ### Select data inputs for model fitting

prev_survey_ids  <- c("DEMO2016PHIA", "DEMO2015DHSA")
artcov_survey_ids  <- "DEMO2016PHIA"
vls_survey_ids <- NULL
recent_survey_ids <- "DEMO2016PHIA"

artnum_calendar_quarter_t1 <- "CY2016Q1"
artnum_calendar_quarter_t2 <- "CY2018Q3"

anc_clients_year2 <- 2018
anc_clients_year2_num_months <- 9

anc_prevalence_year1 <- 2016
anc_prevalence_year2 <- 2018

anc_art_coverage_year1 <- 2016
anc_art_coverage_year2 <- 2018


#' ## Setup the model

naomi_mf <- naomi_model_frame(areas,
                              pop,
                              spec,
                              scope = scope,
                              level = level,
                              calendar_quarter_t1,
                              calendar_quarter_t2,
                              calendar_quarter_t3,
                              artattend = TRUE,
                              artattend_t2 = TRUE)


#' Prepare data inputs

naomi_data <- select_naomi_data(naomi_mf,
                                survey,
                                anc,
                                art,
                                prev_survey_ids,
                                artcov_survey_ids,
                                recent_survey_ids,
                                vls_survey_ids,
                                artnum_calendar_quarter_t1,
                                artnum_calendar_quarter_t2,
                                anc_clients_year_t2 = anc_clients_year2,
                                anc_clients_year_t2_num_months = anc_clients_year2_num_months,
                                anc_prevalence_year1,
                                anc_prevalence_year2,
                                anc_art_coverage_year1,
                                anc_art_coverage_year2)


#' Prepare model inputs and initial parameters

tmb_inputs <- prepare_tmb_inputs(naomi_data)


#' ## Fit the TMB model

set.seed(54486050)
  
fit <- fit_tmb(tmb_inputs)
fit <- sample_tmb(fit)
outputs <- output_package(fit, naomi_mf)

indicators <- add_output_labels(outputs) %>%
  left_join(
    select(outputs$meta_area, area_level, area_id, center_x, center_y),
    by = c("area_level", "area_id")
  ) %>%
  sf::st_as_sf() %>%
  mutate(
    area_level_label = fct_reorder(area_level_label, area_level),
    age_group_label = fct_reorder(age_group_label, as.integer(factor(age_group)))
  )

#' ## Posterior predictive sample
#'
#' Note: this is custom for paper analysis; not in main version of model.

fit_pred <- fit
fit_pred$obj <- naomi:::make_tmb_obj(tmb_inputs$data, tmb_inputs$par_init,
                                     calc_outputs = 2L, inner_verbose = FALSE)
fit_pred <- sample_tmb(fit_pred)


#' Prevalence 15-49 by district

prev_pred_x <- rbinom(length(fit_pred$sample$rho_obs_t1),
                      round(naomi_data$prev_dat$n_eff),
                      fit_pred$sample$rho_obs_t1)
dim(prev_pred_x) <- dim(fit_pred$sample$rho_obs_t1)


#' ART coverage 15-64 by district

artcov_pred_x <- rbinom(length(fit_pred$sample$alpha_obs_t1),
                        round(naomi_data$artcov_dat$n_eff),
                        fit_pred$sample$alpha_obs_t1)
dim(artcov_pred_x) <- dim(fit_pred$sample$alpha_obs_t1)


#' ANC prevalence predictive distribution

ancprev_pred_x <- rbinom(length(fit_pred$sample$anc_rho_obs_t1),
                         naomi_data$anc_prev_t1_dat$anc_prev_n,
                         fit_pred$sample$anc_rho_obs_t1)
dim(ancprev_pred_x) <- dim(fit_pred$sample$anc_rho_obs_t1)


#' ANC ART coverage predictive distribution

ancartcov_pred_x <- rbinom(length(fit_pred$sample$anc_alpha_obs_t1),
                           naomi_data$anc_artcov_t1_dat$anc_artcov_n,
                           fit_pred$sample$anc_alpha_obs_t1)
dim(ancartcov_pred_x) <- dim(fit_pred$sample$anc_alpha_obs_t1)


#' ## Calibrated model outputs

outputs_calib <- calibrate_outputs(outputs,
                                   naomi_mf,
                                   spectrum_plhiv_calibration_level = "national",
                                   spectrum_plhiv_calibration_strat = "sex_age_coarse",
                                   spectrum_artnum_calibration_level = "national", 
                                   spectrum_artnum_calibration_strat = "sex_age_coarse",
                                   spectrum_aware_calibration_level = "national", 
                                   spectrum_aware_calibration_strat = "sex_age_coarse",
                                   spectrum_infections_calibration_level = "national", 
                                   spectrum_infections_calibration_strat = "sex_age_coarse",
                                   calibrate_method = "logistic")


indicators_calib <- add_output_labels(outputs_calib) %>%
  left_join(
    select(outputs_calib$meta_area, area_level, area_id, center_x, center_y),
    by = c("area_level", "area_id")
  ) %>%
  sf::st_as_sf() %>%
  mutate(
    area_level_label = fct_reorder(area_level_label, area_level),
    age_group_label = fct_reorder(age_group_label, as.integer(factor(age_group)))
  )


#' ## Fit model with no ART attendance

naomi_mf0 <- naomi_model_frame(areas,
                               pop,
                               spec,
                               scope = scope,
                               level = level,
                               calendar_quarter_t1,
                               calendar_quarter_t2,
                               calendar_quarter_t3,
                               artattend = FALSE,
                               artattend_t2 = FALSE)

naomi_data0 <- select_naomi_data(naomi_mf0,
                                 survey,
                                 anc,
                                 art,
                                 prev_survey_ids,
                                 artcov_survey_ids,
                                 recent_survey_ids,
                                 vls_survey_ids,
                                 artnum_calendar_quarter_t1,
                                 artnum_calendar_quarter_t2,
                                 anc_clients_year_t2 = anc_clients_year2,
                                 anc_clients_year_t2_num_months = anc_clients_year2_num_months,
                                 anc_prevalence_year1,
                                 anc_prevalence_year2,
                                 anc_art_coverage_year1,
                                 anc_art_coverage_year2)

tmb_inputs0 <- prepare_tmb_inputs(naomi_data0)

fit0 <- fit_tmb(tmb_inputs0)
fit0 <- sample_tmb(fit0)
outputs0 <- output_package(fit0, naomi_mf0)

indicators0 <- add_output_labels(outputs0) %>%
  mutate(
    area_level_label = fct_reorder(area_level_label, area_level),
    age_group_label = fct_reorder(age_group_label, as.integer(factor(age_group)))
  )



#' ## Figure 3
#'
#' Summary of model results at district level
#'
#' * HIV prevalence 15-49 at each level
#' * HIV prevalence, ART coverage, incidence rate
#' * Bubble plot PLHIV
#'
#'

## Figure 3A (180 x 85 mm)

fig3adat <- indicators %>%
  filter(
    indicator == "prevalence",
    age_group == "Y015_049",
    calendar_quarter == "CY2018Q3",
    sex == "both"
  )

fig3a <- fig3adat %>%
  ggplot(aes(fill = mean)) +
  geom_sf(size = 0.1) +
  facet_wrap(~area_level_label, nrow = 1) +
  scale_fill_viridis_c(option = "C", direction = -1,
                       begin = 0.1, end = 0.9,
                       labels = label_percent(1)) +
  expand_limits(fill = 0) +
  labs(tag = "A",
       fill = "HIV prevalence\n15-49 years\nSept 2018") +
  coord_sf(expand = FALSE) +  
  theme_minimal(10) +
  theme(
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines"),
    axis.text = element_blank(),
    strip.text = element_text(face = "bold", size = rel(0.85)),
    legend.title = element_text(size = rel(0.7)),
    legend.text = element_text(size = rel(0.7)),
    plot.tag = element_text(face = "bold")
  )

fig3bdat <- indicators %>%
  filter(
    indicator == "art_coverage",
    age_group == "Y015_999",
    calendar_quarter %in% c("CY2016Q1", "CY2018Q3"),
    sex == "both",
    area_level == 4
  )

fig3b <- fig3bdat %>%
  ggplot(aes(fill = mean)) +
  geom_sf(size = 0.1) +
  facet_wrap(~quarter_label, nrow = 1) +
  scale_fill_viridis_c(option = "D", direction = -1,
                       begin = 0.05, end = 0.9,
                       labels = label_percent(1),
                       limits = c(0.6, 0.853)) +
  labs(tag = "B",
       fill = "ART coverage\n15+ years") +
  coord_sf(expand = FALSE) +  
  theme_minimal(10) +
  theme(
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines"),    
    axis.text = element_blank(),
    strip.text = element_text(face = "bold", size = rel(0.85)),
    legend.title = element_text(size = rel(0.7)),
    legend.text = element_text(size = rel(0.7)),
    plot.tag = element_text(face = "bold")
  )

fig3cdat <- indicators %>%
  filter(
    indicator == "incidence" & age_group == "Y015_049" |
    indicator == "infections" & age_group == "Y015_999",
    calendar_quarter == "CY2018Q3",
    sex == "both",
    area_level == 4
  )

fig3c <- fig3cdat %>%
  st_drop_geometry() %>%
  pivot_wider(c(area_id, area_name), names_from = indicator, values_from = mean) %>%
  left_join(select(areas, area_id, center_x, center_y),
            by = "area_id") %>%
  st_as_sf() %>%
  ggplot(aes(x = center_x, y = center_y, color = incidence, size = infections)) +
  geom_sf(size = 0.1, color = "grey30") +
  geom_point(alpha = 0.7) +
  scale_color_viridis_c(option = "B", direction = -1,
                        begin = 0.05, end = 0.9,
                        labels = label_number(scale = 1000)) +
  scale_size_area(max_size = 8) +
  labs(tag = "C",
       color = "Incidence per 1000\n15-49 years\nSept 2018",
       size = "Annual new infection\n15+ years\nSept 2018",
       x = element_blank(),
       y = element_blank()) +
  expand_limits(color = 0) +
  coord_sf(expand = FALSE) +    
  theme_minimal(10) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    strip.text = element_text(face = "bold", size = rel(0.85)),
    legend.title = element_text(size = rel(0.7)),
    legend.text = element_text(size = rel(0.7)),
    plot.tag = element_text(face = "bold")
  )


## quartz(width = 180 / 25.4, height = 180 / 25.4)

fig3 <- grid.arrange(fig3a, fig3b, fig3c,
                     layout_matrix = rbind(1, c(NA, 2, 3, NA)),
                     widths = c(0.1, 1, 0.8, 0.1))

ggsave("figure3.pdf", fig3, width = 180, height = 180, units = "mm")
ggsave("figure3.png", fig3, width = 180, height = 180, units = "mm")


#' ## Figure 4
#'
#' Sex/age group pyramid for indicators:
#' * Population
#' * PLHIV
#' * New HIV infections
#' * Number on ART
#' * Number untreated
#' * Number undiagnosed
#'

fig4dat <- indicators %>%
  filter(
    indicator %in% c("population", "plhiv", "art_current_residents",
                     "unaware_plhiv_num", "untreated_plhiv_num", "infections"),
    calendar_quarter == "CY2018Q3",
    sex %in% c("male", "female"),
    area_level == 0
  ) %>%
  semi_join(
    filter(naomi::get_age_groups(), age_group_span == 5 | age_group == "Y080_999"),
    by = "age_group"
  ) %>%
  mutate(
    sign = recode(sex, "male" = -1, "female" = 1),
    mean = mean * sign,
    lower = lower * sign,
    upper = upper * sign,
    indicator_label = indicator_label %>%
      recode("PLHIV" = "People living with HIV",
             "ART number (residents)" = "PLHIV on ART",
             "Number PLHIV unaware" = "PLHIV unaware of HIV status") %>%
      factor(c("Population", "People living with HIV", "New infections",
               "PLHIV on ART", "PLHIV not on ART", "PLHIV unaware of HIV status"))
  )

fig4lims <- c("Population" = 1.5e6,
              "People living with HIV" = 120e3,
              "New infections" = 25e3,
              "PLHIV on ART" = 82e3,
              "PLHIV not on ART" = 40e3,
              "PLHIV unaware of HIV status" = 40e3) %>%
  data.frame(indicator_label = names(.), max = .) %>%
  mutate(indicator_label = fct_inorder(indicator_label))

fig4 <- fig4dat %>%
  ggplot(aes(age_group_label, mean, ymin = lower, ymax = upper,
             fill = sex, color = sex)) +
  geom_blank(aes(yintercept = max), data = fig4lims, inherit.aes = FALSE) +
  geom_blank(aes(yintercept = -max), data = fig4lims, inherit.aes = FALSE) +
  geom_col(alpha = 0.7, color = NA) +
  geom_linerange() +
  annotate("text", x = 17, y = -Inf, label = "Male", fontface = "bold",
           size = 2.5, hjust = -0.2, vjust = 1) +
  annotate("text", x = 17, y = Inf, label = "Female", fontface = "bold",
           size = 2.5, hjust = 1.2, vjust = 1) +
  scale_y_continuous(labels = function(x) number(abs(x), scale = 1e-3)) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  coord_flip() +
  labs(y = "Number in thousands", x = element_blank()) +
  facet_wrap(~indicator_label, scales = "free_x", ncol = 3) +
  theme_bw(10) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = rel(1.0))
  )


ggsave("figure4.pdf", fig4, width = 180, height = 130, units = "mm")
ggsave("figure4.png", fig4, width = 180, height = 130, units = "mm")


#' ## Figure 5

#' Posterior predictive distribution
#'
#' * Survey prevalence 15-49
#' * Survey ART coverage 15-64
#' * ANC prevalence x2
#' * ANC ART coverage x2
#' * ART number x2

prev_pred <- naomi_data$prev_dat %>%
  mutate(
    nsim = round(naomi_data$prev_dat$n_eff)
  ) %>%
  bind_cols(as.data.frame(prev_pred_x))

prev_pred <- prev_pred %>%
  semi_join(
    get_age_groups() %>%
    filter(age_group_span == 5,
           age_group_start >= 15,
           age_group_start + age_group_span <= 50)
  ) %>%
  group_by(area_id, survey_id) %>%
  summarise(
    nsim = sum(nsim),
    estsim = sum(x_eff) / sum(n_eff),
    across(V1:V1000, sum),
    .groups = "drop"
  ) %>%
  mutate(
    across(V1:V1000, `/`, nsim)
  ) %>%
  pivot_longer(V1:V1000, names_to = "sample") %>%
  group_by(area_id, survey_id, estsim) %>%
  summarise(
    pp_mean = mean(value),
    pp_median = quantile(value, 0.5, names = FALSE),
    pp_lower = quantile(value, 0.025, names = FALSE),
    pp_upper = quantile(value, 0.975, names = FALSE),
    pp_lower80 = quantile(value, 0.1, names = FALSE),
    pp_upper80 = quantile(value, 0.9, names = FALSE),
    .groups = "drop"
  ) %>%
  full_join(
    indicators %>%
    filter(
      indicator == "prevalence",
      area_level == 4,
      calendar_quarter == "CY2016Q1",
      sex == "both",
      age_group == "Y015_049"
    )
  ) %>%
  left_join(
    survey %>%
    filter(indicator == "prevalence",
           sex == "both",
           age_group == "Y015_049") %>%
    select(indicator, survey_id, area_id, sex, age_group, estimate)
  ) %>%
  mutate(
    farea_name = fct_reorder(area_name, mean, .desc = TRUE)
  )


artcov_pred <- naomi_data$artcov_dat %>%
  mutate(
    nsim = round(naomi_data$artcov_dat$n_eff)
  ) %>%
  bind_cols(as.data.frame(artcov_pred_x))

artcov_pred <- artcov_pred %>%
  semi_join(
    get_age_groups() %>%
    filter(age_group_span == 5,
           age_group_start >= 15,
           age_group_start + age_group_span <= 65)
  ) %>%
  group_by(area_id, survey_id) %>%
  summarise(
    nsim = sum(nsim),
    across(V1:V1000, sum),
    .groups = "drop"
  ) %>%
  mutate(
    across(V1:V1000, `/`, nsim)
  ) %>%
  pivot_longer(V1:V1000, names_to = "sample") %>%
  group_by(area_id, survey_id) %>%
  summarise(
    pp_mean = mean(value),
    pp_median = quantile(value, 0.5, names = FALSE),
    pp_lower = quantile(value, 0.025, names = FALSE),
    pp_upper = quantile(value, 0.975, names = FALSE),
    pp_lower80 = quantile(value, 0.1, names = FALSE),
    pp_upper80 = quantile(value, 0.9, names = FALSE),
    .groups = "drop"
  ) %>%
  full_join(
    indicators %>%
    filter(
      indicator == "art_coverage",
      area_level == 4,
      calendar_quarter == "CY2016Q1",
      sex == "both",
      age_group == "Y015_064"
    )
  ) %>%
  left_join(
    survey %>%
    filter(indicator == "art_coverage",
           sex == "both",
           age_group == "Y015_064") %>%
    select(indicator, survey_id, area_id, sex, age_group, estimate)
  ) %>%
  mutate(
    farea_name = fct_reorder(area_name, mean, .desc = TRUE)
  )

ancprev_pred <- naomi_data$anc_prev_t1_dat %>%
  mutate(
    estimate = anc_prev_x / anc_prev_n,
  ) %>%
  bind_cols(as.data.frame(ancprev_pred_x))

ancprev_pred <- ancprev_pred %>%
  mutate(
    across(V1:V1000, `/`, anc_prev_n)
  ) %>%
  pivot_longer(V1:V1000, names_to = "sample") %>%
  group_by(area_id, anc_prev_x, anc_prev_n, estimate) %>%
  summarise(
    pp_mean = mean(value),
    pp_median = quantile(value, 0.5, names = FALSE),
    pp_lower = quantile(value, 0.025, names = FALSE),
    pp_upper = quantile(value, 0.975, names = FALSE),
    pp_lower80 = quantile(value, 0.1, names = FALSE),
    pp_upper80 = quantile(value, 0.9, names = FALSE),
    .groups = "drop"
  ) %>%
  full_join(
    indicators %>%
    filter(
      indicator == "anc_prevalence",
      area_level == 4,
      calendar_quarter == "CY2016Q1",
      sex == "female",
      age_group == "Y015_049"
    )
  ) %>%
  mutate(
    farea_name = factor(area_name, levels(prev_pred$farea_name))
  )

ancartcov_pred <- naomi_data$anc_artcov_t1_dat %>%
  mutate(
    estimate = anc_artcov_x / anc_artcov_n,
  ) %>%
  bind_cols(as.data.frame(ancartcov_pred_x))

ancartcov_pred <- ancartcov_pred %>%
  mutate(
    across(V1:V1000, `/`, anc_artcov_n)
  ) %>%
  pivot_longer(V1:V1000, names_to = "sample") %>%
  group_by(area_id, anc_artcov_x, anc_artcov_n, estimate) %>%
  summarise(
    pp_mean = mean(value),
    pp_median = quantile(value, 0.5, names = FALSE),
    pp_lower = quantile(value, 0.025, names = FALSE),
    pp_upper = quantile(value, 0.975, names = FALSE),
    pp_lower80 = quantile(value, 0.1, names = FALSE),
    pp_upper80 = quantile(value, 0.9, names = FALSE),
    .groups = "drop"
  ) %>%
  full_join(
    indicators %>%
    filter(
      indicator == "anc_art_coverage",
      area_level == 4,
      calendar_quarter == "CY2016Q1",
      sex == "female",
      age_group == "Y015_049"
    )
  ) %>%
  mutate(
    farea_name = factor(area_name, levels(artcov_pred$farea_name))
  )


fig5a <- prev_pred %>%
  ggplot(aes(x = farea_name, y = mean)) +
  geom_linerange(aes(ymin = lower, ymax = upper), size = 2, color = "grey25") +
  geom_crossbar(ymin = NA, ymax = NA, width = 0.6, size = 0.4, color = "grey25") +
  geom_linerange(aes(ymin = pp_lower80, ymax = pp_upper80, group = survey_id), color = "lightblue3", position = position_dodge(width = 0.35)) +
  geom_point(aes(y = estimate, shape = survey_id), color = "red", position = position_dodge(width = 0.35)) +
  theme_bw(10) +
  labs(tag = "A", y = element_blank(), x = element_blank(), shape = element_blank()) +
  scale_shape_manual(
    element_blank(),
    values = c("DEMO2015DHSA" = 17, "DEMO2016PHIA" = 15),
    labels = c("DEMO2015DHSA" = "MDHS 2015-16", "DEMO2016PHIA" = "MPHIA 2015-16*")
  ) +
  scale_y_continuous("HIV prevalence, 15-49y", labels = label_percent(1)) +
  expand_limits(y = 0) +
  theme(
    legend.position = c(1, 1),
    legend.just = c(1, 0.85),
    legend.background = element_rect(fill = NA),
    legend.key.size = unit(1, "lines"),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    plot.tag = element_text(face = "bold"),
    plot.tag.position = c(0, 1.0),    
    plot.background = element_rect(fill = NA),
    plot.margin = margin(12, 2, -3, 10, "pt")
  )

fig5b <- ancprev_pred %>%
  ggplot(aes(x = farea_name, y = mean)) +
  geom_point(size = 1.5, color = "grey25", alpha = 0.2,
             data = filter(prev_pred, survey_id == "DEMO2015DHSA")) +
  geom_linerange(aes(ymin = lower, ymax = upper), size = 2, color = "grey25") +
  geom_crossbar(ymin = NA, ymax = NA, width = 0.6, size = 0.4, color = "grey25") +
  geom_linerange(aes(ymin = pp_lower80, ymax = pp_upper80), color = "lightblue3", position = position_dodge(width = 0.35)) +
  geom_point(aes(y = estimate), color = "red", position = position_dodge(width = 0.35)) +
  theme_bw(10) +
  labs(tag = "B", y = element_blank(), x = element_blank(), shape = element_blank()) +
  scale_y_continuous("ANC HIV prevalence", labels = label_percent(1)) +
  expand_limits(y = 0) +
  theme(
    legend.position = c(1, 1),
    legend.just = c(1, 0.85),
    legend.background = element_rect(fill = NA),
    legend.key.size = unit(1, "lines"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.tag = element_text(face = "bold"),
    plot.tag.position = c(0, 1.0),
    plot.background = element_rect(fill = NA),
    plot.margin = margin(0, 2, 0, 10, "pt")
  )

fig5c <- artcov_pred %>%
  ggplot(aes(x = farea_name, y = mean)) +
  geom_linerange(aes(ymin = lower, ymax = upper), size = 2, color = "grey25") +
  geom_crossbar(ymin = NA, ymax = NA, width = 0.6, size = 0.4, color = "grey25") +
  geom_linerange(aes(ymin = pp_lower80, ymax = pp_upper80, group = survey_id), color = "lightblue3", position = position_dodge(width = 0.35)) +
  geom_point(aes(y = estimate, shape = survey_id), color = "red", position = position_dodge(width = 0.35)) +
  theme_bw(10) +
  labs(tag = "C", y = element_blank(), x = element_blank(), shape = element_blank()) +
  scale_shape_manual(
    element_blank(),
    values = c("DEMO2016PHIA" = 15),
    labels = c("DEMO2016PHIA" = "MPHIA 2015-16*")
  ) +
  scale_y_continuous("ART coverage, 15-64y", labels = label_percent(1), limits = c(0, 1)) +
  theme(
    legend.position = c(1, 1),
    legend.just = c(1, 0.75),
    legend.background = element_rect(fill = NA),
    legend.key.size = unit(1, "lines"),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    plot.tag = element_text(face = "bold"),
    plot.tag.position = c(0, 1.0),    
    plot.background = element_rect(fill = NA),
    plot.margin = margin(12, 2, -3, 10, "pt")
  )

fig5d <- ancartcov_pred %>%
  ggplot(aes(x = farea_name, y = mean)) +
  geom_point(size = 1.5, color = "grey25", alpha = 0.2,
             data = artcov_pred) +
  geom_linerange(aes(ymin = lower, ymax = upper), size = 2, color = "grey25") +
  geom_crossbar(ymin = NA, ymax = NA, width = 0.6, size = 0.4, color = "grey25") +
  geom_linerange(aes(ymin = pp_lower80, ymax = pp_upper80), color = "lightblue3", position = position_dodge(width = 0.35)) +
  geom_point(aes(y = estimate), color = "red", position = position_dodge(width = 0.35)) +
  theme_bw(10) +
  labs(tag = "D", y = element_blank(), x = element_blank(), shape = element_blank()) +
  scale_y_continuous("ANC ART coverage", labels = label_percent(1), limits = c(0, 1)) +
  theme(
    legend.position = c(1, 1),
    legend.just = c(1, 0.85),
    legend.background = element_rect(fill = NA),
    legend.key.size = unit(1, "lines"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.tag = element_text(face = "bold"),
    plot.tag.position = c(0, 1.0),
    plot.background = element_rect(fill = NA),
    plot.margin = margin(0, 2, 0, 10, "pt")
  )

fig5 <- grid.arrange(fig5a, fig5b, fig5c, fig5d, heights = c(1, 1.15, 1, 1.15))

ggsave("figure5.pdf", fig5, width = 180, height = 220, units = "mm")
ggsave("figure5.png", fig5, width = 180, height = 220, units = "mm")


#' ## Figure 6
#'
#' Summary of results about ART attendance for Lilongwe City (MWI_4_14_demo),
#' Lilongwe (Rural) (MWI_4_15_demo), and Dowa (MWI_4_11_demo).
#'
#' * A: Number of residents on ART and ART clients at facilities
#' * B: Distribution of where residents seek treatment
#' * E: Distribution of where ART clients reside

fig6areas <- c("MWI_4_14_demo", "MWI_4_15_demo", "MWI_4_11_demo")

area_order <- c("Lilongwe City", "Lilongwe", "Dowa", "Salima", "Dedza", "Mchinji", "Kasungu", "Ntchisi")

fig6adat <- indicators %>%
  filter(
    area_name %in% area_order,
    indicator %in% c("art_current", "art_current_residents"),
    age_group == "Y000_999",
    sex == "both",
    calendar_quarter == "CY2018Q3",
    area_level == 4
  ) %>%
  mutate(
    area_name = area_name %>%
      fct_relevel(!!area_order) %>%
      fct_recode("Lilongwe\n(Rural)" = "Lilongwe",
                 "Lilongwe\nCity" = "Lilongwe City"),
    indicator = fct_relevel(indicator, "art_current_residents", "art_current")
  )

fig6a <- fig6adat %>%
  ggplot(aes(area_name, mean, ymin = lower, ymax = upper, fill = indicator, color = indicator)) +
  geom_col(position = "dodge", alpha = 0.7, color = NA) +
  geom_linerange(position = position_dodge(width = 0.9), size = 1.0) +
  scale_color_brewer(palette = "Set1", guide = "none") +
  scale_fill_brewer(palette = "Set1",
                    labels = c("art_current_residents" = "Residents on ART",
                               "art_current" = "ART clients at facilities")) +  
  labs(tag = "A", x = element_blank(), y = element_blank(),
       color = element_blank(), fill = element_blank()) +
  theme_bw(10) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = rel(0.85)),
    plot.tag = element_text(face = "bold"),
    legend.position = c(1, 1.05),
    legend.background = element_blank(),
    legend.just = c(1, 1)
  )


fig6bcdat <- outputs$art_attendance %>%
  left_join(
    areas %>%
    st_drop_geometry() %>%
    select(area_id, area_name) %>% rename_all(~paste0("reside_", .)),
    by = "reside_area_id"
  ) %>%
    left_join(
    areas %>%
    st_drop_geometry() %>%
    select(area_id, area_name) %>% rename_all(~paste0("attend_", .)),
    by = "attend_area_id"
    ) %>%
  filter(calendar_quarter == "CY2018Q3") %>%
  mutate(
    reside_area_name = fct_relevel(reside_area_name, !!area_order) %>%
      fct_recode("Lilongwe (Rural)" = "Lilongwe"),
    attend_area_name = fct_relevel(attend_area_name, !!area_order) %>%
      fct_recode("Lilongwe (Rural)" = "Lilongwe")
  )

  
fig6b <- fig6bcdat %>%
  filter(
    reside_area_name %in% c("Lilongwe City", "Lilongwe (Rural)", "Dowa")
  ) %>%
  arrange(
    reside_area_name,
    attend_area_name != reside_area_name,
    attend_area_name
  ) %>%
  mutate(
    reside_area_name = paste0("Reside: ", reside_area_name) %>%
      fct_inorder(),
    attend_area_name = fct_inorder(paste(attend_area_name, row_number()))
  ) %>%
  ggplot(aes(attend_area_name, prop_residents_mean,
             ymin = prop_residents_lower, ymax = prop_residents_upper)) +
  geom_col(fill = "grey45") +
  geom_linerange(size = 1.0) +
  facet_grid(~reside_area_name, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = label_percent()) +
  scale_x_discrete(labels = function(x) sub(" [0-9]+$", "", x)) +
  labs(tag = "B", x = "District attending ART", y = element_blank()) +
  theme_bw(10) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = rel(0.8)),
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.tag = element_text(face = "bold")
  )

fig6c <- fig6bcdat %>%
  filter(
    attend_area_name %in% c("Lilongwe City", "Lilongwe (Rural)", "Dowa")
  ) %>%
  arrange(
    attend_area_name,
    reside_area_name != attend_area_name,
    reside_area_name
  ) %>%
  mutate(
    attend_area_name = paste0("Attend: ", attend_area_name) %>%
      fct_inorder(),
    reside_area_name = fct_inorder(paste(reside_area_name, row_number()))
  ) %>%
  ggplot(aes(reside_area_name, prop_attendees_mean,
             ymin = prop_attendees_lower, ymax = prop_attendees_upper)) +
  geom_col(fill = "grey45") +
  geom_linerange(size = 1.0) +
  facet_grid(~attend_area_name, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = label_percent()) +
  scale_x_discrete(labels = function(x) sub(" [0-9]+$", "", x)) +
  labs(tag = "C", x = "District of residence", y = element_blank()) +
  theme_bw(10) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = rel(0.8)),
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.tag = element_text(face = "bold")
  )


fig6 <- grid.arrange(fig6a, fig6b, fig6c)

ggsave("figure6.pdf", fig6, width = 180, height = 180, units = "mm")
ggsave("figure6.png", fig6, width = 180, height = 180, units = "mm")

#' ## Supplementary Figure S1
#'
#' Comparison of raw and calibrated results.
#' 
#' S1A: National by sex
#' *  HIV prevalence 15-49
#' *  ART coverage 15+
#' *  Incidence rate 15-49
#'
#' S1B: HIV prevalence by district
#'
#' S1C: ART coverage by district


figS1dat <- indicators %>%
  mutate(version = "Raw") %>%
  bind_rows(
    indicators_calib %>%
    mutate(version = "Calibrated")
  ) %>%
  filter(
    calendar_quarter == "CY2018Q3",
    age_group == "Y015_049" & indicator %in% c("prevalence", "incidence") |
    age_group == "Y015_999" & indicator %in% "art_coverage"
  ) %>%
  mutate(
    version = factor(version, c("Raw", "Calibrated")),
    sex = factor(sex, c("both", "female", "male"), c("Both", "Female", "Male"))
  )


figS1ai <- figS1dat %>%
  filter(
    area_id == "MWI",
    indicator == "prevalence"
  ) %>%
  ggplot(aes(sex, mean, ymin = lower, ymax = upper, fill = version)) +
  geom_col(position = position_dodge(0.9)) +
  geom_linerange(position = position_dodge(0.9)) +
  geom_text(aes(y = upper + 0.006, label = percent(mean, 0.1)),
            position = position_dodge(0.9), fontface = "bold", size = 2.3) +
  scale_y_continuous(label = label_percent(1.0), expand = expansion(c(0, 0.05)),
                     limits = c(0, 0.15)) +
  scale_fill_brewer(palette = "Set1") +
  labs(tag = "A", title = "HIV prevalence, 15-49y", x = NULL, y = NULL, fill = NULL) +
  theme_classic(10) +
  theme(    
    legend.position = "none",
    plot.title = element_text(size = rel(0.9), hjust = 0.5, face = "bold"),
    plot.tag.position = c(0, 1.0),
    plot.tag = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )

figS1aii <- figS1dat %>%
  filter(
    area_id == "MWI",
    indicator == "art_coverage"
  ) %>%
  ggplot(aes(sex, mean, ymin = lower, ymax = upper, fill = version)) +
  geom_col(position = position_dodge(0.9)) +
  geom_linerange(position = position_dodge(0.9)) +
  geom_text(aes(y = upper + 0.04, label = percent(mean, 1.0)),
            position = position_dodge(0.9), fontface = "bold", size = 2.3) +
  scale_y_continuous(label = label_percent(), expand = expansion(c(0, 0.05)),
                     limits = c(0, 0.95)) +
  scale_fill_brewer(palette = "Set1") +
  labs(tag = "", title = "ART coverage, 15+y", x = NULL, y = NULL, fill = NULL) +
  theme_classic(10) +
  theme(    
    legend.position = "none",
    plot.title = element_text(size = rel(0.9), hjust = 0.5, face = "bold"),
    plot.tag.position = c(0, 1.0),
    plot.tag = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )

figS1aiii <- figS1dat %>%
  filter(
    area_id == "MWI",
    indicator == "incidence"
  ) %>%
  ggplot(aes(sex, mean, ymin = lower, ymax = upper, fill = version)) +
  geom_col(position = position_dodge(0.9)) +
  geom_linerange(position = position_dodge(0.9)) +
  geom_text(aes(y = upper + 0.001, label = number(mean, 0.1, scale = 1e3)),
            position = position_dodge(0.9), fontface = "bold", size = 2.3) +
  scale_y_continuous(label = label_number(scale = 1e3), expand = expansion(c(0, 0.05))) +
  scale_fill_brewer(palette = "Set1") +
  labs(tag = "", title = "Incidence per 1000, 15-49y", x = NULL, y = NULL, fill = NULL) +
  theme_classic(10) +
  theme(    
    legend.position = "none",
    plot.title = element_text(size = rel(0.9), hjust = 0.5, face = "bold"),
    plot.tag.position = c(0, 1.0),
    plot.tag = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )


figS1b <- figS1dat %>%
  filter(
    area_level == 4,
    indicator == "prevalence",
    sex == "Both"
  ) %>%
  mutate(
    farea_name = fct_reorder(area_name, mean, .desc = TRUE)
  ) %>%
  ggplot(aes(farea_name, mean, ymin = lower, ymax = upper, fill = version)) +
  geom_col(position = position_dodge(0.9)) +
  geom_linerange(position = position_dodge(0.9)) +
  scale_y_continuous(label = label_percent(1.0), expand = expansion(c(0, 0.05))) +
  scale_fill_brewer(palette = "Set1") +
  labs(tag = "B", y = "HIV prevalence, 15-49y", x = NULL, fill = NULL) +
  theme_classic(10) +
  theme(    
    legend.position = "none",
    plot.title = element_text(size = rel(0.9), hjust = 0.5, face = "bold"),
    plot.tag.position = c(0, 1.08),
    plot.tag = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1.0),
    plot.margin = margin(20, 5.5, 5.5, 5.5, "pt")
  )

figS1c <- figS1dat %>%
  filter(
    area_level == 4,
    indicator == "art_coverage",
    sex == "Both"
  ) %>%
  mutate(
    farea_name = fct_reorder(area_name, mean, .desc = TRUE)
  ) %>%
  ggplot(aes(farea_name, mean, ymin = lower, ymax = upper, fill = version)) +
  geom_col(position = position_dodge(0.9)) +
  geom_linerange(position = position_dodge(0.9)) +
  scale_y_continuous(label = label_percent(1.0), expand = expansion(c(0, 0.05))) +
  scale_fill_brewer(palette = "Set1") +
  labs(tag = "C", y = "ART coverage, 15+y", x = NULL, fill = NULL) +
  coord_cartesian(ylim = c(0.6, 0.9)) +
  theme_classic(10) +
  theme(    
    legend.position = "none",
    plot.title = element_text(size = rel(0.9), hjust = 0.5, face = "bold"),
    plot.tag.position = c(0, 1.08),
    plot.tag = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1.0),
    plot.margin = margin(20, 5.5, 5.5, 5.5, "pt")
  )


figS1leg <- cowplot::get_legend(figS1ai + theme(legend.position = "bottom"))
  
figS1a <- grid.arrange(figS1ai, figS1aii, figS1aiii, nrow = 1)
figS1 <- grid.arrange(figS1a, figS1b, figS1c, figS1leg, ncol = 1, heights = c(1, 1, 1, 0.1))

ggsave("figureS1.pdf", figS1, width = 180, height = 190, units = "mm")
ggsave("figureS1.png", figS1, width = 180, height = 190, units = "mm")



#' ## Supplementary Figure S2
#'
#' Comparison of results with ART attendance = FALSE
#'
#' * HIV prevalence 15-49
#' * ART coverage 15+
#' * Incidence rate
#' * PLHIV 15+
#' * Residents on ART
#' * Attending ART
#' * Unmet need for ART

figS2dat <- indicators %>%
  st_drop_geometry() %>%
  mutate(version = "full") %>%
  bind_rows(
    mutate(indicators0, version = "no_artattend")
  ) %>%
  filter(
    area_id %in% fig6areas,
    calendar_quarter == "CY2018Q3",
    sex == "both",
    indicator %in% c("prevalence", "incidence") & age_group == "Y015_049" |
    indicator %in% c("plhiv", "art_coverage", "plhiv",
                     "art_current", "art_current_residents", "untreated_plhiv_num") &
    age_group == "Y015_999"
  ) %>%
  mutate(
    area_name = area_name %>%
      fct_relevel("Lilongwe City", "Lilongwe", "Dowa") %>%
      fct_recode("Lilongwe\nCity" = "Lilongwe City",
                 "Lilongwe\n(Rural)" = "Lilongwe")
  )


figS2a <- figS2dat %>%
  filter(indicator == "plhiv") %>%
  ggplot(aes(area_name, mean, ymin = lower, ymax = upper, fill = version, color = version)) +
  geom_col(position = "dodge", color = NA, alpha = 1.0) +
  geom_linerange(position = position_dodge(0.9), size = 1.0) +
  labs(tag = "A",
       title = "PLHIV (15+; 000s)",
       x = element_blank(),
       y = element_blank(),
       fill = element_blank()) +
  scale_fill_brewer(palette = "Set2",
                    labels = c("full" = "With ART attend",
                               "no_artattend" = "No ART attend")) +
  scale_color_brewer(palette = "Dark2", guide = FALSE) +
  scale_y_continuous(labels = label_number(scale=1e-3)) +
  expand_limits(y = 1e5) +  
  theme_bw(10) + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = rel(0.9)),
    plot.tag = element_text(face = "bold"),
    legend.position = c(1, 1.05),
    legend.just = c(1, 1),
    legend.background = element_blank(),
    legend.text = element_text(size = rel(0.7)),
    legend.key.size = unit(1, "lines")
  )

figS2b <- figS2dat %>%
  filter(indicator == "art_current") %>%
  ggplot(aes(area_name, mean, ymin = lower, ymax = upper, fill = version, color = version)) +
  geom_col(position = "dodge", color = NA, alpha = 1.0) +
  geom_linerange(position = position_dodge(0.9), size = 1.0) +
  labs(tag = "B",
       title = "ART attending (15+; 000s)",
       x = element_blank(),
       y = element_blank()) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(labels = label_number(scale=1e-3)) +
  expand_limits(y = 1e5) +
  theme_bw(10) + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = rel(0.9)),
    plot.tag = element_text(face = "bold"),
    legend.position = "none"
  )

figS2c <- figS2dat %>%
  filter(indicator == "art_current_residents") %>%
  ggplot(aes(area_name, mean, ymin = lower, ymax = upper, fill = version, color = version)) +
  geom_col(position = "dodge", color = NA, alpha = 1.0) +
  geom_linerange(position = position_dodge(0.9), size = 1.0) +
  labs(tag = "C",
       title = "Residents on ART (15+; 000s)",
       x = element_blank(),
       y = element_blank()) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(labels = label_number(scale=1e-3)) +
  expand_limits(y = 1e5) +
  theme_bw(10) + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = rel(0.9)),
    plot.tag = element_text(face = "bold"),
    legend.position = "none"
  )



figS2d <- figS2dat %>%
  filter(indicator == "prevalence") %>%
  ggplot(aes(area_name, mean, ymin = lower, ymax = upper, fill = version, color = version)) +
  geom_col(position = "dodge", color = NA, alpha = 1.0) +
  geom_linerange(position = position_dodge(0.9), size = 1.0) +
  labs(tag = "D",
       title = "HIV prevalence (15-49)",
       x = element_blank(),
       y = element_blank()) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(labels = label_percent(1)) +
  theme_bw(10) + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = rel(0.9)),
    plot.tag = element_text(face = "bold"),
    legend.position = "none"
  )

figS2e <- figS2dat %>%
  filter(indicator == "art_coverage") %>%
  ggplot(aes(area_name, mean, ymin = lower, ymax = upper, fill = version, color = version)) +
  geom_col(position = "dodge", color = NA, alpha = 1.0) +
  geom_linerange(position = position_dodge(0.9), size = 1.0) +
  labs(tag = "E",
       title = "ART coverage (15+)",
       x = element_blank(),
       y = element_blank()) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(labels = label_percent(1)) +
  expand_limits(y = 1.0) +
  theme_bw(10) + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = rel(0.9)),
    plot.tag = element_text(face = "bold"),
    legend.position = "none"
  )

figS2f <- figS2dat %>%
  filter(indicator == "incidence") %>%
  ggplot(aes(area_name, mean, ymin = lower, ymax = upper, fill = version, color = version)) +
  geom_col(position = "dodge", color = NA, alpha = 1.0) +
  geom_linerange(position = position_dodge(0.9), size = 1.0) +
  labs(tag = "F",
       title = "HIV incidence (15-49; per 1000)",
       x = element_blank(),
       y = element_blank(),
       fill = element_blank()) +
  scale_fill_brewer(palette = "Set2",
                    labels = c("full" = "With ART attendance",
                               "no_artattend" = "No ART attendance"),
                    guide = FALSE) +
  scale_color_brewer(palette = "Dark2", guide = FALSE) +
  scale_y_continuous(labels = label_number(scale = 1e3)) +
  theme_bw(10) + 
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = rel(0.9)),
    plot.tag= element_text(face = "bold"),
    legend.position = c(1, 1.05),
    legend.just = c(1, 1),
    legend.background = element_blank(),
    legend.text = element_text(size = rel(0.7)),
    legend.key.size = unit(1, "lines")
  )


figS2 <- grid.arrange(figS2a, figS2b, figS2c, figS2d, figS2e, figS2f,
                     nrow = 2)


ggsave("figureS2.pdf", figS2, height = 130, width = 180, units = "mm")
ggsave("figureS2.png", figS2, height = 130, width = 180, units = "mm")



#' ## For text

#' Figure 3

fig3adat %>%
  st_drop_geometry() %>%
  select(-contains("label")) %>%
  filter(area_level == 0)

fig3adat %>%
  st_drop_geometry() %>%
  select(-contains("label"), -contains("center")) %>%
  filter(area_level == 4) %>%
  arrange(mean) %>%
  print(n = Inf)

fig3adat %>%
  st_drop_geometry() %>%
  select(-contains("label"), -contains("center")) %>%
  filter(area_level == 4) %>%
  mutate(rse = se / mean) %>%  
  arrange(rse) %>%
  print(n = Inf)

indicators %>%
  filter(
    indicator == "art_coverage",
    age_group == "Y015_999",
    calendar_quarter %in% c("CY2016Q1", "CY2018Q3"),
    sex == "both",
    area_level == 0
  ) %>%
  st_drop_geometry() %>%
  select(-contains("label"), -contains("center"))
  
fig3bdat %>%
  st_drop_geometry() %>%
  select(-contains("label"), -contains("center")) %>%
  arrange(calendar_quarter, mean) %>%
  print(n = Inf)

fig3bdat %>%
  st_drop_geometry() %>%
  pivot_wider(c(area_id, area_name), names_from = calendar_quarter, values_from = mean) %>%
  mutate(diff = CY2018Q3 - CY2016Q1) %>%
  arrange(diff) %>%
  print(n = Inf)

fig3cdat %>%
  st_drop_geometry() %>%
  select(-contains("label"), -contains("center")) %>%
  filter(indicator == "incidence") %>%
  mutate(rse = se / mean) %>%
  arrange(rse) %>%
  print(n = Inf)

fig3cdat %>%
  st_drop_geometry() %>%
  select(-contains("label"), -contains("center")) %>%
  filter(indicator == "infections") %>%
  mutate(rse = se / mean) %>%
  arrange(mean) %>%
  print(n = Inf)

#' Figure 4

fig4dat %>%
  st_drop_geometry() %>%
  select(-contains("label"), -contains("center")) %>%
  group_by(indicator) %>%
  summarise(median(abs(se / mean), na.rm = TRUE))

#' Figure 5


bind_rows(
  mutate(prev_pred, indicator = "prevalence"),
  mutate(artcov_pred, indicator = "art_coverage"),
  mutate(ancprev_pred, indicator = "anc_prevalence"),
  mutate(ancartcov_pred, indicator = "anc_art_coverage")
) %>%
  mutate(indicator = fct_inorder(indicator)) %>%
  filter(!is.na(estimate)) %>%
  group_by(indicator) %>%
  summarise(
    pp80x = sum(estimate >= pp_lower80 & estimate <= pp_upper80),
    pp95x = sum(estimate >= pp_lower & estimate <= pp_upper),
    n = n()
  ) %>%
  mutate(
    pp80cov = pp80x / n,
    pp95cov = pp95x / n
  )
    

#' Figure 6

fig6adat %>%
  st_drop_geometry() %>%
  select(-contains("label"), -contains("center"))

fig6bcdat %>%
  filter(reside_area_name %in% c("Lilongwe City", "Lilongwe (Rural)", "Dowa")) %>%
  select("reside_area_name", "attend_area_name",
         starts_with("prop_residents"), -contains("mode"))

fig6bcdat %>%
  filter(attend_area_name == "Lilongwe City") %>%
  select("reside_area_name", "attend_area_name",
         starts_with("prop_attendees"), -contains("mode"))


#' Figure S2

figS2dat %>%
  select(indicator, area_id, area_name, age_group_label, version, mean, lower, upper) %>%
  arrange(indicator, area_name, version) %>%
  print(n = Inf)

