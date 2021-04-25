#' Analysis conducted using naomi v2.3.14. 
#'
#' Install the Git commit used for analysis with:
#'
#' > devtools::install_github("mrc-ide/naomi@93dc3b1")
#' 

library(naomi)
library(tidyverse)
library(scales)
library(sf)
library(gridExtra)

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
calendar_quarter_t3 <- "CY2019Q4"

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


quartz(width = 180 / 25.4, height = 180 / 25.4)

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
#'
#' * Survey prevalence vs. estimates
#' * Survey ART coverage vs. estimates
#' * ANC prevalence vs. estimates
#' * ANC ART coverage vs. estimates

d <- tmb_inputs$data
p <- fit$par.full %>% split(names(.))

pop15to49 <- naomi_data$mf_model$age15to49 * d$population_t1
pop15plus <- naomi_data$mf_model$age15plus * d$population_t1

mu_asfr <- d$X_asfr %*% p$beta_asfr +
  d$Z_asfr_x %*% p$ui_asfr_x
anc_clients <- as.vector(d$population_t1 * exp(d$log_asfr_t1_offset + mu_asfr))

mu_rho <- d$X_rho %*% p$beta_rho +
  d$logit_rho_offset +
  d$Z_rho_a %*% p$u_rho_a +
  d$Z_rho_as %*% p$u_rho_as

mu_anc_rho <- mu_rho +
  d$logit_anc_rho_t1_offset +
  d$X_ancrho %*% p$beta_anc_rho

mu_alpha <- d$X_alpha %*% p$beta_alpha +
  d$logit_alpha_offset +
  d$Z_alpha_a %*% p$u_alpha_a +
  d$Z_alpha_as %*% p$u_alpha_as

mu_anc_alpha <- mu_alpha +
  d$logit_anc_alpha_t1_offset +
  d$X_ancalpha %*% p$beta_anc_alpha


sum(plogis(as.matrix(mu_rho)) * pop15to49) / sum(pop15to49)
sum(plogis(as.matrix(mu_rho)) * plogis(as.matrix(mu_alpha)) * pop15to49) /
  sum(plogis(as.matrix(mu_rho)) * pop15to49)

ux_seq <- seq(-2.5, 1.5, 0.01)
pred_prev15to49 <- colSums(plogis(outer(mu_rho, ux_seq, "+")) * pop15to49) / sum(pop15to49)
pred_ancprev <- colSums(plogis(outer(mu_anc_rho, ux_seq, "+")) * anc_clients) / sum(anc_clients)

pred_artcov15plus <- colSums(plogis(outer(mu_alpha, ux_seq, "+")) * plogis(as.vector(mu_anc_rho)) * pop15plus) /
  sum(plogis(as.vector(mu_anc_rho)) * pop15plus)

pred_ancartcov <- colSums(plogis(outer(mu_anc_alpha, ux_seq, "+")) *
                          plogis(as.vector(mu_anc_rho)) * anc_clients) /
  sum(plogis(as.vector(mu_anc_rho)) * anc_clients)

df_pred <- data.frame(pred_prev15to49, pred_ancprev, pred_artcov15plus, pred_ancartcov)


fig5a <- indicators %>%
  filter(
    indicator == "prevalence",
    sex == "both",
    age_group == "Y015_049",
    calendar_quarter == "CY2016Q1",
    area_level == 4
  ) %>%
  left_join(
    survey %>%
    filter(
      indicator == "prevalence",
      sex == "both",
      age_group == "Y015_049",
      survey_id %in% c("DEMO2015DHSA", "DEMO2016PHIA")
    )%>%
    select(area_id, survey_id, estimate),
    by = "area_id"
  ) %>%
  ggplot(aes(estimate, mean, ymin = lower, ymax = upper, shape = survey_id,
             name = area_name)) +
  geom_abline(linetype = "solid", color = "grey90") +
  geom_point(color = "grey25") +
  geom_linerange(color = "grey25") +
  scale_x_continuous(labels = label_percent(1), limits = c(0.02, 0.23)) +
  scale_y_continuous(labels = label_percent(1), limits = c(0.02, 0.23)) +
  scale_shape_manual(
    element_blank(),
    values = c("DEMO2015DHSA" = 16, "DEMO2016PHIA" = 15),
    labels = c("DEMO2015DHSA" = "MDHS 2015-16*", "DEMO2016PHIA" = "MPHIA 2015-16**")
  ) +
  labs(tag = "A",
       x = "Data: Survey HIV prevalence 15-49 y",
       y = "Model: Prevalence 15-49 y") +
  coord_fixed() +
  theme_bw(10) +
  theme(
    legend.position = c(1, 0),
    legend.just = c(1, 0),
    legend.background = element_rect(fill = NA),
    panel.grid = element_blank(),
    plot.tag = element_text(face = "bold")
  )

fig5b <- indicators %>%
  filter(
    indicator == "art_coverage",
    sex == "both",
    age_group == "Y015_999",
    calendar_quarter == "CY2016Q1",
    area_level == 4
  ) %>%
  left_join(
    survey %>%
    filter(
      indicator == "art_coverage",
      sex == "both",
      age_group == "Y015_064",
      survey_id %in% c("DEMO2015DHSA", "DEMO2016PHIA")
    )%>%
    select(area_id, survey_id, estimate),
    by = "area_id"
  ) %>%
  filter(!is.na(survey_id)) %>%
  ggplot(aes(estimate, mean, ymin = lower, ymax = upper, shape = survey_id,
             name = area_name)) +
  geom_abline(linetype = "solid", color = "grey90") +
  geom_point(color = "grey25") +
  geom_linerange(color = "grey25") +
  scale_x_continuous(labels = label_percent(1), limits = c(0.4, 0.88)) +
  scale_y_continuous(labels = label_percent(1), limits = c(0.4, 0.88)) +
  scale_shape_manual(
    element_blank(),
    values = c("DEMO2016PHIA" = 15),
    labels = c("DEMO2016PHIA" = "MPHIA 2015-16**")
  ) +
  labs(tag = "B",
       x = "Data: Survey ART coverage 15-64 y",
       y = "Model: ART coverage 15+ y") +
  coord_fixed() +
  theme_bw(10) +
  theme(
    legend.position = c(1, 0),
    legend.just = c(1, 0),
    legend.background = element_rect(fill = NA),
    panel.grid = element_blank(),
    plot.tag = element_text(face = "bold")
  )


fig5c <- indicators %>%
  filter(
    indicator == "prevalence",
    sex == "both",
    age_group == "Y015_049",
    calendar_quarter == "CY2016Q1",
    area_level == 4
  ) %>%
  left_join(
    naomi_data$anc_prev_t1_dat %>%
    mutate(anc_prev_obs = anc_prev_x / anc_prev_n) %>%
    select(area_id, anc_prev_obs),
    by = "area_id"
  ) %>%
  ggplot(aes(anc_prev_obs, mean, ymin = lower, ymax = upper,
             name = area_name,)) +
  geom_abline(linetype = "solid", color = "grey90") +
  geom_line(aes(pred_ancprev, pred_prev15to49), data = df_pred,
            color = "red4", linetype = "dashed", size = 0.8,
            inherit.aes = FALSE) +
  geom_point(color = "grey25") +
  geom_linerange(color = "grey25") +
  scale_x_continuous(labels = label_percent(1), limits = c(0.02, 0.23)) +
  scale_y_continuous(labels = label_percent(1), limits = c(0.02, 0.23)) +
  labs(tag = "C",
       x = "Data: ANC HIV prevalence",
       y = "Model: Prevalence 15-49 y") +
  coord_fixed() +
  theme_bw(10) +
  theme(
    panel.grid = element_blank(),
    plot.tag = element_text(face = "bold")
  )

fig5d <- indicators %>%
  filter(
    indicator == "art_coverage",
    sex == "both",
    age_group == "Y015_999",
    calendar_quarter == "CY2016Q1",
    area_level == 4
  ) %>%
  left_join(
    naomi_data$anc_artcov_t1_dat %>%
    mutate(anc_artcov_obs = anc_artcov_x / anc_artcov_n) %>%
    select(area_id, anc_artcov_obs),
    by = "area_id"
  ) %>%
  ggplot(aes(anc_artcov_obs, mean, ymin = lower, ymax = upper,
             name = area_name,)) +
  geom_abline(linetype = "solid", color = "grey90") +
  geom_line(aes(pred_ancartcov, pred_artcov15plus), data = df_pred,
            color = "red4", linetype = "dashed", size = 0.8,
            inherit.aes = FALSE) +
  geom_point(color = "grey25") +
  geom_linerange(color = "grey25") +
  scale_x_continuous(labels = label_percent(1), limits = c(0.4, 0.88)) +
  scale_y_continuous(labels = label_percent(1), limits = c(0.4, 0.88)) +
  labs(tag = "D",
       x = "Data: ANC ART coverage",
       y = "Model: ART coverage 15+ y") +
  coord_fixed() +
  theme_bw(10) +
  theme(
    panel.grid = element_blank(),
    plot.tag = element_text(face = "bold")
  )

fig5 <- grid.arrange(fig5a, fig5b, fig5c, fig5d)

ggsave("figure5.pdf", fig5, width = 180, height = 180, units = "mm")
ggsave("figure5.png", fig5, width = 180, height = 180, units = "mm")


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
  scale_color_brewer(palette = "Set1", guide = FALSE) +
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


quartz(w = 180 / 25.4, h = 180 / 25.4)

fig6 <- grid.arrange(fig6a, fig6b, fig6c)

ggsave("figure6.pdf", fig6, width = 180, height = 180, units = "mm")
ggsave("figure6.png", fig6, width = 180, height = 180, units = "mm")


#' ## Figure 7
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

fig7dat <- indicators %>%
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


fig7a <- fig7dat %>%
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

fig7b <- fig7dat %>%
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

fig7c <- fig7dat %>%
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



fig7d <- fig7dat %>%
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

fig7e <- fig7dat %>%
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

fig7f <- fig7dat %>%
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


fig7 <- grid.arrange(fig7a, fig7b, fig7c, fig7d, fig7e, fig7f,
                     nrow = 2)


ggsave("figure7.pdf", fig7, height = 130, width = 180, units = "mm")
ggsave("figure7.png", fig7, height = 130, width = 180, units = "mm")



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


fig7dat %>%
  select(indicator, area_id, area_name, age_group_label, version, mean, lower, upper) %>%
  arrange(indicator, area_name, version) %>%
  print(n = Inf)
