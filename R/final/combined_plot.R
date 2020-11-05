plot_df <- mutate(read_rds(here("R/final/01_noninformative-censoring-sim/noninf-cens-plot_df.RDS")),
                  Censoring = "non-Inf") %>%
    bind_rows(., mutate(read_rds("R/simple_sim/plotdf_inform_cens.RDS"),
                        Censoring = "Inf"))

plot_df <- plot_df %>% filter(estimand == "Relative Risk", Hazard != "Correct") %>%
    unite(Estimator, Hazard, col = "Estimator") %>%
    filter(Estimator %in% c("Kaplan-Meier_SuperLearner",
                            "G-Comp_Misspecified",
                            "TMLE_Misspecified",
                            "TMLE_SuperLearner")) %>%
    mutate(Estimator = case_when(Estimator == "Kaplan-Meier_SuperLearner" ~ "Kaplan-Meier",
                                 Estimator == "G-Comp_Misspecified" ~ "Misspecified G-Comp",
                                 Estimator == "TMLE_Misspecified" ~ "Misspecified TMLE",
                                 Estimator == "TMLE_SuperLearner" ~ "SuperLearner TMLE",
                                 T ~ Estimator))

inf_df <- plot_df %>% filter(Censoring == "Inf")
noninf_df <- plot_df %>% filter(Censoring == "non-Inf")

# plots -------------------------------------------------------------------


noninf_plot <- noninf_df %>%
    ggplot() + theme_minimal() +
    geom_boxplot(aes(y = estimate, x = Estimator), outlier.shape = NA) +
    geom_hline(aes(yintercept = truth),
               data = distinct(dplyr::select(noninf_df, estimand, truth)),
               colour = "blue", alpha = .7) +
    labs(title = "Estimator Performance with non-Informative Censoring",
         y = expr("Relative Risk\n")) +
    geom_label(aes(x = Estimator, y = estimate,
                   label = if_else(is.na(Cover), "NA",
                                   paste0(round(Cover, 2)*100, "%"))),
               data = summarise_all(group_by_if(noninf_df, ~!is.numeric(.)),
                                    ~quantile(., .95, na.rm = T)),
               size = 6, na.rm = T) +
    ylim(c(0.54, 0.83)) + xlab(expr("\nEstimators")) +
    theme(title = element_text(size = 20, face = "bold"),
          axis.title = element_text(face = "bold", size = 16),
          axis.text = element_text(size = 16))
ggsave(filename = "noninf_plot.png", path = "R/final/",
       device = "png", width = 11, height = 6, units = "in", noninf_plot)

inf_plot <- inf_df %>%
    ggplot() + theme_minimal() +
    geom_boxplot(aes(y = estimate, x = Estimator), outlier.shape = NA) +
    geom_hline(aes(yintercept = truth),
               data = distinct(dplyr::select(inf_df, estimand, truth)),
               colour = "blue", alpha = .7) +
    labs(title = "Estimator Performance with Informative Censoring",
         y = expr("Relative Risk\n")) +
    geom_label(aes(x = Estimator, y = estimate,
                   label = if_else(is.na(Cover), "NA",
                                   paste0(round(Cover, 2)*100, "%"))),
               data = summarise_all(group_by_if(inf_df, ~!is.numeric(.)),
                                    ~quantile(., .95, na.rm = T)),
               size = 6, na.rm = T) +
    ylim(c(0.8, 1.29)) + xlab(expr("\nEstimators")) +
    theme(title = element_text(size = 20, face = "bold"),
          axis.title = element_text(face = "bold", size = 16),
          axis.text = element_text(size = 16))
ggsave(filename = "inf_plot.png", path = "R/final/",
       device = "png", width = 11, height = 6, units = "in", inf_plot)


combined_plot <- gridExtra::arrangeGrob(noninf_plot, inf_plot)
ggsave(filename = "combined_plot.png", path = "R/final/",
       device = "png", width = 6, height = 9, units = "in", combined_plot)
