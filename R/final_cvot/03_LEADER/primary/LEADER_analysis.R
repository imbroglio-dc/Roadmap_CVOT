library(tidyverse); library(here); library(survival); library(mice);
library(SuperLearner); library(survtmle); library(MOSS)
source(file = here("R/functions/my_sl_functions.R"))
source(file = here("R/functions/my_MOSS_hazard_methods.R"))
source(file = here("R/functions/my_survtmle_functions.R"))
i_am("R/final_cvot/03_LEADER/primary/LEADER_analysis.R")

set.seed(0)

## time coarsening ----
timescale <- 4 # time coarsened to 4 month intervals

## define target times ----
targets <- 1:4 * 12 / timescale # target end of years 1-4

## data cleaning ----
### cleaned covariates ----
source(file = here("R/functions/LEADER_W_clean_new.R"))

obs <- dplyr::select(outcomes, USUBJID, AVAL_MACE, EVENT_MACE) %>%
    left_join(., W_imputed) %>%
    rename(TIME = AVAL_MACE, EVENT = EVENT_MACE) %>%
    dplyr::select(-USUBJID) %>%
    mutate(TIME = ceiling(TIME / timescale))

# estimation ----
results <- list(estimates = NULL)
adjust_vars <- model.matrix(~. + -1, dplyr::select(obs, -c(TIME, EVENT, ARM))) %>%
    as_tibble()
adjust_vars_names <- colnames(adjust_vars)
adjust_vars <- adjust_vars %>% rename_all(~paste0("L", 1:ncol(adjust_vars)))
names(adjust_vars_names) <- colnames(adjust_vars)

## km ----
km_fit <- survival::survfit(survival::Surv(obs$TIME,
                                           obs$EVENT,
                                           type = "right") ~ ARM,
                            type = "kaplan-meier", data = obs)
km_est <- summary(km_fit, times = targets)

results$estimates <- tibble(t = km_est$time) %>%
    mutate(A = rep(c(0, 1), each = length(t)/2),
           s = km_est$surv, # s0, s1
           se = km_est$std.err) %>%
    pivot_wider(names_from = "A", values_from = c('s', 'se'), names_sep = "") %>%
    cbind(Estimator = "Kaplan-Meier", .)


# cox RR estimate -----------------------------------------------------------------------------

cox_fit <- coxph(Surv(obs$TIME, obs$EVENT) ~ ARM, data = mutate(obs, ARM = factor(ARM)))
summary(cox_fit)


# ldr_cox <- survfit(cox_fit,
#                    newdata = data.frame("ARM" = unique(tmpobs$ARM)))
#
# ldr_cox_rr <- summary(ldr_cox, times = 1:4 * 12)$surv %>% as_tibble() %>%
#     mutate("t" = as.character(1:4 * 12),
#            Estimator = "G-Comp",
#            RR = (1 - `2`) / (1 - `1`))
#
# res_plot <- result_tbl %>% filter(Estimand  == "RR")
# res_plot[res_plot$Estimator == "G-Comp", "Estimate"] <- ldr_cox_rr$RR
# res_plot <- res_plot %>%
#     ggplot(aes(x=`t`, y = Estimate, colour = Estimator)) +
#     geom_errorbar(aes(ymin = Estimate - 1.96*se,
#                       ymax = Estimate + 1.96*se), width = .5,
#                   position = position_dodge(.5)) +
#     geom_point(position = position_dodge(.5)) + theme_minimal() +
#     labs(x = "Months", y = "Relative Risk", title = "LEADER re-Analysis")
#
# ggsave(filename = "leader.png", path = "R/final_cvot/03_LEADER/",
#        device = "png", width = 9, height = 6, units = "in", result_plot)

## TMLE ----
### Define superlearner learners and libraries ----
glmnets <- create.Learner(base_learner = "SL.glmnet",
                          tune = list("alpha" = c(0, 0.5, 1)),
                          detailed_names = T)

leader_cens_glm <- function(Y, X, newX, family, obsWeights, model = TRUE, ...)
{   # indicators for == 14, 15, 16 are for timescale = 3, i.e. 3month time intervals
    # indicators for == 11, 12 are for timescale = 4, i.e. 4 month time intervals
    if (is.matrix(X)) {
        X = as.data.frame(X)
    }
    fit.glm <- glm(Y ~ I(t==11) + I(t==12) + ., data = X,
                   family = family, weights = obsWeights, model = model)
    if (is.matrix(newX)) {
        newX = as.data.frame(newX)
    }
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred = pred, fit = fit)
    return(out)
}
environment(leader_cens_glm) <- asNamespace('SuperLearner')

leader_cens_bayesglm <- function(Y, X, newX, family, obsWeights, ...)
{   # indicators for == 14, 15, 16 are for timescale = 3, i.e. 3 month time intervals
    # indicators for == 11, 12 are for timescale = 4, i.e. 4 month time intervals
    .SL.require("arm")
    if (is.matrix(X)) {
        X = as.data.frame(X)
    }
    if (is.matrix(newX)) {
        newX = as.data.frame(newX)
    }
    fit.glm <- arm::bayesglm(Y ~ I(t==11) + I(t==12) + .,
                             data = X, family = family, weights = obsWeights)
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("SL.bayesglm")
    return(out)
}
environment(leader_cens_bayesglm) <- asNamespace("SuperLearner")


sl_lib_g <- c("SL.glm", glmnets$names, "SL.bayesglm")
sl_lib_censor <- c("SL.glm", glmnets$names, "SL.bayesglm",
                   "leader_cens_glm", "leader_cens_bayesglm")
sl_lib_failure <- c("SL.glm", glmnets$names, # "SL.ranger",
                    "SL.xgboost", "SL.bayesglm")


### hazard estimation ----
# get SuperLearner estimate for failure hazard and cumulative censoring probability.
# See novo.nordisk:::init_sl_fit for other returned values, but the failure event
# hazards and censoring probabilities are needed for the surv_tmle function)
set.seed(256)
if (file.exists(here("R/03_LEADER/primary/sl_fit-fullcov.RDS"))) {
    sl_fit <- readRDS(here("R/03_LEADER/primary/sl_fit-fullcov.RDS"))
} else {
    if (file.exists(here("R/final_cvot/03_LEADER/primary/output/sl_fit_primary.RDS"))) {
        sl_fit <- readRDS(here("R/final_cvot/03_LEADER/primary/output/sl_fit_primary.RDS"))
    } else {
        sl_fit <- my_init_sl_fit(
            T_tilde = obs$TIME,
            Delta = as.numeric(obs$EVENT),
            A = as.numeric(obs$ARM),
            W = adjust_vars,
            t_max = max(targets),
            sl_failure = sl_lib_failure,
            sl_censoring = sl_lib_censor,
            sl_treatment = "SL.glm",
            cv.Control = list(V = 10))
        saveRDS(sl_fit, here("R/final_cvot/03_LEADER/primary/output/sl_fit_primary.RDS"), compress = F)
    }

    haz_sl <- list(sl_fit$density_failure_1$clone(),
                   sl_fit$density_failure_0$clone())
    haz_sl[[1]]$haz2surv()
    haz_sl[[2]]$haz2surv()
    names(haz_sl) <- c("A = 1", "A = 0")

    saveRDS(sl_fit, here("R/03_LEADER/primary/sl_fit-fullcov.RDS"), compress = F)
}


    # tmle --------------------------------------------------------------------

    SL_ftime <- sl_fit$models$Y
    sl_G_dC <- sl_fit$G_dC
    glm_trt <- paste0(colnames(adjust_vars), collapse = " + ")
    rm(sl_fit)

    # TMLE

    tmle_sl <- surv_tmle(ftime = obs$TIME,
                         ftype = obs$EVENT,
                         targets = targets,
                         trt = obs$ARM,
                         t0 = max(targets), adjustVars = adjust_vars,
                         SL.ftime = SL_ftime, SL.ctime = sl_G_dC,
                         SL.trt = sl_lib_g,
                         returnIC = T, returnModels = T,
                         ftypeOfInterest = 1, trtOfInterest = c(1, 0),
                         maxIter = 25, method = "hazard")

    tmle_sl_out <- suppressWarnings(
        t(1 - tmle_sl$est) %>% unname() %>% as_tibble() %>%
            cbind(Estimator = "TMLE", t = targets, .) %>%
            rename("s0" = V1, "s1" = V2))
    results$estimates <- suppressWarnings(
        matrix(sqrt(diag(tmle_sl$var)),
               nrow = length(targets), byrow = F,
               dimnames = list(NULL, c("se0", "se1"))) %>%
            as.data.frame() %>% cbind(tmle_sl_out, .)) %>%
        bind_rows(results$estimates, .)

    # IC-based variance estimation --------------------------------------------

    ic <- t(tmle_sl$ic)


    se <- tibble(s0_se = apply(ic[1:(nrow(ic)/2), ], 1, var),
                 s1_se = apply(ic[(nrow(ic)/2 + 1):nrow(ic), ], 1, var),
                 RD_se = sapply(1:(nrow(ic)/2), function(i) {
                     var(colSums(c(1, -1) * ic[c(i, i + (nrow(ic)/2)), ]))
                 }),
                 RR_se = sapply(1:(nrow(ic)/2), function(i) {
                     s0 = tmle_sl_out[i, "s0"] ; s1 = tmle_sl_out[i, "s1"]
                     var(colSums(c((1 - s1)/(1 - s0)^2, -1/(1 - s0)) *
                                     ic[c(i, i + (nrow(ic)/2)), ]))
                 }),
                 SR_se = sapply(1:(nrow(ic)/2), function(i) {
                     s0 = tmle_sl_out[i, "s0"] ; s1 = tmle_sl_out[i, "s1"]
                     var(colSums(c(-(s1)/(s0)^2, 1/s0) * ic[c(i, i + (nrow(ic)/2)), ]))
                 })) %>% apply(., 2, function(j) sqrt(j / ncol(ic))) %>%
        as_tibble() %>% bind_cols(t = targets, .) %>%
        pivot_longer(-`t`, names_to = c("Estimand", "tmp"), names_sep = "_",
                     values_to = "se") %>% dplyr::select(-tmp) %>%
        cbind(Estimator = "TMLE", .)
    rm(tmle_sl)


se <- tibble(s0_se = apply(ic[1:(nrow(ic)/2), ], 1, var),
             s1_se = apply(ic[(nrow(ic)/2 + 1):nrow(ic), ], 1, var),
             RD_se = sapply(1:(nrow(ic)/2), function(i) {
                 var(colSums(c(1, -1) * ic[c(i, i + (nrow(ic)/2)), ]))
             }),
             RR_se = sapply(1:(nrow(ic)/2), function(i) {
                 s0 = tmle_sl_out[i, "s0"] ; s1 = tmle_sl_out[i, "s1"]
                 var(colSums(c((1 - s1)/(1 - s0)^2, -1/(1 - s0)) *
                                 ic[c(i, i + (nrow(ic)/2)), ]))
             }),
             SR_se = sapply(1:(nrow(ic)/2), function(i) {
                 s0 = tmle_sl_out[i, "s0"] ; s1 = tmle_sl_out[i, "s1"]
                 var(colSums(c(-(s1)/(s0)^2, 1/s0) * ic[c(i, i + (nrow(ic)/2)), ]))
             })) %>% apply(., 2, function(j) sqrt(j / ncol(ic))) %>%
    as_tibble() %>% bind_cols(t = targets, .) %>%
    pivot_longer(-`t`, names_to = c("Estimand", "tmp"), names_sep = "_",
                 values_to = "se") %>% dplyr::select(-tmp) %>%
    cbind(Estimator = "TMLE", .)
rm(tmle_sl)

# plot results ------------------------------------------------------------

if (!file.exists(here("R/03_LEADER/primary/LEADER-estimates-tbl-fullcov.RDS"))) {
    result_tbl <- results$estimates %>%
        mutate(RD = (1 - s1) - (1 - s0),
               RR = (1 - s1) / (1 - s0),
               SR = s1 / s0,
               s0_se = se0, s1_se = se1,
               RR_se = sqrt(se1^2 / (1 - s0)^2 + se0^2 * ((1 - s1) / (1 - s0)^2)^2),
               RD_se = sqrt(se1^2 + se0^2),
               SR_se = sqrt(se1^2 / s0^2 + se0^2 * s1^2 / s0^4),
               Estimator = as.character(Estimator)) %>% select(-c(se0, se1)) %>%
        pivot_longer(cols = c(`s0`, `s1`, `RD`, `RR`, `SR`), names_to = "Estimand",
                     values_to = "Estimate") %>%
        pivot_longer(cols = contains("se"), names_to = c("se_est", "tmp"),
                     names_sep = "_", values_to = "se") %>% dplyr::select(-`tmp`) %>%
        filter(se_est == Estimand) %>% dplyr::select(-se_est) %>%
        full_join(., se, by = c("Estimator", "t", "Estimand")) %>%
        mutate(se = case_when(is.na(se.y) ~ se.x,
                              T ~ se.y)) %>%
        dplyr::select(-c(se.x, se.y)) %>% distinct() %>%
        pivot_wider(names_from = Estimator, names_sep = "_",
                    values_from = c(se, Estimate)) %>%
        mutate(Eff = `se_Kaplan-Meier`^2 / `se_TMLE`^2) %>%
        pivot_longer(-c(`t`, Estimand, Eff), names_sep = "_",
                     names_to = c(".value", "Estimator")) %>%
        mutate(Eff = case_when(Estimator == "TMLE" ~ Eff,
                               Estimator == "G-Comp" ~ NaN,
                               T ~ 1),
               t = t * timescale) %>%
        dplyr::select(Estimator, t, Estimand, Estimate, se, Eff)
    saveRDS(object = result_tbl, here("R/03_LEADER/primary/LEADER-estimates-tbl-fullcov.RDS"))
} else {
    result_tbl <- readRDS(here("R/03_LEADER/primary/LEADER-estimates-tbl-fullcov.RDS"))
}


result_plot <- result_tbl %>% filter(Estimand == "RR", Estimator != "G-Comp") %>%
    ggplot(aes(x = as.character(t), y = Estimate, colour = Estimator)) +
    geom_errorbar(aes(ymin = Estimate - 1.96*se,
                      ymax = Estimate + 1.96*se), width = .5,
                  position = position_dodge(.5)) +
    geom_point(position = position_dodge(.5)) + theme_minimal() +
    labs(x = "Months", y = "Relative Risk", title = "LEADER MACE Analysis") +
    geom_label(aes(y = Estimate - 1.96 * `se`,
                   label = paste0("Eff = ", round(Eff, 2)*100, "%")), nudge_y = -.025,
               data = filter(result_tbl, Estimand  == "RR", Estimator == "TMLE"),
               colour = 'black')
ggsave(filename = "leader.png", path = "R/03_LEADER/primary/",
       plot = result_plot, device = "png", width = 9, height = 6, units = "in")


# result_plot <- result_tbl %>% filter(Estimand == "RD", Estimator != "G-Comp") %>%
#     ggplot(aes(x = as.character(`t`), y = Estimate, colour = Estimator)) +
#     geom_errorbar(aes(ymin = Estimate - 1.96*se,
#                       ymax = Estimate + 1.96*se), width = .5,
#                   position = position_dodge(.5)) +
#     geom_point(position = position_dodge(.5)) + theme_minimal() +
#     labs(x = "\nTime (Months)", y = "Risk Difference\n",
#          title = "LEADER Analysis: Composite Event (MACE)") +
#     theme(axis.text=element_text(size=12),
#           axis.title=element_text(size=16)) +
# geom_label(aes(x = as.character(`t`), y = Estimate - 1.96 * `se` - 0.001,
#                label = paste0("Eff = ", round(Eff, 3)*100, "%")),
#            data = filter(result_tbl, Estimand  == "RD", Estimator == "TMLE"),
#            colour = 'black')
# result_plot
#
# ggsave(filename = "leader_MACE.png", path = "R/competing_risks/LEADER/nfMI/",
#        plot = result_plot, device = "png", width = 9, height = 6, units = "in")


tmp <- result_tbl %>% filter(Estimand == "RR", Estimator != "G-Comp", `t` == 48) %>%
    mutate("ci" = 1.96*2*se) %>%
    rbind(., data.frame("Estimator" = "Cox HR", 't' = NA, 'Estimand' = "HR",
                        "Estimate" = exp(cox_fit$coefficients[1]), "se" = exp(sqrt(cox_fit$var[1,1])),
                        "Eff" = 1, "ci" = exp(cox_fit$coefficients[1] + 1.96 * sqrt(cox_fit$var[1,1])) -
                            exp(cox_fit$coefficients[1] - 1.96 * sqrt(cox_fit$var[1,1])))) %>%
    mutate(ci = ci/tail(ci, 1)) %>%
    mutate(p.val = pnorm(Estimate, 1, se)*2, upper = round(Estimate + 1.96*se, 2), lower = round(Estimate - 1.96*se, 2))

result_tbl %>% filter(Estimand == "RR", Estimator != "G-Comp") %>% group_by(`t`) %>%
    mutate(relKM_ci = 1.96*se*2,
           relKM_ci = relKM_ci / head(relKM_ci, 1),
           rel_se = sqrt(Eff))



combined_result_plot <- result_tbl %>%
    filter(Estimand == "RR", Estimator != "G-Comp", `t` == 48) %>%
    mutate(l = Estimate - 1.96*se, u = Estimate + 1.96*se) %>%
    dplyr::select(Estimator, Estimand, Estimate, l, u) %>%
    rbind(., data.frame("Estimator" = "Cox", 'Estimand' = "HR",
                        "Estimate" = exp(cox_fit$coefficients[1]),
                        "l" = exp(cox_fit$coefficients[1] - 1.96 * sqrt(cox_fit$var[1,1])),
                        "u" = exp(cox_fit$coefficients[1] + 1.96 * sqrt(cox_fit$var[1,1]))
    )) %>% mutate(relCI = (u - l) / tail(u - l, 1)) %>%
    ggplot(aes(x = Estimator, y = Estimate, colour = Estimator)) +
    geom_errorbar(aes(ymin = l, ymax = u), width = 0.8) +
    geom_point(size = 2) + theme_minimal() +
    geom_hline(aes(yintercept = 1), alpha = 0.3, colour = "black") +
    geom_label(aes(y = l, label = paste0("Rel CI = ", round(relCI, 2)*100, "%")),
               nudge_y = -.03, colour = 'black') +
    scale_y_continuous(name = "Hazard Ratio\n",
                       sec.axis = sec_axis(~., name="Year 4 Relative Risk\n")) +
    labs(x = "\nEstimator",  title = "LEADER MACE Analysis") +
    # theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
    theme(legend.position="none")
combined_result_plot
ggsave(filename = "leader-primary+cox.png", path = "R/final_cvot/03_LEADER/primary/output/",
       plot = combined_result_plot, device = "png", width = 6, height = 4, units = "in")



res_plot <- result_tbl %>% filter(Estimand  == "RR")
res_plot[res_plot$Estimator == "G-Comp", "Estimate"] <- ldr_cox_rr$RR
res_plot <- res_plot %>% mutate(Estimator = case_when(Estimator == "G-Comp" ~ "G-Comp: Cox",
                                                      T ~ Estimator)) %>%
    ggplot(aes(x = as.character(`t`), y = Estimate, colour = Estimator)) +
    geom_errorbar(aes(ymin = Estimate - 1.96*se,
                      ymax = Estimate + 1.96*se), width = .5,
                  position = position_dodge(.5)) +
    geom_point(position = position_dodge(.5)) + theme_minimal() +
    labs(x = "Months", y = "Relative Risk") #, title = "LEADER re-Analysis"

# ggsave(filename = "leader.png", path = "R/03_LEADER/",
#        device = "png", width = 9, height = 6, units = "in", result_plot)



# relative CIs --------------------------------------------------------------------------------

rr_tbl <- result_tbl %>% filter(Estimand == "RR", Estimator != "G-Comp") %>%
    dplyr::select(-Estimand) %>%
    mutate(ci = 2 * 1.96 * se) %>%
    rbind(data.frame("Estimator" = "Cox", "t" = NA,
                     "Estimate" = exp(ldr_cox$coef),
                     "se" = NA, "Eff" = NA,
                     "ci" = diff(exp(ldr_cox$coef +
                                         sqrt(ldr_cox$var[1, 1]) * 1.96 * c(-1, 1))))) %>%
    mutate(relCI = ci / tail(ci, 1))

