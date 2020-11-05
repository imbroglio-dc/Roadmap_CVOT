
# setup & data ------------------------------------------------------------


library(tidyverse); library(mice); library(survminer); library(survival)
library(SuperLearner); library(survtmle); library(MOSS); library(here)
source(file = here("R/functions/my_sl_functions.R"))
source(file = here("R/functions/my_MOSS_hazard_methods.R"))
source(file = here("R/functions/my_survtmle_functions.R"))

set.seed(0)

covar_names <- head(colnames(readxl::read_excel("data/test_leader.xlsx")), -3)
covar_names[c(1, 5)] <- c("ISREGION", "EGFREPB")

obs <- haven::read_sas("data/ADaM/adtte.sas7bdat") %>%
    filter(PARAMCD == "MACEEVTM") %>%
    left_join(., haven::read_sas("data/ADaM/adsl.sas7bdat")) %>%
    select(AVAL, CNSR, all_of(covar_names)) %>%
    mutate_if(is_character, as_factor) %>%
    mutate_if(~length(levels(.)) == 2, ~as.logical(as.numeric(.)-1)) %>%
    mutate(ARM = as.numeric(ARM), CNSR = 1-CNSR, AVAL = ceiling(AVAL)) %>%
    rename(TIME = AVAL, EVENT = CNSR) %>%
    select(TIME, EVENT, ARM, everything()) %>%
    mutate(BMI_imp = is.na(BMIBL),  # indicators for missing values
           AGE_imp = is.na(AGE))

obs <- complete(mice(obs, m = 10, maxit = 10))

# km ----------------------------------------------------------------------

results <- list(estimates = NULL)
targets <- (1:4)*12

km_fit <- survival::survfit(survival::Surv(obs$TIME,
                                           obs$EVENT,
                                           type = "right") ~ ARM,
                            type = "kaplan-meier", data = obs)
km_est <- summary(km_fit, times = targets)

results$estimates <- tibble(t = round(km_est$time, 1)) %>%
    mutate(A = rep(c(0, 1), each = length(t)/2),
           s = km_est$surv, # s0, s1
           se = km_est$std.err) %>%
    pivot_wider(names_from = "A", values_from = c('s', 'se'), names_sep = "") %>%
    cbind(Estimator = "Kaplan-Meier", .)


# hazard estim ------------------------------------------------------------


sl_lib_g <- c("SL.glm")
sl_lib_censor <- c("SL.mean", "sl_glm", "SL.glm.interaction", "SL.glmnet",
                   "SL.ranger", "sl_xgboost", "sl_bayesglm")
sl_lib_failure <- c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.glmnet",
                    "SL.ranger", "sl_xgboost", "SL.bayesglm")

# create a baseline covariates data.frame
W <- dplyr::select(obs, -c(TIME, EVENT, ARM))

# coarsen the time scale of the obs dataset
# timescale <- 4 # time coarsening factor
# obs <- obs %>%
#     mutate(time = ceiling(time/timescale))

# get SuperLearner estimate for failure hazard and cumulative censoring probability.
# See novo.nordisk:::init_sl_fit for other returned values, but the failure event
# hazards and censoring probabilities are needed for the surv_tmle function)
set.seed(256)



sl_fit <- my_init_sl_fit(
    T_tilde = obs$TIME,
    Delta = as.numeric(obs$EVENT),
    A = as.numeric(obs$ARM),
    W = W,
    t_max = max(obs$TIME),
    sl_failure = sl_lib_failure,
    sl_censoring = sl_lib_censor,
    sl_treatment = sl_lib_g,
    cv.Control = list(V = 10))
saveRDS(sl_fit, here("R/final/03_LEADER/sl_fit.RDS"), compress = T)

haz_sl <- list(sl_fit$density_failure_1$clone(),
               sl_fit$density_failure_0$clone())
haz_sl[[1]]$haz2surv()
haz_sl[[2]]$haz2surv()
names(haz_sl) <- c("A = 1", "A = 0")

results$estimates <- data.frame(s1 = colMeans(haz_sl[[1]]$survival),
                                s0 = colMeans(haz_sl[[2]]$survival)) %>%
    mutate(t = (1:length(s1))) %>% filter(t %in% targets) %>%
    cbind(Estimator = "G-Comp: SuperLearner", .) %>%
    bind_rows(results$estimates, .)


# tmle --------------------------------------------------------------------

SL_ftime <- sl_fit$models$Y
sl_G_dC <- sl_fit$G_dC
glm_trt <- paste0(colnames(W), collapse = " + ")

# TMLE

tmle_sl <- surv_tmle(ftime = obs$TIME,
                     ftype = obs$EVENT,
                     targets = targets,
                     trt = obs$ARM,
                     t0 = targets, adjustVars = W,
                     SL.ftime = SL_ftime, SL.ctime = sl_G_dC,
                     glm.trt = glm_trt,
                     returnIC = T, returnModels = T,
                     ftypeOfInterest = 1, trtOfInterest = c(1, 0),
                     maxIter = 25, method = "hazard")

tmle_sl_out <- suppressWarnings(
    t(1 - tmle_sl$est) %>% unname() %>% as_tibble() %>%
        cbind(Estimator = "TMLE: Hazard - SL", t = targets, .) %>%
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
                 var(colSums(c(1, -1) * ic[c(i, i+(nrow(ic)/2)), ]))
             }),
             RR_se = sapply(1:(nrow(ic)/2), function(i) {
                 s0 = tmle_sl_out[i, "s0"] ; s1 = tmle_sl_out[i, "s1"]
                 var(colSums(c((1-s1)/(1-s0)^2, -1/(1-s0)) *
                                 ic[c(i, i+(nrow(ic)/2)), ]))
             }),
             SR_se = sapply(1:(nrow(ic)/2), function(i) {
                 s0 = tmle_sl_out[i, "s0"] ; s1 = tmle_sl_out[i, "s1"]
                 var(colSums(c(-(s1)/(s0)^2, 1/s0) * ic[c(i, i+(nrow(ic)/2)), ]))
             })) %>% apply(., 2, function(j) sqrt(j / ncol(ic))) %>%
    as_tibble() %>% bind_cols(t = targets, .) %>%
    pivot_longer(-`t`, names_to = c("Estimand", "tmp"), names_sep = "_",
                 values_to = "se") %>% dplyr::select(-tmp) %>%
    cbind(Estimator = "TMLE", .)

# plot results ------------------------------------------------------------

result_tbl <- results$estimates %>%
    mutate(RD = (1-s1) - (1-s0),
           RR = (1-s1) / (1-s0),
           SR = s1 / s0,
           s0_se = se0, s1_se = se1,
           RR_se = sqrt(se1^2 / (1-s0)^2 + se0^2 * ((1-s1) / (1-s0)^2)^2),
           RD_se = sqrt(se1^2 + se0^2),
           SR_se = sqrt(se1^2 / s0^2 + se0^2 * s1^2 / s0^4),
           Estimator = as.character(Estimator),
           Estimator = case_when(Estimator == "G-Comp: SuperLearner" ~ "G-Comp",
                                 Estimator == "TMLE: Hazard - SL" ~ "TMLE",
                                 T ~ Estimator)) %>% select(-c(se0, se1)) %>%
    pivot_longer(cols = c(`s0`, `s1`, `RD`, `RR`, `SR`), names_to = "Estimand",
                 values_to = "Estimate") %>%
    pivot_longer(cols = contains("se"), names_to = c("se_est", "tmp"),
                 names_sep = "_", values_to = "se") %>% dplyr::select(-`tmp`) %>%
    filter(se_est == Estimand) %>% dplyr::select(-se_est) %>%
    full_join(., se, by = c("Estimator", "t", "Estimand")) %>%
    mutate(se = case_when(is.na(se.y) ~ se.x,
                          T ~ se.y)) %>%
    dplyr::select(-c(se.x, se.y)) %>%
    pivot_wider(names_from = Estimator, names_sep = "_",
                values_from = c(se, Estimate)) %>%
    mutate(Eff = `se_Kaplan-Meier`^2 / `se_TMLE`^2) %>%
    pivot_longer(-c(`t`, Estimand, Eff), names_sep = "_",
                 names_to = c(".value", "Estimator")) %>%
    mutate(Eff = case_when(Estimator == "TMLE" ~ Eff,
                           Estimator == "G-Comp" ~ NaN,
                           T ~ 1)) %>%
    dplyr::select(Estimator, t, Estimand, Estimate, se, Eff) %>%
    mutate(t = as.character(`t`))
saveRDS(object = result_tbl, here("R/final/03_LEADER/LEADER-estimates.RDS"))


result_plot <- result_tbl %>% filter(Estimand  == "RR") %>%
    ggplot(aes(x=`t`, y = Estimate, colour = Estimator)) +
    geom_errorbar(aes(ymin = Estimate - 1.96*se,
                      ymax = Estimate + 1.96*se), width = .5,
                  position = position_dodge(.5)) +
    geom_point(position = position_dodge(.5)) + theme_minimal() +
    labs(x = "Months", y = "Relative Risk", title = "LEADER re-Analysis") +
    geom_label(aes(x =`t`, y = Estimate - 1.96*`se` - .02,
                   label = paste0("Eff = ", round(Eff, 3)*100, "%")),
               data = filter(result_tbl, Estimand  == "RR", Estimator == "TMLE"),
               colour = 'black')

ggsave(filename = "leader.png", path = "R/final/03_LEADER/",
       device = "png", width = 9, height = 6, units = "in", result_plot)

results$estimates %>%
    mutate(RR = (1-s1) / (1-s0),
           RR_se = sqrt(se1^2 / (1-s0)^2 + se0^2 * ((1-s1) / (1-s0)^2)^2),
           RR_se_ic = c(head(RR_se, -4), se$RR_se),
           Estimator = as.character(Estimator),
           Estimator = case_when(Estimator == "G-Comp: SuperLearner" ~ "G-Comp",
                                 Estimator == "TMLE: Hazard - SL" ~ "TMLE",
                                 T ~ Estimator),
           `t` = as.character(`t`)) %>%
    pivot_longer(cols = c(`s0`, `s1`, `RR`), names_to = "Estimand",
                 values_to = "Estimate") %>%
    pivot_longer(cols = c(`se0`, `se1`, `RR_se_ic`), names_to = "se_est",
                 values_to = "se") %>%
    mutate(se_est = case_when(se_est == "se0" ~ "s0",
                              se_est == "se1" ~ "s1",
                              T ~ "RR")) %>%
    filter(Estimand == se_est) %>% select(-se_est) -> tmp

tmp %>%
    ggplot(aes(x=`t`, y = Estimate, colour = Estimator)) +
    facet_wrap(~Estimand, scales = "free", ncol = 1) +
    geom_errorbar(aes(ymin = Estimate - 1.96*se,
                      ymax = Estimate + 1.96*se), width = .5,
                  position = position_dodge(.5)) +
    geom_point(position = position_dodge(.5)) + theme_minimal() +
    labs(x = "Months", y = "Estimates", title = "LEADER re-Analysis")

ggsave(filename = "leader+.png", path = "R/final/03_LEADER/",
       device = "png", width = 6, height = 9, units = "in")




# # coxph  ------------------------------------------------------------------
#
# library(survival)
# cox_fit <- coxph(Surv(time = TIME, event = EVENT, type = "right") ~ ., data = obs)
#
# surv1 <- t(survival:::survfit.coxph(cox_fit, newdata = dplyr::mutate(ARM = 1, obs),
#                                     stype = 1)$surv)
# surv1 <- survival_curve$new(t = 1:ncol(surv1), survival = surv1)
# surv1$surv2haz()
#
# surv0 <- t(survival:::survfit.coxph(cox_fit, newdata =
#                                       dplyr::mutate(ARM = 0, obs), stype = 1)$surv)
# surv0 <- survival_curve$new(t = 1:ncol(surv0), survival = surv0)
# surv0$surv2haz()
#
# tmp <- glmnet:::glmnet(x = model.matrix(~ -1 + ., obs[, 3:16]),
#                        y = Surv(obs$TIME, obs$EVENT, type = "right"),
#                        family = "cox")
#
