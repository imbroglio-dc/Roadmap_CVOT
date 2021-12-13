
# setup & data ----
# devtools::install_github("wilsoncai1992/MOSS")

## load libraries & helper functions ----
library(mice); library(survival); library(SuperLearner); library(survtmle); library(MOSS)
library(tidyverse); library(here)
source(file = here("R/functions/my_sl_functions.R"))
source(file = here("R/functions/my_MOSS_hazard_methods.R"))
source(file = here("R/functions/my_survtmle_functions.R"))

set.seed(0)

## clean and combine LEADER data ----
source(file = here("R/functions/LEADER_W_clean.R"))

## time coarsening ----
timescale <- 4 # time coarsened to 4 month intervals

## define target times ----
targets <- 1:4 * 12 / timescale # target end of years 1-4

obs <- haven::read_sas("data/ADaM/adtte.sas7bdat") %>%
    filter(PARAMCD == "MACEEVTM") %>%
    dplyr::select(AVAL, CNSR, USUBJID) %>% left_join(., W) %>%
    dplyr::select(AVAL, CNSR, ARM, everything(), -USUBJID) %>%
    mutate_if(is_character, as_factor) %>%
    mutate_if(~length(levels(.)) == 2, ~as.logical(as.numeric(.) - 1)) %>%
    mutate(ARM = as.numeric(ARM), CNSR = 1 - CNSR, AVAL = ceiling(AVAL/timescale)) %>%
    rename(TIME = AVAL, EVENT = CNSR) %>%
    dplyr::select(-c(NYHACLAS, RETINSEV)) # removing NYHACLAS and RETINSEV for missingness

# imputation
obs <- dplyr::select(obs, -c(TIME, EVENT, ARM, INSNVFL, CHFFL)) %>%
    mice::mice(m = 1, maxit = 10) %>% mice::complete(10L) %>%
    cbind(., dplyr::select(obs, c(TIME, EVENT, ARM, INSNVFL, CHFFL))) %>%
    dplyr::select(TIME, EVENT, ARM, everything())

obs <- W %>% select_if(anyNA) %>% rename_all(~paste0(., ".imp")) %>%
    mutate_all(is.na) %>% cbind(obs, .)


# estimation ----
results <- list(estimates = NULL)
adjust_vars <- model.matrix(~-1 + ., dplyr::select(obs, -c(TIME, EVENT, ARM))) %>%
    t() %>% as.data.frame() %>% distinct() %>% t() %>% as_tibble()
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


sl_lib_g <- c("SL.glm")
sl_lib_censor <- c("SL.glm", glmnets$names, "SL.bayesglm",
                   "leader_cens_glm", "leader_cens_bayesglm")
sl_lib_failure <- c("SL.glm", glmnets$names, # "SL.ranger", "SL.xgboost",
                    "SL.bayesglm")


### hazard estimation ----
# get SuperLearner estimate for failure hazard and cumulative censoring probability.
# See novo.nordisk:::init_sl_fit for other returned values, but the failure event
# hazards and censoring probabilities are needed for the surv_tmle function)
set.seed(256)
if (file.exists(here("R/final_cvot/03_LEADER/sl_fit.RDS"))) {
    sl_fit <- readRDS(here("R/final_cvot/03_LEADER/sl_fit.RDS"))
} else {
    sl_fit <- my_init_sl_fit(
        T_tilde = obs$TIME,
        Delta = as.numeric(obs$EVENT),
        A = as.numeric(obs$ARM),
        W = adjust_vars,
        t_max = max(targets),
        sl_failure = sl_lib_failure,
        sl_censoring = sl_lib_censor,
        sl_treatment = sl_lib_g,
        cv.Control = list(V = 10))
    saveRDS(sl_fit, here("R/final_cvot/03_LEADER/sl_fit.RDS"), compress = F)
}

haz_sl <- list(sl_fit$density_failure_1$clone(),
               sl_fit$density_failure_0$clone())
haz_sl[[1]]$haz2surv()
haz_sl[[2]]$haz2surv()
names(haz_sl) <- c("A = 1", "A = 0")

results$estimates <- data.frame(s1 = colMeans(haz_sl[[1]]$survival),
                                s0 = colMeans(haz_sl[[2]]$survival)) %>%
    mutate(t = 1:length(s1)) %>% filter(t %in% targets) %>%
    cbind(Estimator = "G-Comp", .) %>%
    bind_rows(results$estimates, .)


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
                     glm.trt = glm_trt,
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

# plot results ------------------------------------------------------------

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
# saveRDS(object = result_tbl, here("R/final_cvot/03_LEADER/LEADER-estimates-tbl.RDS"))
result_tbl <- readRDS(here("R/final_cvot/03_LEADER/LEADER-estimates-tbl.RDS"))

result_plot <- result_tbl %>% filter(Estimand == "RR") %>%
    ggplot(aes(x = as.character(t), y = Estimate, colour = Estimator)) +
    geom_errorbar(aes(ymin = Estimate - 1.96*se,
                      ymax = Estimate + 1.96*se), width = .5,
                  position = position_dodge(.5)) +
    geom_point(position = position_dodge(.5)) + theme_minimal() +
    labs(x = "Months", y = "Relative Risk", title = "LEADER re-Analysis") # +
# geom_label(aes(x =`t`, y = Estimate - 1.96*`se` - .02,
#                label = paste0("Eff = ", round(Eff, 3)*100, "%")),
#            data = filter(result_tbl, Estimand  == "RR", Estimator == "TMLE"),
#            colour = 'black')
# ggsave(filename = "leader.png", path = "R/final_cvot/03_LEADER/",
#        plot = result_plot, device = "png", width = 9, height = 6, units = "in")

result_plot <- result_tbl %>% filter(Estimand == "RD", Estimator != "G-Comp") %>%
    ggplot(aes(x = as.character(`t`), y = Estimate, colour = Estimator)) +
    geom_errorbar(aes(ymin = Estimate - 1.96*se,
                      ymax = Estimate + 1.96*se), width = .5,
                  position = position_dodge(.5)) +
    geom_point(position = position_dodge(.5)) + theme_minimal() +
    labs(x = "\nTime (Months)", y = "Risk Difference\n",
         title = "LEADER Analysis: Composite Event (MACE)") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=16)) +
geom_label(aes(x = as.character(`t`), y = Estimate - 1.96 * `se` - 0.001,
               label = paste0("Eff = ", round(Eff, 3)*100, "%")),
           data = filter(result_tbl, Estimand  == "RD", Estimator == "TMLE"),
           colour = 'black')
result_plot

ggsave(filename = "leader_MACE.png", path = "R/competing_risks/LEADER/",
       plot = result_plot, device = "png", width = 9, height = 6, units = "in")





# cox RR estimate -----------------------------------------------------------------------------

tmpobs <- mutate(obs, ARM = factor(ARM))

ldr_cox <- survfit(coxph(Surv(obs$TIME, obs$EVENT) ~ ARM, data = tmpobs),
                   newdata = data.frame("ARM" = unique(tmpobs$ARM)))

ldr_cox_rr <- summary(ldr_cox, times = 1:4 * 12)$surv %>% as_tibble() %>%
    mutate("t" = as.character(1:4 * 12),
           Estimator = "G-Comp",
           RR = (1 - `2`) / (1 - `1`))

res_plot <- result_tbl %>% filter(Estimand  == "RR")
res_plot[res_plot$Estimator == "G-Comp", "Estimate"] <- ldr_cox_rr$RR
res_plot <- res_plot %>%
    ggplot(aes(x=`t`, y = Estimate, colour = Estimator)) +
    geom_errorbar(aes(ymin = Estimate - 1.96*se,
                      ymax = Estimate + 1.96*se), width = .5,
                  position = position_dodge(.5)) +
    geom_point(position = position_dodge(.5)) + theme_minimal() +
    labs(x = "Months", y = "Relative Risk", title = "LEADER re-Analysis")

# ggsave(filename = "leader.png", path = "R/final_cvot/03_LEADER/",
#        device = "png", width = 9, height = 6, units = "in", result_plot)
