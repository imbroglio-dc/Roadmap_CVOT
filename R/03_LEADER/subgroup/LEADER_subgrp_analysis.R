
# setup & data ------------------------------------------------------------


library(tidyverse); library(mice); library(survminer); library(survival)
library(SuperLearner); library(survtmle); library(MOSS); library(here)
library(foreach); library(doParallel); library(doRNG)
source(file = here("R/functions/my_sl_functions.R"))
source(file = here("R/functions/my_MOSS_hazard_methods.R"))
source(file = here("R/functions/my_survtmle_functions.R"))

set.seed(123456789)
n_cores <- 16

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

# subgroup analysis -------------------------------------------------------

subgroups <- obs %>%
    transmute(Sex = as.factor(case_when(SEX == "M" ~ "Male",
                                        SEX == "F" ~ "Female")),
              Age = factor(case_when(AGE.imp == T ~ ">=60 yr",
                                     AGE < 60 ~ "<60 yr",
                                     AGE >= 60 ~ ">=60 yr")),
              Region = as_factor(
                  case_when(ISREGION == "Western Europe" ~ "Europe",
                            ISREGION == "Eastern Europe" ~ "Europe",
                            T ~ as.character(ISREGION))),
              Race = as_factor(RACE),
              EthnicGroup = as_factor(
                  case_when(ETHNIC == "HISPANIC OR LATINO" ~ "Hispanic",
                            ETHNIC == "NOT HISPANIC OR LATINO" ~ "Non-Hispanic")),
              BMI = factor(case_when(BMIBL.imp == T ~ "NA",
                                     BMIBL > 30 ~ ">30",
                                     BMIBL <= 30 ~ "<=30"),
                           exclude = "NA"),
              HbA1c = as_factor(case_when(HBA1CBL <= 8.3 ~ "<=8.3%",
                                          HBA1CBL > 8.3 ~ ">8.3%")),
              diabDur = factor(case_when(DIABDUR.imp == T ~ "NA",
                                         DIABDUR > 11 ~ ">11 yrs",
                                         DIABDUR <= 11 ~ "<=11 yrs"),
                               exclude = "NA"),
              CVDRisk = as_factor(
                  case_when(CV.HIRISK == T ~ ">=50 yr + established CVD",
                            CV.HIRISK == F ~ ">=60 yr + risk factors for CVD")),
              CHF = as_factor(
                  case_when(CHFFL == T ~ "Yes",
                            CHFFL == F ~ "No")),
              antiDiabTherapy = factor(ANTDBFL,
                                       levels = c("1 OAD", "> 1 OADs", "Insulin+OAD(s)",
                                                  "Insulin-OAD", "None"),
                                       labels = c("1 Oral antidiabetic agent",
                                                  ">1 Oral antidiabetic agent",
                                                  "Insulin with oral antidiabetic agent",
                                                  "Insulin without oral antidiabetic agent",
                                                  "None")),
              renDis_SevMod = factor(
                  case_when(RENFSEV.MDRD.BL %in% c("Severe (EGFR<30)",
                                                   "Moderate (EGFR<60)") ~ "<60 ml/min/1.73m^2",
                            T ~ ">=60 ml/min/1.73m^2")),
              renDis_Sev = factor(
                  case_when(RENFSEV.MDRD.BL %in% c("Severe (EGFR<30)") ~ "<30 ml/min/1.73m^2",
                            T ~ ">=30 ml/min/1.73m^2"))) %>%
    as_tibble()

# superlearner libraries

leader_cens_glm <- function(Y, X, newX, family, obsWeights, model = TRUE, ...)
{ # indicators for > 14, 15, 16 are for timescale = 3, i.e. 3month time intervals
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
{ # indicators for > 14, 15, 16 are for timescale = 3, i.e. 3month time intervals
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

leader_cens_glmnet <- function (Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10,
                                nlambda = 100, useMin = TRUE, loss = "deviance", ...)
{
    .SL.require("glmnet")
    if (!is.matrix(X)) {
        X <- model.matrix(~-1 + ., X)
        newX <- model.matrix(~-1 + ., newX)
    }

    X <- cbind(X, "t11" = as.numeric(X[, "t"] == 11))
    X <- cbind(X, "t12" = as.numeric(X[, "t"] == 12))

    newX <- cbind(newX, "t11" = as.numeric(newX[, "t"] == 11))
    newX <- cbind(newX, "t12" = as.numeric(newX[, "t"] == 12))

    fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
                               lambda = NULL, type.measure = loss, nfolds = nfolds,
                               family = family$family, alpha = alpha, nlambda = nlambda,
                               ...)
    pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin,
                                                                      "lambda.min", "lambda.1se"))
    fit <- list(object = fitCV, useMin = useMin)
    class(fit) <- "SL.glmnet"
    out <- list(pred = pred, fit = fit)
    return(out)
}
environment(leader_cens_glmnet) <- asNamespace("SuperLearner")

screen_cor <- create.Learner(base_learner = "screen.corRank",
                             tune = list("rank" = c(15, 40, 60)),
                             detailed_names = T)

glmnets <- create.Learner(base_learner = "SL.glmnet",
                          tune = list("alpha" = c(0, 0.5, 1)),
                          detailed_names = T)

cens_glmnets <- create.Learner(base_learner = "leader_cens_glmnet",
                               tune = list("alpha" = c(0, 0.5, 1)),
                               detailed_names = T)

# estimation ----------------------------------------------------------------------------------

if (!file.exists(here(paste0("R/03_LEADER/subgroup/output/",
                             "LEADER-subgrp-estimates-fullcov.RDS")))) {

    iter <- do.call(rbind, lapply(1:ncol(subgroups), function(j) {
        levels <- levels(dplyr::pull(subgroups, colnames(subgroups)[j]))
        index <- t(sapply(1:length(levels), function(i)
            c("subgroup" = colnames(subgroups)[j], "level" = levels[i])))
        return(index)
    }))[c(1:4, 8, 5, 7, 6, 12, 10:9, 11, 14, 13, 15:16, 18, 17,
          19:20, 22, 21, 24, 23, 25:33), ] %>%
        as_tibble()

    subgroup_results <- vector("list", length = nrow(iter))
    names(subgroup_results) <- rownames(iter)

    B <- 1
    while (B <= nrow(iter)) {
        indices <- seq(B, min(nrow(iter), B + n_cores - 1))
        B <- B + n_cores
        cat("Estimating Subgroup Levels", min(indices), "-",
            max(indices), "\n")
        subgroup_results[indices] <- foreach(k = indices,
                                             .inorder = TRUE) %dopar%
            {
                subgroup <- dplyr::pull(subgroups, iter$subgroup[k])
                level <- iter$level[k]

                subgroup_obs <- obs[subgroup %in% level, ] %>%
                    select_if(~length(unique(.)) > 1)

                subgroup_W <- dplyr::select(subgroup_obs, -c(TIME, EVENT, ARM))
                subgroup_W <- as_tibble(model.matrix(~-1 + ., subgroup_W))
                subgroup_W_names <- colnames(subgroup_W)
                subgroup_W <- rename_all(subgroup_W, ~paste0("L", 1:ncol(subgroup_W)))
                names(subgroup_W_names) <- colnames(subgroup_W)

                results <- list("estimates" = NULL,
                                "subgroup_W" = subgroup_W_names)

                ## Cox --------------------------------------------------------------
                cox_fit <- coxph(Surv(subgroup_obs$TIME,
                                      subgroup_obs$EVENT) ~ subgroup_obs$ARM)
                results$hr <- tibble(hr = exp(cox_fit$coefficients),
                                     se = sqrt(diag(cox_fit$var)),
                                     CI = paste0(round(exp(cox_fit$coefficients
                                                           - 1.96*se), 4), " - ",
                                                 round(exp(cox_fit$coefficients +
                                                               1.96*se), 4)),
                                     CI_width = abs(diff(exp(cox_fit$coefficients +
                                                                 c(-1, 1)*1.96*se))))

                ## kaplan-meier ------------------------------------------------------
                km_fit <- survival::survfit(survival::Surv(subgroup_obs$TIME,
                                                           subgroup_obs$EVENT,
                                                           type = "right") ~ ARM,
                                            type = "kaplan-meier",
                                            data = subgroup_obs)
                km_est <- summary(km_fit, times = targets)

                results$estimates <- tibble(time = round(km_est$time, 1)) %>%
                    mutate(A = rep(c(0, 1), each = length(time)/2),
                           s = km_est$surv, # s0, s1
                           se = km_est$std.err) %>%
                    pivot_wider(names_from = "A",
                                values_from = c('s', 'se'), names_sep = "") %>%
                    cbind(Estimator = "Kaplan-Meier", .)


                # TMLE setup ----------------------------------------------------------
                if (nrow(subgroup_obs) <= 500) {
                    screeners <- screen_cor$names[1] # 15 vars
                } else if (nrow(subgroup_obs) <= 1000) {
                    screeners <- screen_cor$names[2] # 40 vars
                } else if (nrow(subgroup_obs) <= 1500) {
                    screeners <- screen_cor$names[3] # 60 vars
                } else {
                    screeners <- "All" # 40, 60, or all vars
                }

                sl_trt <- expand.grid(c("SL.glm", glmnets$names, "SL.bayesglm"), "All")
                sl_trt <- lapply(1:nrow(sl_trt),
                                 function(i) as.character(unlist(sl_trt[i, ])))
                sl_censor <- expand.grid(c("SL.glm", cens_glmnets$names, "SL.bayesglm",
                                           "leader_cens_glm", "leader_cens_bayesglm"),
                                         screeners)
                sl_censor <- lapply(1:nrow(sl_censor),
                                    function(i) as.character(unlist(sl_censor[i, ])))
                sl_events <- expand.grid(c("SL.glm", "SL.xgboost", "SL.rpart",
                                           glmnets$names, "SL.bayesglm"), screeners)
                sl_events <- lapply(1:nrow(sl_events),
                                    function(i) as.character(unlist(sl_events[i, ])))


                ## TMLE: initial fit --------------------------------------------
                suppressWarnings(suppressMessages(
                    sl_fit <- my_init_sl_fit(
                        T_tilde = subgroup_obs$TIME,
                        Delta = as.numeric(subgroup_obs$EVENT),
                        A = as.numeric(subgroup_obs$ARM),
                        W = subgroup_W,
                        t_max = max(targets),
                        sl_failure = sl_events,
                        sl_censoring = sl_censor,
                        sl_treatment = "SL.glm",
                        cv.Control = list(V = 10))
                ))

                haz_sl <- list(sl_fit$density_failure_1$clone(),
                               sl_fit$density_failure_0$clone())
                haz_sl[[1]]$haz2surv()
                haz_sl[[2]]$haz2surv()
                names(haz_sl) <- c("A = 1", "A = 0")

                results$estimates <- data.frame(s1 = colMeans(haz_sl[[1]]$survival),
                                                s0 = colMeans(haz_sl[[2]]$survival)) %>%
                    mutate(time = 1:length(s1)) %>% filter(time %in% targets) %>%
                    cbind(Estimator = "G-Comp", .) %>%
                    bind_rows(results$estimates, .)

                results$sl_fit_Y <- sl_fit$models$Y
                results$sl_fit_C <- sl_fit$models$C

                sl_Y <- sl_fit$models$Y
                sl_GdC <- sl_fit$G_dC
                saveRDS(object = sl_fit,
                        file = paste0("R/03_LEADER/subgroup/sl_fits/",
                                              iter$subgroup[k], "-",
                                              k, "_fits.RDS"))

                # TMLE: update and output ----------------------------------------
                suppressWarnings(suppressMessages(
                    tmle_sl <- surv_tmle(ftime = subgroup_obs$TIME,
                                         ftype = as.numeric(subgroup_obs$EVENT),
                                         targets = targets,
                                         trt = as.numeric(subgroup_obs$ARM),
                                         t0 = max(targets), adjustVars = subgroup_W,
                                         SL.ftime = sl_Y, SL.ctime = sl_GdC,
                                         SL.trt = sl_trt,
                                         returnIC = T, returnModels = F,
                                         ftypeOfInterest = 1, trtOfInterest = c(1, 0),
                                         maxIter = 30, method = "hazard")
                ))

                tmle_sl_out <- suppressWarnings(
                    t(1 - tmle_sl$est) %>% unname() %>% as_tibble() %>%
                        cbind(Estimator = "TMLE", time = targets, .) %>%
                        rename("s0" = V1, "s1" = V2))
                results$estimates <- suppressWarnings(
                    matrix(sqrt(diag(tmle_sl$var)),
                           nrow = length(targets), byrow = F,
                           dimnames = list(NULL, c("se0", "se1"))) %>%
                        as.data.frame() %>% cbind(tmle_sl_out, .)) %>%
                    bind_rows(results$estimates, .)

                ic <- t(tmle_sl$ic)

                results$se <- try(
                    tibble(s0_se = apply(ic[1:(nrow(ic)/2), ], 1, var),
                           s1_se = apply(ic[(nrow(ic)/2 + 1):nrow(ic), ], 1, var),
                           RD_se = sapply(1:(nrow(ic)/2), function(i) {
                               var(colSums(c(1, -1) * ic[c(i, i+(nrow(ic)/2)), ]))
                           }),
                           RR_se = sapply(1:(nrow(ic)/2), function(i) {
                               s0 = tmle_sl_out[i, "s0"]
                               s1 = tmle_sl_out[i, "s1"]
                               var(colSums(c((1-s1)/(1-s0)^2, -1/(1-s0)) *
                                               ic[c(i, i+(nrow(ic)/2)), ]))
                           }),
                           SR_se = sapply(1:(nrow(ic)/2), function(i) {
                               s0 = tmle_sl_out[i, "s0"]
                               s1 = tmle_sl_out[i, "s1"]
                               var(colSums(c(-(s1)/(s0)^2, 1/s0) *
                                               ic[c(i, i+(nrow(ic)/2)), ]))
                           })) %>% apply(., 2,
                                         function(j) sqrt(j / ncol(ic))) %>%
                        as_tibble() %>% bind_cols(t = targets, .) %>%
                        pivot_longer(-`t`, names_to = c("Estimand", "tmp"),
                                     names_sep = "_",
                                     values_to = "se") %>%
                        dplyr::select(-tmp) %>%
                        cbind(Estimator = "TMLE", .))

                results$tmle_fit <- tmle_sl

                results$estimates <- cbind(results$estimates,
                                           "n" = nrow(subgroup_obs))

                # Output results ----------------------------------------------
                return(results)
            }
    # }
    saveRDS(subgroup_results,
            file = here(paste0("R/03_LEADER/subgroup/output/",
                               "LEADER-subgrp-estimates-fullcov.RDS")),
            compress = F)
} else {
    subgroup_results <-
        read_rds(here(paste0("R/03_LEADER/subgroup/output/",
                             "LEADER-subgrp-estimates-fullcov.RDS")))
}


# subgroup plots -------------------------------------------------------------------

if (!file.exists(here("R/03_LEADER/subgroup/output/LEADER-subgrp-est-tbl.csv"))) {
    result_tbl <- lapply(1:length(subgroup_results), function(k) {
        subgroup <- dplyr::pull(subgroups, colnames(subgroups)[iter[k, 1]])
        subgrp <- paste0(colnames(subgroups)[iter[k, 1]], ": ",
                         unique(subgroup)[iter[k, 2]])
        cbind("subgroup" = subgrp, subgroup_results[[k]]$estimates)
    }) %>% bind_rows() %>%
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
        filter(se_est == Estimand) %>% dplyr::select(-se_est) %>% distinct() %>%
        pivot_wider(names_from = Estimator, names_sep = "_",
                    values_from = c(se, Estimate)) %>%
        mutate(Eff = `se_Kaplan-Meier`^2 / `se_TMLE`^2,
               rel_CI = `se_TMLE` / `se_Kaplan-Meier`) %>%
        pivot_longer(-c(subgroup, `time`, Estimand, Eff, rel_CI, n), names_sep = "_",
                     names_to = c(".value", "Estimator")) %>%
        mutate(Eff = case_when(Estimator == "TMLE" ~ Eff,
                               Estimator == "G-Comp" ~ NaN,
                               T ~ 1),
               rel_CI =  case_when(Estimator == "TMLE" ~ rel_CI,
                                   Estimator == "G-Comp" ~ NaN,
                                   T ~ 1)) %>%
        dplyr::select(subgroup, Estimator, `time`, Estimand, Estimate, se, Eff, rel_CI, n) %>%
        mutate(`time` = as.character(`time` * timescale), CI_width = 2*1.96*se)
    write_csv(result_tbl, here("R/03_LEADER/subgroup/output/LEADER-subgrp-est-tbl.csv"))
} else {
    result_tbl <- read_csv(here("R/03_LEADER/subgroup/output/LEADER-subgrp-est-tbl.csv"))
}

result_plot <- result_tbl %>% filter(Estimand  == "RR") %>%
    ggplot(aes(y=as.character(`time`), x = Estimate, colour = Estimator)) +
    facet_wrap(~subgroup, scales = "free") +
    geom_errorbar(aes(xmin = Estimate - 1.96*se,
                      xmax = Estimate + 1.96*se), width = .5,
                  position = position_dodge(.5)) +
    geom_vline(xintercept = 1, alpha = 0.7) +
    geom_point(position = position_dodge(.5)) + theme_bw() +
    labs(x = "Months", y = "Relative Risk", title = "LEADER re-Analysis") +
    # geom_label(aes(y = as.character(`time`), x = Estimate - 1.96*`se` - .02,
    #                label = paste0("Eff = ", round(Eff, 2)*100, "%")),
    #            data = filter(result_tbl, Estimand  == "RR", Estimator == "TMLE"),
    #            colour = 'black')

ggsave(filename = "leader.png", path = "R/03_LEADER/subgroup/output",
       device = "png", width = 9, height = 6, units = "in", result_plot)

# Year 4 RR  ----------------------------------------------------------------------------------

tmp <- result_tbl %>% filter(Estimand == "RR", time == 12*timescale) %>%
    left_join(., bind_rows(lapply(1:ncol(subgroups), function(i) {
        subgroups[, i] %>% group_by_all %>% tally() %>%
            rename_all(~c("subgroup", "n")) %>%
            mutate(subgroup = as_factor(paste0(colnames(subgroups)[i], ": ", subgroup)))
    }))) %>%
    mutate(CI = paste0(round(Estimate, 2), " (",
                       round(Estimate - 1.96*se, 2), " - ",
                       round(Estimate + 1.96*se, 2), ")")) %>%
    dplyr::select(subgroup, Estimator, Estimate, CI, se, n) %>%
    group_by(subgroup) %>% mutate(rel_CI = se / head(se, 1))

# hypothesis tests
rr_pvals <- lapply(c("TMLE", "Kaplan-Meier"), function(estimator) {
    tmp_est <- filter(tmp, Estimator == estimator) %>% dplyr::select(-Estimator)
    out <- lapply(colnames(subgroups), function(sg) {
        data <- filter(tmp_est, str_detect(subgroup, paste0(sg, ":")))
        if (nrow(data) == 2) {
            tstat <- diff(data$Estimate) / sqrt(data$se[1]^2 + data$se[2]^2)
            df <- (data$se[1]^2 + data$se[2]^2)^2 / (data$se[1]^4/(data$n[1] - 1) +
                                                         data$se[2]^4 / (data$n[2] - 1))
            pvalue <- 2*pt(abs(tstat), df, lower.tail = FALSE)
            names(pvalue) <- sg
            return(pvalue)
        }
        else {
            k <- nrow(data)
            w_j <- 1 / data$se^2
            w <- sum(w_j)
            x_bar <- sum(data$Estimate * w_j) / w
            num <- 1 / (k - 1) * sum(w_j * (data$Estimate - x_bar)^2)
            denom <- 1 + (2 * (k - 2)) / (k^2 - 1) * sum(1 / (data$n - 1) * (1 - w_j / w)^2)
            W <- num / denom
            df = c(k - 1, (k^2 - 1) / (3 * sum((1 - w_j / w)^2 / (data$n - 1))))
            pvalue <- pf(W, df[1], df[2], lower.tail = FALSE)
            names(pvalue) <- sg
            return(pvalue)
        }
    }) %>% unlist() %>% data.frame(.)
    colnames(out) <- estimator
    return(out)
}) %>% bind_cols()



# hazard ratios -------------------------------------------------------------------------------

hr_tbl <- lapply(subgroup_results, function(sg) sg$hr) %>% bind_rows() %>%
    cbind("subgroup" = names(subgroup_results), .) %>%
    left_join(., bind_rows(lapply(1:ncol(subgroups), function(i) {
        subgroups[, i] %>% group_by_all %>% tally() %>%
            rename_all(~c("subgroup", "n")) %>%
            mutate(subgroup = as_factor(paste0(colnames(subgroups)[i], ": ", subgroup)))
    }))) %>% dplyr::select("subgroup", everything())

# hypothesis tests
hr_pvals <- lapply(colnames(subgroups), function(sg) {
    data <- filter(hr_tbl, str_detect(subgroup, paste0(sg, ":")))
    if (nrow(data) == 2) {
        tstat <- diff(log(data$hr)) / sqrt(data$se[1]^2 + data$se[2]^2)
        df <- (data$se[1]^2 + data$se[2]^2)^2 / (data$se[1]^4/(data$n[1] - 1) +
                                                     data$se[2]^4 / (data$n[2] - 1))
        pvalue <- 2*pt(abs(tstat), df, lower.tail = FALSE)
        names(pvalue) <- sg
        return(pvalue)
    }
    else {
        k <- nrow(data)
        w_j <- 1 / data$se^2
        w <- sum(w_j)
        x_bar <- sum(log(data$hr) * w_j) / w
        num <- 1 / (k - 1) * sum(w_j * (log(data$hr) - x_bar)^2)
        denom <- 1 + (2 * (k - 2)) / (k^2 - 1) * sum(1 / (data$n - 1) * (1 - w_j / w)^2)
        W <- num / denom
        df = c(k - 1, (k^2 - 1) / (3 * sum((1 - w_j / w)^2 / (data$n - 1))))
        pvalue <- pf(W, df[1], df[2], lower.tail = FALSE)
        names(pvalue) <- sg
        return(pvalue)
    }
}) %>% unlist() %>% data.frame("p_val" = .)


# compare relative CI -----------------------------------------------------

rel_ci <- tmp %>% filter(Estimator != "G-Comp") %>%
    dplyr::select(subgroup, Estimator, Estimate, se) %>%
    mutate(Estimate = paste0(round(Estimate, 2), ": (",
                             round(Estimate - 1.96*se, 4), " - ",
                             round(Estimate + 1.96*se, 4), ")"),
           se = se*1.96*2) %>%
    rename("ci" = se) %>%
    pivot_wider(names_from = Estimator, values_from = c(Estimate, ci)) %>% left_join(hr_tbl) %>%
    mutate(hr = paste0(round(hr, 2), ": (", CI, ")")) %>%
    dplyr::select(subgroup, `Estimate_Kaplan-Meier`, Estimate_TMLE,
                  `ci_Kaplan-Meier`, ci_TMLE, hr, CI_width) %>%
    rename("ci_HR" = CI_width, Estimate_HR = hr) %>%
    pivot_longer(cols = c(`ci_Kaplan-Meier`, ci_TMLE, ci_HR),
                 names_to = "Estimator", values_to = "ci", names_prefix = "ci_") %>%
    group_by(subgroup) %>% mutate(rel_CI = ci / tail(ci, 1)) %>%
    pivot_longer(cols = c(`Estimate_Kaplan-Meier`, Estimate_TMLE, Estimate_HR),
                 names_to = "estimator", values_to = "Estimate", names_prefix = "Estimate_") %>%
    filter(estimator == Estimator) %>% dplyr::select(subgroup, Estimator, Estimate, rel_CI)


filter(result_tbl,
       str_detect(subgroup, "Race"), time == 48,
       Estimand == "RR", Estimator == "TMLE")$se %>%
    sapply(., function(se) pnorm(1 - 1.96*se, 0.87, se, lower.tail = T))



