

# setup & data ------------------------------------------------------------


library(tidyverse); library(mice); library(survminer); library(survival)
library(SuperLearner); library(survtmle); library(MOSS); library(here)
library(foreach); library(doParallel)
source(file = here("R/functions/my_sl_functions.R"))
source(file = here("R/functions/my_MOSS_hazard_methods.R"))
source(file = here("R/functions/my_survtmle_functions.R"))

set.seed(0)
n_cores <- 6
registerDoParallel(n_cores)
source(file = here("R/functions/LEADER_W_clean.R"))
timescale <- 4 # time coarsened to 4 month intervals
times <- 1:4 * 12 / timescale

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
    mice::mice(m = 1, maxit = 10) %>% mice::complete() %>%
    cbind(., dplyr::select(obs, c(TIME, EVENT, ARM, INSNVFL, CHFFL))) %>%
    dplyr::select(TIME, EVENT, ARM, everything())

obs <- W %>% select_if(anyNA) %>% rename_all(~paste0(., ".imp")) %>%
    mutate_all(is.na) %>% cbind(obs, .)

# subgroup analysis -------------------------------------------------------

obs_subgrp <- obs %>%
    mutate(RACE = case_when(RACE == "WHITE" ~ "WHITE",
                            RACE == "BLACK OR AFRICAN AMERICAN" ~ "BLACK",
                            RACE == 'ASIAN' ~ "ASIAN",
                            T ~ "OTHER"))
# dplyr::select(obs, "TIME", "EVENT", "ARM", "AGE", "MALE", "EGFREPB",
#                         "CREATBL", "BMIBL", "HBA1CBL", "HISPANIC", "MIFL",
#                         "SMOKER", "STROKEFL", "BMIBL.imp", "AGE.imp",
#                         "RACE", "ANTDBFL", "CVRISK") %>% as_tibble()

subgroups <- transmute(obs,
                       Sex = as.factor(case_when(MALE ~ "Male",
                                                 T ~ "Female")),
                       Age = as.factor(case_when(AGE >= 60 ~ ">=60 yr",
                                                 T ~ "<60 yr")),
                       Region = as_factor(
                           case_when(ISREGION == "Western Europe" ~ "Europe",
                                     ISREGION == "Eastern Europe" ~ "Europe",
                                     T ~ as.character(ISREGION))),
                       Race = as_factor(
                           case_when(RACE == "WHITE" ~ "WHITE",
                                     RACE == "BLACK OR AFRICAN AMERICAN" ~ "BLACK",
                                     RACE == 'ASIAN' ~ "ASIAN",
                                     T ~ "OTHER")),
                       EthnicGroup = as_factor(case_when(HISPANIC ~ "Hispanic",
                                                         T ~ "Non Hispanic")),
                       BMI = as_factor(case_when(BMIBL <= 30 ~ "<=30",
                                                 T ~ ">30")),
                       HbA1c = as.factor(case_when(HBA1CBL > 8.3 ~ ">8.3",
                                                   T ~ "<=8.3")),
                       CVDRisk = as.factor(
                           case_when(CVRISK ~ ">=60 yr + risk factors for CVD",
                                     T ~ ">=50 yr + established CVD")),
                       CHF = CHFFL,
                       antiDiabTherapy = ANTDBFL,
                       'diabDur' = as.factor(case_when(DIABDUR > 11 ~ ">11 yrs",
                                                       T ~ "<=11 yrs")),
                       renDis_SevMod = (RENFSEV %in% c("Moderate (EGFR<60)",
                                                       "Severe (EGFR<30)")),
                       renDis_SevMod = as.factor(case_when(renDis_SevMod ~ "<60",
                                                           T ~ ">=60")),
                       renDis_Sev = RENFSEV == "Severe (EGFR<30)",
                       renDis_Sev = as.factor(case_when(renDis_Sev ~ "<30",
                                                        T ~ ">=30"))) %>%
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
                             tune = list("rank" = c(15, 40)),
                             detailed_names = T)

glmnets <- create.Learner(base_learner = "SL.glmnet",
                          tune = list("alpha" = c(0, 0.5, 1)),
                          detailed_names = T)

cens_glmnets <- create.Learner(base_learner = "leader_cens_glmnet",
                               tune = list("alpha" = c(0, 0.5, 1)),
                               detailed_names = T)

# estimation ----------------------------------------------------------------------------------

if (!file.exists(here("R/final_cvot/03_LEADER/LEADER-subgrp-estimates.RDS"))) {

    iter <- do.call(rbind, lapply(1:ncol(subgroups), function(j) {
        levels <- unique(dplyr::pull(subgroups, colnames(subgroups)[j]))
        index <- t(sapply(1:length(levels), function(i) c(j, i)))
        rownames(index) <- paste0(colnames(subgroups)[j], ": ", levels)
        return(index)
    }))
    colnames(iter) <- c("j", "i")

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
                j <- iter[k, 1]
                i <- iter[k, 2]
                subgroup <- dplyr::pull(subgroups, colnames(subgroups)[j])
                levels <- unique(subgroup)

                subgroup_obs <- obs_subgrp[subgroup == levels[i], ] %>%
                    select_if(~length(unique(.)) > 1)

                subgroup_W <- dplyr::select(subgroup_obs, -c(TIME, EVENT, ARM))
                subgroup_W <- model.matrix(~-1 + ., subgroup_W) %>%
                    as_tibble()
                subgroup_W_names <- colnames(subgroup_W)
                subgroup_W <- subgroup_W %>% rename_all(~paste0("L", 1:ncol(subgroup_W)))
                names(subgroup_W_names) <- colnames(subgroup_W)

                results <- list("estimates" = NULL,
                                "subgroup_W" = subgroup_W_names)

                ## Cox --------------------------------------------------------------
                cox_fit <- coxph(Surv(subgroup_obs$TIME,
                                      subgroup_obs$EVENT) ~ subgroup_obs$ARM)
                results$hr <- tibble(hr = exp(cox_fit$coefficients),
                                     se = sqrt(diag(cox_fit$var)),
                                     CI = paste0(round(exp(cox_fit$coefficients
                                                           - 1.96*se), 2), " - ",
                                                 round(exp(cox_fit$coefficients +
                                                               1.96*se), 2)),
                                     CI_width = abs(diff(exp(cox_fit$coefficients +
                                                                 c(-1, 1)*1.96*se))))

                ## kaplan-meier ------------------------------------------------------
                km_fit <- survival::survfit(survival::Surv(subgroup_obs$TIME,
                                                           subgroup_obs$EVENT,
                                                           type = "right") ~ ARM,
                                            type = "kaplan-meier",
                                            data = subgroup_obs)
                km_est <- summary(km_fit, times = times)

                results$estimates <- tibble(time = round(km_est$time, 1)) %>%
                    mutate(A = rep(c(0, 1), each = length(time)/2),
                           s = km_est$surv, # s0, s1
                           se = km_est$std.err) %>%
                    pivot_wider(names_from = "A",
                                values_from = c('s', 'se'), names_sep = "") %>%
                    cbind(Estimator = "Kaplan-Meier", .)


                # TMLE setup ----------------------------------------------------------
                screeners <- c(screen_cor$names, "All") #, screen.rf$names)
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
                        t_max = max(times),
                        sl_failure = sl_events,
                        sl_censoring = sl_censor,
                        sl_treatment = "SL.mean",
                        cv.Control = list(V = 10))
                ))

                haz_sl <- list(sl_fit$density_failure_1$clone(),
                               sl_fit$density_failure_0$clone())
                haz_sl[[1]]$haz2surv()
                haz_sl[[2]]$haz2surv()
                names(haz_sl) <- c("A = 1", "A = 0")

                results$estimates <- data.frame(s1 = colMeans(haz_sl[[1]]$survival),
                                                s0 = colMeans(haz_sl[[2]]$survival)) %>%
                    mutate(time = 1:length(s1)) %>% filter(time %in% times) %>%
                    cbind(Estimator = "G-Comp", .) %>%
                    bind_rows(results$estimates, .)

                SL_ftime <- sl_fit$models$Y
                sl_G_dC <- sl_fit$G_dC
                saveRDS(sl_fit, file = paste0("R/final_cvot/03_LEADER/subgroup/",
                                              colnames(subgroups)[j], "-",
                                              i, "_fits.RDS"))
                rm(sl_fit)


                # TMLE: update and output ----------------------------------------
                suppressWarnings(suppressMessages(
                    tmle_sl <- surv_tmle(ftime = subgroup_obs$TIME,
                                         ftype = as.numeric(subgroup_obs$EVENT),
                                         targets = times,
                                         trt = as.numeric(subgroup_obs$ARM),
                                         t0 = max(times), adjustVars = subgroup_W,
                                         SL.ftime = SL_ftime, SL.ctime = sl_G_dC,
                                         SL.trt = sl_trt,
                                         returnIC = T, returnModels = F,
                                         ftypeOfInterest = 1, trtOfInterest = c(1, 0),
                                         maxIter = 25, method = "hazard")
                ))

                tmle_sl_out <- suppressWarnings(
                    t(1 - tmle_sl$est) %>% unname() %>% as_tibble() %>%
                        cbind(Estimator = "TMLE", time = times, .) %>%
                        rename("s0" = V1, "s1" = V2))
                results$estimates <- suppressWarnings(
                    matrix(sqrt(diag(tmle_sl$var)),
                           nrow = length(times), byrow = F,
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
                        as_tibble() %>% bind_cols(t = times, .) %>%
                        pivot_longer(-`t`, names_to = c("Estimand", "tmp"),
                                     names_sep = "_",
                                     values_to = "se") %>%
                        dplyr::select(-tmp) %>%
                        cbind(Estimator = "TMLE", .))

                results$tmle_fit <- tmle_sl


                # Output results ----------------------------------------------
                return(results)
            }
    }
    saveRDS(subgroup_results, file = here("R/final_cvot/03_LEADER/LEADER-subgrp-estimates-fullcov.RDS"),
            compress = F)
} else {
    subgroup_results <- read_rds(here("R/final_cvot/03_LEADER/LEADER-subgrp-estimates.RDS"))
}


# subgroup plots -------------------------------------------------------------------

if (!file.exists(here("R/final_cvot/03_LEADER/subgroup/LEADER-subgrp-estimates.RDS"))) {
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
        pivot_longer(-c(subgroup, `time`, Estimand, Eff, rel_CI), names_sep = "_",
                     names_to = c(".value", "Estimator")) %>%
        mutate(Eff = case_when(Estimator == "TMLE" ~ Eff,
                               Estimator == "G-Comp" ~ NaN,
                               T ~ 1),
               rel_CI =  case_when(Estimator == "TMLE" ~ rel_CI,
                                   Estimator == "G-Comp" ~ NaN,
                                   T ~ 1)) %>%
        dplyr::select(subgroup, Estimator, `time`, Estimand, Estimate, se, Eff, rel_CI) %>%
        mutate(`time` = as.character(`time` * timescale), CI_width = 2*1.96*se)
    saveRDS(object = result_tbl, here("R/final_cvot/03_LEADER/subgroup/LEADER-subgrp-result_tbl.RDS"))
} else {
    result_tbl <- read_rds(here("R/final_cvot/03_LEADER/subgroup/LEADER-subgrp-result_tbl.RDS"))
}

result_plot <- result_tbl %>% filter(Estimand  == "RR") %>%
    ggplot(aes(y=as.character(`time`), x = Estimate, colour = Estimator)) +
    facet_wrap(~subgroup, scales = "free", ncol = ) +
    geom_errorbar(aes(xmin = Estimate - 1.96*se,
                      xmax = Estimate + 1.96*se), width = .5,
                  position = position_dodge(.5), ) +
    geom_vline(xintercept = 1, alpha = 0.7) +
    geom_point(position = position_dodge(.5)) + theme_bw() +
    labs(x = "Months", y = "Relative Risk", title = "LEADER re-Analysis") +
    geom_label(aes(x =as.character(`time`), y = Estimate - 1.96*`se` - .02,
                   label = paste0("Eff = ", round(Eff, 2)*100, "%")),
               data = filter(result_tbl, Estimand  == "RR", Estimator == "TMLE"),
               colour = 'black')

ggsave(filename = "leader.png", path = "R/final_cvot/03_LEADER/subgroup/",
       device = "png", width = 9, height = 6, units = "in", result_plot)

# Year 4 RR ----------------------------------------------------------------------------------

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

rel_ci <- tmp %>% filter(Estimator != "G-Comp") %>% dplyr::select(subgroup, Estimator, se) %>%
    mutate(se = se*1.96*2) %>% rename("ci" = se) %>%
    pivot_wider(names_from = Estimator, values_from = ci) %>% left_join(hr_tbl) %>%
    dplyr::select(subgroup, `Kaplan-Meier`, TMLE, CI_width) %>%
    rename("HR" = CI_width) %>% pivot_longer(cols = c(`Kaplan-Meier`, TMLE, HR),
                                             names_to = "Estimator", values_to = "ci") %>%
    group_by(subgroup) %>% mutate(rel_CI = ci / tail(ci, 1))

