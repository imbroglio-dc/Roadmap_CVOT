
# packages ----------------------------------------------------------------


library(survival); library(survminer);
library(MOSS); library(SuperLearner); library(survtmle)
library(parallel); library(foreach); library(doParallel)
library(tidyverse); library(skimr); library(here)
source(file = here("R/functions/my_sl_functions.R"))
source(file = here("R/functions/my_MOSS_hazard_methods.R"))
source(file = here("R/functions/my_survtmle_functions.R"))
options(na.action = na.omit)
set.seed(0)


# test_leader -------------------------------------------------------------


test_leader <- readxl::read_excel(here("data/test_leader.xlsx"))
test_leader <- test_leader %>%
    mutate_if(is_character, as_factor) %>%
    mutate_if(~length(levels(.)) == 2, ~as.logical(as.numeric(.)-1)) %>%
    mutate(ARM = as.numeric(ARM)) %>%
    select(ARM, event, time_days, everything(), -subjid) %>%
    # filter_all(all_vars(!is.na(.))) # remove 9 rows missing BMIBL
    mutate(BMI_missing = is.na(BMIBL),
           BMIBL = case_when(is.na(BMIBL) ~ mean(BMIBL, na.rm = T),
                             T ~ BMIBL)) # impute missing BMIBL w/ mean

cens_km <- survfit(Surv(time = test_leader$time_days, event = 1 - test_leader$event,
                        type = 'right') ~ 1, type = "kaplan-meier")
cens_km <- stepfun(cens_km$time / 7, c(1, cens_km$surv))



# Simulating Outcomes -----------------------------


set.seed(0)
sim_leader <- test_leader %>%
    mutate(ARM = rbinom(nrow(test_leader), 1, 0.5),
           AGE = as.numeric(scale(AGE)),
           SMOKER = case_when(SMOKER == "CURRENT SMOKER" ~ 1,
                              SMOKER == "PREVIOUS SMOKER" ~ .5,
                              T ~ 0)) %>%
    sample_n(1000)

time_span <- seq(1, 300, 1)
d_cens <- survival_curve$new(t = time_span,
                             survival = cens_km(time_span))
d_cens$surv2haz()
d_cens$hazard[is.na(d_cens$hazard)] <- 0
d_cens$hazard <- matrix(pmax(.001, as.numeric(d_cens$hazard)),
                        nrow = 6, ncol = ncol(d_cens$hazard), byrow = T)
d_cens$hazard[, 1:48] <- (diag(c(3^c(0, .5, 1), 16^c(0, .5, 1))) %*%
    d_cens$hazard)[, 1:48]
d_cens$haz2surv()
cens_surv <- list("0" = d_cens$survival[1:3, ],
                  "1" = d_cens$survival[4:6, ])
cens_surv <- lapply(cens_surv, function(s) {
    rownames(s) <- sort(unique(sim_leader$SMOKER))
    return(s)
    })

cens_survivals <- lapply(cens_surv, function(s) {
    s <- count(dplyr::select(sim_leader, SMOKER), SMOKER)$n * s
    s <- colSums(s) / nrow(sim_leader)
})

as.data.frame(d_cens$survival) %>%
    mutate(A = rep(c("0", "1"), each = 3),
           W = rep(c("0", "0.5", "1"), times = 2)) %>%
    pivot_longer(-c(A, W), names_to = "t", names_prefix = "V", values_to = "s") %>%
    mutate("t" = as.numeric(`t`)) %>% ggplot() +
    geom_line(aes(x = `t`, y = `s`, colour = W, linetype = A)) + theme_minimal() +
    labs(title = "Censoring `Survival` Curves")



# Sim Hazards and Survivals ---------------------------------------


sim_haz_fn <- function(t, A, W, lambda = c(1.5e-3), betas = c(16, 1.5, 1)) {
    t1 <- log(betas[1]) * W # multiplies base hazard by beta^W
    t2 <- log(betas[2]) * A * W # multiplies base hazard by beta^(A*W)
    t3 <- log(betas[3]) * A * t # dummy var. necessary to return vector along t
    haz_AW <- lambda * exp(t1 + t2 + t3)
    return(haz_AW)
}


hazards <- lapply(0:1, function(arm) {
    h <- do.call(rbind, lapply(c(0, 0.5, 1), function(w) {
        sim_haz_fn(t = time_span, A = arm, W = w)
    }))
    rownames(h) <- c(0, 0.5, 1)
    return(h)
})
names(hazards) <- c("0", "1")

do.call(rbind, lapply(c("1", "0"), function(a) {
    cbind(t = time_span, A = a, as.data.frame(t(hazards[[a]])))
})) %>% pivot_longer(cols = c(`0`, `0.5`, `1`), names_to = "W", values_to = "haz") %>%
    ggplot(aes(x=t, y=haz, colour=W, linetype=A)) + geom_line() + theme_minimal() +
    labs(title = "Conditional Hazards")

## CONDITIONAL SURVIVALS | A, W
s_AW <- lapply(hazards, function(haz) {
    s_aw <- t(apply(haz, 1, function(haz_w) {
        tail(cumprod(1 - c(0, haz_w)), -1)
    }))
    colnames(s_aw) <- 1:ncol(s_aw)
    return(s_aw)
})

do.call(rbind, lapply(c("1", "0"), function(a) {
    cbind(t = time_span, A = a, as.data.frame(t(s_AW[[a]])))
})) %>% pivot_longer(cols = c(`0`, `0.5`, `1`), names_to = "W", values_to = "haz") %>%
    ggplot(aes(x=t, y=haz, colour=W, linetype=A)) + geom_line() + theme_minimal() +
    labs(title = "Conditional Survivals")


## TREATMENT-SPECIFIC SURVIVAL CURVES
survivals <- lapply(s_AW, function(s) {
    s <- count(dplyr::select(sim_leader, SMOKER), SMOKER)$n * s
    s <- colSums(s) / nrow(sim_leader)
})

ic_plot <- tibble(t = rep(time_span, 2),
       A = rep(names(survivals), each = length(time_span)),
       Failure = as.vector(do.call(c, survivals)),
       Censoring = as.vector(do.call(c, cens_survivals))) %>%
    pivot_longer(c(Failure, Censoring), names_to = "Type", values_to = "y") %>%
    mutate(Type = factor(Type, levels = c("Failure", "Censoring"))) %>%
    ggplot(aes(x = t, y = y, colour = A, linetype = Type)) +
    geom_line(size = 1.2) + theme_minimal() +
    geom_vline(aes(xintercept = 104), alpha = .5, size = 1.2,
               colour = "blue", linetype = 3) +
    labs(title = "Data-Generating Curves for Informative Censoring",
         y = expr("Counterfactual Survival\n"), x= "Time") +
    theme(title = element_text(size = 16),
          plot.title = element_text(size = 16, face = "bold"),
          axis.title = element_text(face = "bold", size = 16),
          axis.text = element_text(size = 16),
          legend.text = element_text(size = 14))
ggsave(ic_plot, filename = "ic-data.png", device = "png",
       path = here("R/02_informative-censoring-sim/"),
       width = 11, height = 6, units = "in")


# TRUE RISKS
targets <- 2*52 # 3 year by weekly interval
true_risks <- tibble(t = rep(time_span, 2),
                     A = rep(names(survivals), each = length(time_span)),
                     s = as.vector(do.call(c, survivals))) %>%
    filter(t %in% targets) %>% group_by(t) %>%
    pivot_wider(names_from = A, values_from = s, names_prefix = "s") %>%
    mutate(RD = s0 - s1, RR = (1 - s1) / (1 - s0))
true_risks

truths <- as.data.frame(t(mutate(true_risks[, -1], SR = s1/s0))) %>%
    rownames_to_column() %>% rename(`stats` = rowname, "truth" = V1)

# simulating data ----------------------------------------------------------

B <- 1e3
n_cores <- max(detectCores() - 4, 1)
registerDoParallel(n_cores)
cl <- makeForkCluster(n_cores)
clusterSetRNGStream(cl = cl, iseed = 0)
sim_data <- foreach(b = 1:B,
                    .errorhandling = "remove",
                    .verbose = T) %dopar%
    {
        sim <- sample_n(tbl = sim_leader, size = nrow(sim_leader), replace = T) %>%
            mutate(ARM = rbinom(nrow(sim_leader), 1, .5))
        sim_w <- as.matrix(summarise(group_by(sim, ARM, SMOKER), n = n()))
        sim <- bind_rows(apply(sim_w, 1, function(r) {
            s_aw <- s_AW[[as.character(r['ARM'])]][as.character(r['SMOKER']), ]
            c_aw <- cens_surv[[as.character(r['ARM'])]][as.character(r['SMOKER']), ]
            f_time <- colSums(outer(s_aw, runif(r['n'], 0, 1), `>`)) + 1
            c_time <- colSums(outer(c_aw, runif(r['n'], 0, 1), `>`)) + 1
            time <- ceiling(pmin(f_time, c_time))
            event <- as.numeric(f_time < c_time)
            tibble("time" = time, "event" = event,
                   "ARM" = r["ARM"], "SMOKER" = r["SMOKER"])
        })) %>%
            cbind(dplyr::select(sim, -c("ARM", "SMOKER", "time_days", "event")))
        # can bind properly to respect covariate covariance

        # summary(coxph(Surv(time = sim$time,
        #                    event = sim$event,
        #                    type = 'right') ~ ARM, data = sim))
        # plot(survfit(Surv(sim$time, sim$event, type = "right") ~ ARM,
        #              type = "kaplan-meier", data = sim), ymin = .2)
        return(sim)
    }
stopCluster(cl)
rm(cl)
saveRDS(sim_data, file = here("R/02_informative-censoring-sim/inf-cens-data.RDS"))
# sim_data <- read_rds(here("R/02_informative-censoring-sim/inf-cens-data.RDS"))

# Estimation --------------------------------------------------------------

n_cores <- max(detectCores() - 4, 1)
registerDoParallel(n_cores)
cl <- makeForkCluster(n_cores)
B <- length(sim_data)
clusterSetRNGStream(cl = cl, iseed = 0)
system.time(
    sim_estimates <- foreach(sim = sim_data[1:B], #
                             .errorhandling = "pass",
                             .verbose = T) %dopar%
        {
            options(warn = -1)
            results <- list()

            # kaplan-meier ------------------------------------------------------------

            km_fit <- try(summary(survfit(Surv(sim$time, sim$event, type = "right") ~
                                              ARM, type = "kaplan-meier", data = sim),
                                  times = targets, scale = 1))

            results$KM <- try(tibble(t = round(km_fit$time, 1)) %>%
                                  mutate(A = rep(c(0, 1), each = length(t)/2),
                                         s = km_fit$surv, # s0, s1
                                         se = km_fit$std.err) %>%
                                  pivot_wider(names_from = A, values_from = c(s, se),
                                              names_sep = "") %>%
                                  mutate("RR" = (1-s1)/(1-s0), RD = s0 - s1) %>%
                                  dplyr::select(t, s0, s1, RD, RR, se0, se1) %>%
                                  cbind(Estimator = "Kaplan-Meier", .))


            # initial haz fit ---------------------------------------------------

            sl_lib_g <- "SL.glm"
            sl_lib_censor <- c("SL.mean", "cens_glm", "SL.glmnet",
                               "sl_bayesglm", "SL.ranger")
            sl_lib_failure <- c("SL.mean", "SL.glm", "SL.glmnet",
                                "SL.bayesglm", "SL.xgboost", "SL.ranger")
            timescale <- 4
            sim <- mutate(sim, time = ceiling(time/timescale))
            W <- dplyr::select(sim, SMOKER, AGE, BMIBL)

            # misspec haz -------------------------------------------------------------

            suppressMessages(suppressWarnings(
                sl_fit_misspec <- my_init_sl_fit( # allows passing in cvControl args
                    T_tilde = sim$time, # fortnights
                    Delta = as.numeric(sim$event),
                    A = as.numeric(sim$ARM),
                    W = W,
                    t_max = max(sim$time),
                    sl_failure = "simple_misspec",
                    sl_censoring = "SL.glm",
                    sl_treatment = sl_lib_g,
                    cv.Control = list(V = 2),
                    cl = NULL) # need to debug initial sl fit parallelization
            ))

            haz_mis <- list(sl_fit_misspec$density_failure_1$clone(deep = TRUE),
                            sl_fit_misspec$density_failure_0$clone(deep = TRUE))
            # haz_mis <- lapply(haz_mis, upscale_haz, timescale = timescale)
            haz_mis[[1]]$haz2surv()
            haz_mis[[2]]$haz2surv()
            names(haz_mis) <- c("A = 1", "A = 0")

            # g-comp misspec ----------------------------------------------

            k_grid <- (1:ncol(haz_mis[[1]]$hazard) - 1) * timescale
            surv_init <- data.frame(s = c(colMeans(haz_mis[[1]]$survival),
                                          colMeans(haz_mis[[2]]$survival)),
                                    A = rep(c("A = 1", "A = 0"), each = length(k_grid)),
                                    t = k_grid)

            results$gcomp <- try(surv_init %>%
                                     filter(t %in% targets) %>% group_by(t) %>%
                                     pivot_wider(names_from = A, values_from = s) %>%
                                     rename(s1 = `A = 1`, s0 = `A = 0`) %>%
                                     mutate(RD = s0 - s1, RR = (1 - s1) / (1 - s0)) %>%
                                     ungroup() %>% dplyr::select(t, s0, s1, RD, RR) %>%
                                     cbind(Estimator = "G-Comp: Misspecified", .))

            # correct haz --------------------------------------------------------------

            suppressMessages(suppressWarnings(
                sl_fit_correct <- my_init_sl_fit( # allows passing in cvControl args
                    T_tilde = sim$time, # fortnights
                    Delta = as.numeric(sim$event),
                    A = as.numeric(sim$ARM),
                    W = W,
                    t_max = max(sim$time),
                    sl_failure = "correct_infcens",
                    sl_censoring = "SL.glm",
                    sl_treatment = sl_lib_g,
                    cv.Control = list(V = 2),
                    cl = NULL) # need to debug initial sl fit parallelization
            ))

            haz_cor <- list(sl_fit_correct$density_failure_1$clone(deep = TRUE),
                            sl_fit_correct$density_failure_0$clone(deep = TRUE))
            # haz_cor <- lapply(haz_cor, upscale_haz, timescale = timescale)
            haz_cor[[1]]$haz2surv()
            haz_cor[[2]]$haz2surv()
            names(haz_cor) <- c("A = 1", "A = 0")

            # g-comp correct ----------------------------------------------

            surv_init <- data.frame(s = c(colMeans(haz_cor[[1]]$survival),
                                          colMeans(haz_cor[[2]]$survival)),
                                    A = rep(c("A = 1", "A = 0"), each = length(k_grid)),
                                    t = k_grid)

            results$gcomp <- try(surv_init %>%
                                     filter(t %in% targets) %>% group_by(t) %>%
                                     pivot_wider(names_from = A, values_from = s) %>%
                                     rename(s1 = `A = 1`, s0 = `A = 0`) %>%
                                     mutate(RD = s0 - s1, RR = (1 - s1) / (1 - s0)) %>%
                                     ungroup() %>% dplyr::select(t, s0, s1, RD, RR) %>%
                                     cbind(Estimator = "G-Comp: Correct", .) %>%
                                     rbind(results$gcomp, .))

            # sl haz ------------------------------------------------------------------

            suppressMessages(suppressWarnings(
                sl_fit <- my_init_sl_fit( # allows passing in cvControl args
                    T_tilde = sim$time, # fortnights
                    Delta = as.numeric(sim$event),
                    A = as.numeric(sim$ARM),
                    W = W,
                    t_max = max(sim$time),
                    sl_failure = sl_lib_failure,
                    sl_censoring = sl_lib_censor,
                    sl_treatment = sl_lib_g,
                    cv.Control = list(V = 5),
                    cl = NULL) # need to debug initial sl fit parallelization
            ))

            haz_sl <- list(sl_fit$density_failure_1$clone(deep = TRUE),
                           sl_fit$density_failure_0$clone(deep = TRUE))
            # haz_sl <- lapply(haz_sl, upscale_haz, timescale = timescale)
            haz_sl[[1]]$haz2surv()
            haz_sl[[2]]$haz2surv()
            names(haz_sl) <- c("A = 1", "A = 0")

            # g-comp sl ----------------------------------------------

            surv_init <- data.frame(s = c(colMeans(haz_sl[[1]]$survival),
                                          colMeans(haz_sl[[2]]$survival)),
                                    A = rep(c("A = 1", "A = 0"), each = length(k_grid)),
                                    t = k_grid)

            results$gcomp <- try(surv_init %>%
                                     filter(t %in% targets) %>% group_by(t) %>%
                                     pivot_wider(names_from = A, values_from = s) %>%
                                     rename(s1 = `A = 1`, s0 = `A = 0`) %>%
                                     mutate(RD = s0 - s1, RR = (1 - s1) / (1 - s0)) %>%
                                     ungroup() %>% dplyr::select(t, s0, s1, RD, RR) %>%
                                     cbind(Estimator = "G-Comp: SuperLearner", .) %>%
                                     rbind(results$gcomp, .))

            # SurvTMLE Haz ---------------------------------------------------------

            # SurvTMLE Haz ---------------------------------------------------------

            ## from g-comp hazard
            results$SL_hazards <- list("misspec" = haz_mis,
                                       "correct" = haz_cor,
                                       "sl" = haz_sl)

            SL_ftime <- list(sl_fit_misspec, sl_fit_correct, sl_fit)
            names(SL_ftime) <- c("misspec", "correct", "sl")
            SL_ftime <- lapply(SL_ftime,
                               function(sl_fits) {
                                   haz_fit <- do.call(rbind, lapply(1:2, function(i) {
                                       h <- t(sl_fits[[i]]$hazard)
                                       colnames(h) <- 1:ncol(h)
                                       h <- cbind(t = 1:nrow(h), trt = 2-i, h)
                                   })) %>% as_tibble() %>%
                                       pivot_longer(cols = -c(t, trt),
                                                    names_to = "id",
                                                    values_to = "Q1PseudoHaz") %>%
                                       mutate(id = as.numeric(id)) %>%
                                       arrange(id, trt, t)
                                   list("J1" = haz_fit)
                               })
            results$SL_ftime <- SL_ftime

            sl_G_dC <- sl_fit$G_dC

            glm_trt <- paste0(colnames(W), collapse = " + ")
            glm_ftime_mis <- "trt + SMOKER"
            glm_ftime_cor <- "SMOKER + SMOKER:trt"


            # tmle misspec ------------------------------------------------------------

            tmle_mis <- try(surv_tmle(ftime = sim$time, ftype = sim$event,
                                      targets = targets/timescale, trt = sim$ARM,
                                      t0 = max(targets/timescale), adjustVars = W,
                                      SL.ftime = SL_ftime[["misspec"]],
                                      SL.ctime = sl_G_dC, glm.trt = glm_trt,
                                      returnIC = T, returnModels = T,
                                      ftypeOfInterest = 1, trtOfInterest = c(1, 0),
                                      maxIter = 10, method = "hazard"))

            results$survtmle_h <- suppressWarnings(
                try(t(1 - tmle_mis$est) %>% unname() %>%
                        cbind(t = targets, .) %>% as_tibble() %>%
                        rename("s0" = V2, "s1" = V3) %>%
                        mutate(RD = s0 - s1, RR = (1-s1)/(1-s0))))
            results$survtmle_h <- suppressWarnings(
                try(matrix(sqrt(diag(tmle_mis$var)),
                           nrow = length(targets), byrow = T,
                           dimnames = list(NULL, c("se0", "se1"))) %>%
                        as.data.frame() %>% cbind(results$survtmle_h, .) %>%
                        cbind(Estimator = "TMLE: Hazard - Misspec", .)))
            results$tmleICmis <- try(tmle_mis$ic)

            # tmle correct ------------------------------------------------------------

            tmle_cor <- try(surv_tmle(ftime = sim$time, ftype = sim$event,
                                      targets = targets/timescale, trt = sim$ARM,
                                      t0 = max(targets/timescale), adjustVars = W,
                                      SL.ftime = SL_ftime[["correct"]],
                                      glm.trt = glm_trt, SL.ctime = sl_G_dC,
                                      returnIC = T, returnModels = T,
                                      ftypeOfInterest = 1, trtOfInterest = c(1, 0),
                                      maxIter = 10, method = "hazard"))

            tmle_cor_out <- suppressWarnings(
                try(t(1 - tmle_cor$est) %>% unname() %>%
                        cbind(t = targets, .) %>% as_tibble() %>%
                        rename("s0" = V2, "s1" = V3) %>%
                        mutate(RD = s0 - s1, RR = (1-s1)/(1-s0))))
            results$survtmle_h <- suppressWarnings(
                try(matrix(sqrt(diag(tmle_cor$var)),
                           nrow = length(targets), byrow = T,
                           dimnames = list(NULL, c("se0", "se1"))) %>%
                        as.data.frame() %>% cbind(tmle_cor_out, .) %>%
                        cbind(Estimator = "TMLE: Hazard - Correct", .) %>%
                        rbind(results$survtmle_h, .)))
            results$tmleICcor <- try(tmle_cor$ic)

            # tmle sl -----------------------------------------------------------------

            tmle_sl <- try(surv_tmle(ftime = sim$time, ftype = sim$event,
                                     targets = targets/timescale, trt = sim$ARM,
                                     t0 = max(targets/timescale), adjustVars = W,
                                     SL.ftime = SL_ftime[["sl"]],
                                     glm.trt = glm_trt, SL.ctime = sl_G_dC,
                                     returnIC = T, returnModels = T,
                                     ftypeOfInterest = 1, trtOfInterest = c(1, 0),
                                     maxIter = 10, method = "hazard"))

            tmle_sl_out <- suppressWarnings(
                try(t(1 - tmle_sl$est) %>% unname() %>%
                        cbind(t = targets, .) %>% as_tibble() %>%
                        rename("s0" = V2, "s1" = V3) %>%
                        mutate(RD = s0 - s1, RR = (1-s1)/(1-s0))))
            results$survtmle_h <- suppressWarnings(
                try(matrix(sqrt(diag(tmle_sl$var)),
                           nrow = length(targets), byrow = T,
                           dimnames = list(NULL, c("se0", "se1"))) %>%
                        as.data.frame() %>% cbind(tmle_sl_out, .) %>%
                        cbind(Estimator = "TMLE: Hazard - SL", .) %>%
                        rbind(results$survtmle_h, .)))
            results$tmleICsl <- try(tmle_sl$ic)

            # results ----------------------------------------------------------------

            return(results)
        })
stopCluster(cl)
rm(cl)

saveRDS(sim_estimates, here("R/02_informative-censoring-sim/inf-cens-estimates.RDS"))
# sim_estimates <- read_rds(here("R/02_informative-censoring-sim/inf-cens-estimates.RDS"))


# estimator performance ----------------------------------------------------------

estimates <- bind_rows(lapply(sim_estimates, function(sim_est) {
    sim_est$KM <- sim_est$KM %>%
        mutate(s0_se = se0, s1_se = se1,
               RR_se = sqrt(se1^2 / (1-s0)^2 + se0^2 * ((1-s1) / (1-s0)^2)^2),
               RD_se = sqrt(se1^2 + se0^2),
               SR_se = sqrt(se1^2 / s0^2 + se0^2 * s1^2 / s0^4))
    ic_ind <- str_detect(names(sim_est), "IC")
    est_ind <- str_detect(names(sim_est), "survtmle")
    ic <- t(do.call(cbind, sim_est[ic_ind]))
    rownames(ic) <- paste0(rep(names(sim_est)[ic_ind], each = 2),
                           c(".0", ".1"))
    est <- do.call(rbind, sim_est[est_ind]) %>%
        dplyr::select(Estimator, s0, s1)

    se <- tibble(s0_se = apply(ic[seq(1, nrow(ic), 2), ], 1, var),
                 s1_se = apply(ic[seq(2, nrow(ic), 2), ], 1, var),
                 RD_se = sapply(1:(nrow(ic)/2), function(i) {
                     var(colSums(c(1, -1) * ic[(2*i-1):(2*i), ]))
                 }),
                 RR_se = sapply(1:(nrow(ic)/2), function(i) {
                     s0 = est[i, "s0"] ; s1 = est[i, "s1"]
                     var(colSums(c((1-s1)/(1-s0)^2, -1/(1-s0)) *
                                     ic[(2*i-1):(2*i), ]))
                 }),
                 SR_se = sapply(1:(nrow(ic)/2), function(i) {
                     s0 = est[i, "s0"] ; s1 = est[i, "s1"]
                     var(colSums(c(-(s1)/(s0)^2, 1/s0) * ic[(2*i-1):(2*i), ]))
                 })) %>% apply(., 2, function(j) sqrt(j / ncol(ic))) %>%
        as_tibble() %>% bind_cols(Estimator = est$Estimator, .)
    sim_est$survtmle_h <- suppressWarnings(
        left_join(sim_est$survtmle_h, se, by = "Estimator"))
    suppressWarnings(intersect(names(sim_est),
                               c("KM", "gcomp", "survtmle_h", "survtmle_m")) %>%
                         sim_est[.] %>% bind_rows() %>%
                         dplyr::select(-c(se0, se1)) %>%
                         mutate(SR = s1/s0,
                                Estimator = factor(Estimator,
                                                   levels = c("Kaplan-Meier",
                                                              "G-Comp: Misspecified",
                                                              "G-Comp: Correct",
                                                              "G-Comp: SuperLearner",
                                                              "TMLE: Hazard - Misspec",
                                                              "TMLE: Hazard - Correct",
                                                              "TMLE: Hazard - SL"),
                                                   labels = c("Kaplan-Meier", "G-Comp: Mis",
                                                              "G-Comp: Corr", "G-Comp: SL",
                                                              "TMLE: Mis", "TMLE: Corr",
                                                              "TMLE: SL"))) %>%
                         separate(col = Estimator, into = c("Estimator", "Hazard"), sep = ": ") %>%
                         mutate(Hazard = case_when(Hazard == "Mis" ~ "Misspecified",
                                                   Hazard == "Corr"  ~ "Correct",
                                                   Hazard == "SL" ~ "SuperLearner")))
})) %>%
    mutate(iter = rep(1:length(sim_estimates),
                      each = n()/length(sim_estimates)))
estimates <- estimates %>%
    filter(Estimator == "Kaplan-Meier") %>% dplyr::select(-Hazard) %>%
    cbind(Hazard = rep(c("Misspecified", "Correct", "SuperLearner"),
                       each = length(sim_estimates))) %>%
    dplyr::select(Estimator, Hazard, everything()) %>%
    rbind(., filter(estimates, Estimator != "Kaplan-Meier")) %>%
    mutate(Estimator = factor(Estimator, levels = c("Kaplan-Meier", "G-Comp",
                                                    "TMLE", "LTMLE"))) %>%
    arrange(iter, Estimator) %>% as_tibble()

estimates <- dplyr::select(true_risks, c("t", "RR", "s0", "s1")) %>%
    rename(RR_true = RR, s0_true = s0, s1_true = s1) %>%
    mutate(RD_true = s0_true - s1_true, SR_true = s1_true / s0_true) %>%
    right_join(estimates, ., by = 't') %>%
    pivot_longer(cols = -c(t, Estimator, Hazard, iter, contains(c("se", "true"))),
                 names_to = "estimand", values_to = "estimate") %>%
    pivot_longer(cols = contains(c("true")), names_to = "estimand_t",
                 names_pattern = "([[:alnum:]]+)", values_to = "truth") %>%
    filter(estimand == estimand_t) %>%
    pivot_longer(cols = contains(c("se")), names_to = "estimand_se",
                 names_pattern = "([[:alnum:]]+)", values_to = "se") %>%
    filter(estimand_t == estimand_se) %>%
    dplyr::select(-c(estimand_t, estimand_se))

coverage <- estimates %>%
    mutate(Bias = estimate - truth,
           MSE = (estimate - truth)^2,
           Cover = abs(estimate-truth) < 1.96*se) %>%
    group_by(Estimator, Hazard, estimand, t) %>%
    mutate(O.se = sqrt(var(estimate))) %>%
    dplyr::select(-c(iter, estimate, truth, se)) %>%
    summarise_all(mean) %>% ungroup  # %>%
# mutate(Hazard = case_when(as.character(Estimator) == "Kaplan-Meier" ~ "NA",
#                           T ~ as.character(Hazard))) %>% distinct()
saveRDS(sim_estimates, here("R/02_informative-censoring-sim/inf-cens-coverage.RDS"))

# plots -------------------------------------------------------------------

plot_df <- estimates %>%
    left_join(., coverage, by = c("Estimator", "Hazard", "t", "estimand")) %>%
    mutate(estimand = factor(estimand,
                             levels = c("s0", "s1", "RD", "RR", "SR"),
                             labels = c("Control Survival", "Treated Survival",
                                        "Risk Difference", "Relative Risk",
                                        "Relative Survival")))
saveRDS(plot_df, here("R/02_informative-censoring-sim/inf-cens-plot_df.RDS"))

plot_df %>% ggplot() +
    geom_boxplot(aes(y = estimate, x = Estimator), outlier.shape = NA) +
    facet_wrap(Hazard~estimand, ncol = 5, scales = "free") + theme_minimal() +
    geom_hline(aes(yintercept = truth),
               data = distinct(dplyr::select(plot_df, estimand, truth)),
               colour = "blue", alpha = .7) +
    labs(title = "Estimator Performance with Informative Censoring") +
    geom_label(aes(x = Estimator, y = estimate,
                   label = if_else(is.na(Cover), "NA",
                                   paste0(round(Cover, 2)*100, "%"))),
               data = summarise_all(group_by_if(plot_df, ~!is.numeric(.)),
                                    ~quantile(., .95, na.rm = T)),
               size = 2.5, na.rm = T)
ggsave(filename = "inf-cens-perf-plot.png",
       path = here("R/02_informative-censoring-sim/"),
       device = "png", width = 16, height = 9, units = "in")


# cox hazard ratio --------------------------------------------------------

hr_est <- lapply(sim_data, function(sim) {
    cox <- coxph(Surv(sim$time, sim$event) ~ sim$ARM)
    return(exp(cox$coefficients))
}) %>% unlist() %>% as_tibble() %>%
    summarise(est = mean(value), se = sqrt(var(value)),
              lower = quantile(value, 0.025), upper = quantile(value, 0.975))

