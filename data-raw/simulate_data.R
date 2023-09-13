# test_leader -------------------------------------------------------------
library(tidyverse)
set.seed(0)

load("../data/test_leader.rda")
cens_km <- survival::survfit(survival::Surv(time = test_leader$time_days,
                                            event = 1 - test_leader$event,
                                            type = 'right') ~ 1,
                             type = "kaplan-meier")
cens_km <- stats::stepfun(cens_km$time / 7, c(1, cens_km$surv)) # convert to weekly


sim_leader <- test_leader %>%
    mutate(ARM = rbinom(n(), 1, 0.5),
           AGE = as.numeric(scale(AGE)),
           SMOKER = case_when(SMOKER == "CURRENT SMOKER" ~ 1,
                              SMOKER == "PREVIOUS SMOKER" ~ .5,
                              T ~ 0),
           OBESEBL = as.numeric(BMIBL > 30))


# Constant Hazard ---------------------------------------

set.seed(0)
sim_W <- group_by(sim_leader, SMOKER, OBESEBL) %>% summarise(n = n()) %>%
    arrange(SMOKER, OBESEBL)

time_span <- seq(1, 300, 1)
sim_haz_fn <- function(t, A, W, lambda = c(2e-3),
                       betas = c(1, 9, 1.5, .5, 2)) {
    t1 <- log(betas[1]) * A * t # time-varying term, multiples hazard by 1 along t
    t2 <- log(betas[2]) * W["SMOKER"] # multiplies base hazard by beta^W
    t3 <- log(betas[3]) * A * W["SMOKER"] # multiplies base hazard by beta^(A*W)
    t4 <- log(betas[4]) * (A == 1) * W["OBESEBL"]
    t5 <- log(betas[5]) * (A == 0) * W["OBESEBL"]
    haz_AW <- lambda * exp(t1 + t2 + t3 + t4 + t5)
    return(haz_AW)
}

hazards <- lapply(0:1, function(arm) {
    h <- t(apply(sim_W, 1, function(w) {
        sim_haz_fn(t = time_span, A = arm, W = w)
    }))
    rownames(h) <- apply(sim_W, 1, function(w)
        paste0("Sm:", w["SMOKER"], ", Ob:", w["OBESEBL"]))
    return(h)
})
names(hazards) <- c("0", "1")

s_AW <- lapply(hazards, function(haz) {
    s_aw <- t(apply(haz, 1, function(haz_w) {
        tail(cumprod(1 - c(0, haz_w)), -1)
    }))
    colnames(s_aw) <- 1:ncol(s_aw)
    return(s_aw)
})

constant_hazard <- sim_leader %>% group_by(ARM, SMOKER, OBESEBL) %>%
    summarise(N = n()) %>% as.matrix()
constant_hazard <- bind_rows(apply(constant_hazard, 1, function(r) {
    s_aw = s_AW[[as.character(r['ARM'])]][
        paste0("Sm:", r["SMOKER"], ", Ob:", r["OBESEBL"]), ]
    f_time = colSums(outer(s_aw, runif(r['N'], 0, 1), `>`)) + 1
    c_time = colSums(outer(cens_km(time_span), runif(r['N'], 0, 1), `>`)) + 1
    time = ceiling(pmin(f_time, c_time))
    event = as.numeric(f_time < c_time)
    tibble("time" = time, "event" = event, "ARM" = r["ARM"],
           "SMOKER" = r["SMOKER"], "OBESEBL" = r["OBESEBL"])
})) %>% sample_frac() %>%
    cbind(., dplyr::select(sim_leader, -c("ARM", "SMOKER", "OBESEBL",
                                          "time_days", "event")))

attr(constant_hazard, "true_survivals") <- lapply(s_AW, function(s) {
    colSums(s * sim_W$n) / nrow(sim_leader)})


# Informative Censoring ---------------------------------------------------

# censoring based on observed, scaled by 3^SMOKER for ctl and 16^SMOKER for txt
set.seed(0)
d_cens <- MOSS::survival_curve$new(
    t = 1:ceiling(max(test_leader$time_days)/7),
    survival = cens_km(1:ceiling(max(test_leader$time_days) / 7))
    )
d_cens$survival_to_hazard()
d_cens$hazard <- matrix(pmax(.001, as.numeric(d_cens$hazard)),
                        nrow = 6, ncol = ncol(d_cens$hazard), byrow = T)
d_cens$hazard[, 1:48] <- (diag(c(3^c(0, .5, 1), 16^c(0, .5, 1))) %*%
                              d_cens$hazard)[, 1:48]
d_cens$hazard_to_survival()
cens_surv <- list("0" = d_cens$survival[1:3, ],
                  "1" = d_cens$survival[4:6, ])
cens_surv <- lapply(cens_surv, function(s) {
    rownames(s) <- sort(unique(sim_leader$SMOKER))
    return(s)
})

sim_W <- group_by(sim_leader, SMOKER) %>% summarise(n = n()) %>% arrange(SMOKER)

sim_haz_fn <- function(t, A, W, lambda = c(1.5e-3), betas = c(16, 1.5, 1)) {
    t1 <- log(betas[1]) * W["SMOKER"] # base hazard x beta^W
    t2 <- log(betas[2]) * A * W["SMOKER"] # base hazard x beta^(A*W)
    t3 <- log(betas[3]) * A * t # dummy var. necessary to return vector along t
    haz_AW <- lambda * exp(t1 + t2 + t3)
    return(haz_AW)
}


hazards <- lapply(0:1, function(arm) {
    h <- t(apply(sim_W, 1, function(w) {
        sim_haz_fn(t = time_span, A = arm, W = w)
    }))
    rownames(h) <- apply(sim_W, 1, function(w) w["SMOKER"])
    return(h)
})
names(hazards) <- c("0", "1")

## CONDITIONAL SURVIVALS | A, W
s_AW <- lapply(hazards, function(haz) {
    s_aw <- t(apply(haz, 1, function(haz_w) {
        tail(cumprod(1 - c(0, haz_w)), -1)
    }))
    colnames(s_aw) <- 1:ncol(s_aw)
    return(s_aw)
})

informative_censoring <- sim_leader %>%
    mutate(ARM = rbinom(nrow(sim_leader), 1, .5)) %>%
    group_by(ARM, SMOKER) %>% summarise(N = n()) %>% as.matrix()

informative_censoring <- bind_rows(apply(informative_censoring, 1, function(r) {
    s_aw <- s_AW[[as.character(r['ARM'])]][as.character(r['SMOKER']), ]
    c_aw <- cens_surv[[as.character(r['ARM'])]][as.character(r['SMOKER']), ]
    f_time <- colSums(outer(s_aw, runif(r['N'], 0, 1), `>`)) + 1
    c_time <- colSums(outer(c_aw, runif(r['N'], 0, 1), `>`)) + 1
    time <- ceiling(pmin(f_time, c_time))
    event <- as.numeric(f_time < c_time)
    tibble("time" = time, "event" = event,
           "ARM" = r["ARM"], "SMOKER" = r["SMOKER"])
})) %>% sample_frac() %>%
    cbind(dplyr::select(sim_leader, -c("ARM", "SMOKER", "time_days", "event")))

attr(informative_censoring, "true_survivals") <- lapply(s_AW, function(s) {
    colSums(s * sim_W$n) / nrow(sim_leader)})
