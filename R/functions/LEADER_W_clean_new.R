
library(here); library(tidyverse); library(mice); library(haven)
data_path <- ("data/ex2211-3748-3/Analysis Ready Datasets/SAS_analysis/")
include_slices <- FALSE # includes slicing/binning of variables (e.g. AGE >= 65, BMI>30, etc)

# read in subject level covariate dataset
W <- read_sas(paste0(data_path, "adsl.sas7bdat"))

# 9341th subject in subject level dataset, not in the time-to-event dataset
#   FASFL = N, so supposed to be taken out?
# remove the subject who is not in the time-to-event dataset
W <- read_sas(paste0(data_path, "adtte.sas7bdat")) %>%
    dplyr::select(USUBJID) %>% distinct() %>%
    left_join(., W)

W <- W %>%
    # REMOVE UNITS, ID VARIABLES, AND POST-BASELINE VARIABLES ---------------
dplyr::select_at(vars(!ends_with("U"))) %>% # remove units
    dplyr::select_at(vars(!ends_with("N"))) %>% # redundant factor numbers
    dplyr::select(-c(ARMCD, ACTARM, ACTARMCD, TRTSEQP, TRTSEQA, TRT01P, TRT01A)) %>% # redundant ARMs
    dplyr::select_at(vars(!starts_with(c("TRDU", "EO", # treatment duration, end of study /
                                         "DCS", "TOTREAT")))) %>% # treatment, etc. (post-baseline)
    dplyr::select_if(~length(unique(.)) > 1) %>% # no variation (no information) in variable
    dplyr::select_at(vars(!ends_with("DT"))) %>%  # dates of ____
    dplyr::select_at(vars(!starts_with("RAND"))) %>%  # > 90% missingness
    dplyr::select(-c("age_category")) %>% # little variability (no information)
    dplyr::select(-c("COMPLFL")) %>% # post-baseline variable: completed trial flag
    dplyr::select(-c("SUBJID")) %>% # redundant w/ USUBJID
    dplyr::select(-c("CVHIFL")) %>% # redundant w/ CVRISK
    dplyr::select(-c("GERDCMFL")) %>% # post-baseline? what does conmed mean?
    dplyr::select(-c("HBC2BL")) %>% # redundant with HBA1CBL - just a unit conversion?
    dplyr::select(-c("HDL2BL")) %>% # redundant with HDL1BL - just a unit conversion?
    dplyr::select(-c("LDL2BL")) %>% # redundant with LDL1BL - just a unit conversion?
    dplyr::select(-c("LSTVISIT")) %>% # post-baseline
    dplyr::select(-c("NONCMPRS")) %>% # post-baseline
    dplyr::select(-c("PREVTFL")) %>% # Primary MACE Event Population Flag - no info, only empty or Y
    dplyr::select(-c("RACE")) %>% # redundant w/ RACEGR
    dplyr::select(-c("EXCFAIL", "INCFAIL", 'SCRNFL')) %>% # > 90% missingness
    dplyr::select(-c("EXMACEFL", 'NYHACLAS', 'RETINSEV')) %>% # > 80% missingness
    dplyr::select(-c("STDURY", 'TRHOLDUR', "TRPROP")) %>% # POST-BASELINE
    dplyr::select(-c("TRIG2BL")) %>%  # redundant with TRIG1BL - just a unit conversion?
    dplyr::select(-c("CHOL2BL")) %>%  # redundant with CHOL1BL - just a unit conversion?
    dplyr::select(-c("SITEID")) %>%
    # change column values / types --------------------------------------------
dplyr::mutate(CVRISK = CVRISK == "High", # ok, cvrisk no missingness
              # SITEID = as.factor(SITEID),
              SEX = as.factor(SEX),
              ARM = as.numeric(ARM == "Liraglutide"),
              SMOKER = factor(SMOKER, ordered = T,
                              levels = c("NEVER SMOKED", "PREVIOUS SMOKER",
                                         "CURRENT SMOKER")),
              RACEGR = as.factor(case_when(RACEGR == "BLACK OR AFRICAN AMERICAN" ~ "BLACK",
                                           T ~ RACEGR)),
              ANTDBFL = as.factor(ANTDBFL),
              AGEBLG2 = factor(AGEBLG2, ordered = T,
                               levels = c("Adults (below 65)", "Adults 65-74",
                                          "Adults 75-84", "Adults 85 or Older"),
                               labels = c("<65", "65-74", "75-84", ">=85"))) %>%
    mutate_all(~na_if(., "")) %>%
    mutate_at(vars(starts_with("RNF"), starts_with("RENF")),
              ~factor(., ordered = T,
                      levels = c("Normal (EGFR>=90)", "Mild (EGFR<90)",
                                 "Moderate (EGFR<60)", "Severe (EGFR<30)"))) %>%
    dplyr::mutate_if(~mean(unique(.) %in% c("Y", "N", "y", "n", "", NA)) == 1,
                     ~case_when(. %in% c("Y", "y") ~ TRUE,
                                . %in% c("N", "n") ~ FALSE,
                                T ~ NA)) %>%
    # mutate_if(is.character, as_factor) %>%
    # rename columns ----------------------------------------------------------
dplyr::rename(AGE.GRP = AGEBLG2,
              AGE.60PLUS = AGEGR1,
              BMI.30PLUS = BMIGRPBL,
              CV.HIRISK = CVRISK,
              DDUR.12PLUS = DDURGR,
              EGFR.CKD.BL = EGFREPB,
              EGFR.CKD.BL.30LESS = EGFC30BF,
              EGFR.CKD.BL.60LESS = EGFC60BF,
              EGFR.CKD.SCRN = EGFREPI,
              EGFR.CKD.SCRN.30LESS = EGFRC30F,
              EGFR.CKD.SCRN.60LESS = EGFRC60F,
              EGFR.MDRD.BL = EGFMDRBC,
              EGFR.MDRD.BL.30LESS = EGFMCBHF,
              EGFR.MDRD.BL.60LESS = EGFMCBLF,
              EGFR.MDRD.IVRS.SCRN = EGFRMDR,
              EGFR.MDRD.IVRS.SCRN.30LESS = EGFRM30F,
              EGFR.MDRD.IVRS.SCRN.60LESS = EGFRM60F,
              EGFR.MDRD.CALC.SCRN = EGFRMDRC,
              EGFR.MDRD.CALC.SCRN.30LESS = EGFRMCHF,
              EGFR.MDRD.CALC.SCRN.60LESS = EGFRMCLF,
              HBA1C.8.3LESS = HBA1CGR,
              RACE = RACEGR,
              RENFSEV.CKD.BL = RENFCKD,
              RENFSEV.MDRD.BL = RENFSEV,
              RENFSEV.CKD.SCRN = RNFSCCKD,
              RENFSEV.MDRD.SCRN = RNFSCSEV) %>% # AGE60PLUS, BMI30PLUS, DDUR12PLUS
    mutate(AGE.60PLUS = AGE >= 60,
           BMI.30PLUS = BMIBL > 30,
           DDUR.12PLUS = DIABDUR >= 12,
           EGFR.CKD.BL.30LESS = EGFR.CKD.BL < 30,
           EGFR.CKD.BL.60LESS = EGFR.CKD.BL < 60,
           EGFR.CKD.SCRN.30LESS = EGFR.CKD.SCRN < 30,
           EGFR.CKD.SCRN.60LESS = EGFR.CKD.SCRN < 60,
           EGFR.MDRD.BL.30LESS = EGFR.MDRD.BL < 30,
           EGFR.MDRD.BL.60LESS = EGFR.MDRD.BL < 60,
           EGFR.MDRD.IVRS.SCRN.30LESS = EGFR.MDRD.IVRS.SCRN < 30,
           EGFR.MDRD.IVRS.SCRN.60LESS = EGFR.MDRD.IVRS.SCRN < 60,
           EGFR.MDRD.CALC.SCRN.30LESS = EGFR.MDRD.CALC.SCRN < 30,
           EGFR.MDRD.CALC.SCRN.60LESS = EGFR.MDRD.CALC.SCRN < 60,
           HBA1C.8.3LESS = HBA1CBL <= 8.3)

W <- W[, order(colnames(W))]
W <- dplyr::select(W, c(USUBJID, ARM, contains("(LESS)|(PLUS)")), everything())

# Match AGE & AGE.GRP
W[is.na(W$AGE), "AGE.GRP"] <- NA


# combine with old W ------------------------------------------------------
W_new <- W

source(file = here("R/functions/LEADER_W_clean.R"))
W_new <- dplyr::select(W, -c("USUBJID")) %>%
    full_join(W_new, .,
              by = c("ARM", "AHYPERFL", "ANTDBFL", "CHDFL", "CHFFL", "CREATBL",
                     "DIABPBL", "GERDBLFL", "H2BLFL", "HBA1CBL", "HYPFL", "IHDFL",
                     "INSNVFL", "KIDFL", "LVSDFL", "MICFL", "MIFL", "NEPSCRFL",
                     "PADFL", "PPIFL", "PULSEBL", "RACE", "RETSCRFL", "REVASFL",
                     "SEX", "SMOKER", "STENFL", "STROKEFL", "SYSBPBL",
                     "EGFR.CKD.BL" = "EGFREPB", "AGE", "BMIBL", "CHOL1BL",
                     "DIABDUR", "HDL1BL", "LDL1BL", "WSTCIRBL", "TRIG1BL",
                     "EGFR.CKD.SCRN" = "EGFREPI",
                     "RENFSEV.CKD.BL" = "RENFCKD",
                     "RENFSEV.MDRD.BL" = "RENFSEV",
                     "RENFSEV.CKD.SCRN" = "RNFSCCKD",
                     "RENFSEV.MDRD.SCRN" = "RNFSCSEV")) %>%
    as_tibble() %>%
    dplyr::select(-c("CVRISK"))

# imputation --------------------------------------------------------------

set.seed(123456)
imputed <- mice::mice(dplyr::select_at(W_new, vars(-c("USUBJID", "ARM",
                                                      "AGE.GRP", "INSNVFL"),
                                                   -contains("PLUS"), -contains("LESS"))),
                      m = 1, maxit = 20, visitSequence = "monotone")

W_imputed <- full_join(dplyr::select_if(W_new, ~!anyNA(.)),
                       complete(imputed)) %>%
    mutate(AGE.60PLUS = AGE >= 60,
           AGE.GRP = factor(case_when(AGE < 65 ~ "<65",
                                      (AGE >= 65) & (AGE < 75) ~ "65-74",
                                      (AGE >= 75) & (AGE < 85) ~ "75-84",
                                      AGE >= 85 ~ ">=85"),
                            ordered = T,
                            levels = c("<65", "65-74", "75-84", ">=85"),
                            labels = c("<65", "65-74", "75-84", ">=85")),
           BMI.30PLUS = BMIBL > 30,
           DDUR.12PLUS = DIABDUR >= 12,
           EGFR.CKD.BL.30LESS = EGFR.CKD.BL < 30,
           EGFR.CKD.BL.60LESS = EGFR.CKD.BL < 60,
           EGFR.CKD.SCRN.30LESS = EGFR.CKD.SCRN < 30,
           EGFR.CKD.SCRN.60LESS = EGFR.CKD.SCRN < 60,
           EGFR.MDRD.BL.30LESS = EGFR.MDRD.BL < 30,
           EGFR.MDRD.BL.60LESS = EGFR.MDRD.BL < 60,
           EGFR.MDRD.IVRS.SCRN.30LESS = EGFR.MDRD.IVRS.SCRN < 30,
           EGFR.MDRD.IVRS.SCRN.60LESS = EGFR.MDRD.IVRS.SCRN < 60,
           EGFR.MDRD.CALC.SCRN.30LESS = EGFR.MDRD.CALC.SCRN < 30,
           EGFR.MDRD.CALC.SCRN.60LESS = EGFR.MDRD.CALC.SCRN < 60,
           HBA1C.8.3LESS = HBA1CBL <= 8.3)
if (!include_slices) {
    W_imputed <- W_imputed %>%
        dplyr::select(-c(AGE.60PLUS, AGE.GRP, BMI.30PLUS, DDUR.12PLUS,
                         EGFR.CKD.BL.30LESS, EGFR.CKD.BL.60LESS,
                         EGFR.CKD.SCRN.30LESS, EGFR.CKD.SCRN.60LESS,
                         EGFR.MDRD.BL.30LESS, EGFR.MDRD.BL.60LESS,
                         EGFR.MDRD.IVRS.SCRN.30LESS,
                         EGFR.MDRD.IVRS.SCRN.60LESS,
                         EGFR.MDRD.CALC.SCRN.30LESS,
                         EGFR.MDRD.CALC.SCRN.60LESS,
                         HBA1C.8.3LESS))
}


W_imputed <- W_new[, colnames(W_new) %in% colnames(W_imputed)] %>%
    dplyr::select_if(anyNA) %>%
    mutate_all(~case_when(is.na(.) ~ T,
                          T ~ F)) %>%
    rename_all(~paste0(., ".imp")) %>%
    cbind(W_imputed, .) %>% as_tibble()

W_imputed <- W_imputed[, order(colnames(W_imputed))]
W_imputed <- dplyr::select(W_imputed, c("USUBJID", "ARM"), everything())

## 25(or 15) imputation flags, 65 covariates, 15 binning/cutoffs, + USUBJID and ARM - 107 columns total ----


outcomes <- haven::read_sas(paste0(data_path, "adtte.sas7bdat")) %>%
    dplyr::select("USUBJID", 'PARAMCD', "AVAL", "PARCAT1") %>%
    dplyr::filter(PARAMCD %in% c("MACEEVTM", "MCECVDTM", "MACEMITM",
                                 "MCENFSTM", "NONCVTM"# , "PRMACETM",
    )) %>%
    rename("EVENT" = "PARCAT1") %>%
    mutate(PARAMCD = case_when(PARAMCD == "MACEEVTM" ~ "MACE",
                               PARAMCD == "MACEMITM" ~ "NFMI",
                               PARAMCD == "MCECVDTM" ~ "CVDEATH",
                               PARAMCD == "MCENFSTM" ~ "NFSTROKE",
                               PARAMCD == "NONCVTM" ~ 'NONCVDEATH' #,
                               # PARAMCD == "PRMACETM" ~ 'PRIMACE'
    )) %>%
    pivot_wider(id_cols = "USUBJID", names_from = "PARAMCD", values_from = c("AVAL", "EVENT")) %>%
    mutate_at(vars(starts_with("EVENT")),
              ~as.numeric(. == "TIME TO EVENT"))

