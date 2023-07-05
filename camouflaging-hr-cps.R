--------------------------------------------------------------------------------
# Hierarchical Regression Analysis to Identify Demographic Predictors of CAT-Q
--------------------------------------------------------------------------------
  
suppressPackageStartupMessages({
  library(tidyverse)
  library(sjmisc)
  library(car)
  library(sensemakr)
  library(emmeans)
})

# Read data
df <- read.csv("survey-final.csv",TRUE,",")

# Cleaning----------------------------------------------------------------------
df <- na.omit(df)
df <- filter(df, age < 400)
view(df)

## Recode demographic variables 
df$gender_coded <- rec(
  df$gender,
  rec = "Man (also referred to as cisgender man) = 1; Woman (also referred to
  as cisgender woman) = 2; else = 3" 
) 
df$diagnosis_recoded <- rec(
  df$diagnosis_coded,
  rec = "none = 1; chronic medical conditions = 2; substance use disorders = 3;  
  neurodivergent = 4; mental health diagnoses = 5; multidiagnosis = 6;
  else = 7" 
) 
df$ethnicity_coded <- rec(
  df$ethnicity, 
  rec = "Caucasian = 1; Black or African American = 2; American or Alaskan 
  Native = 3; Asian = 4; Native Hawaiian or Pacific Islander = 5; else = 6"
)
df$sexual_orient_coded <- rec(
  df$sexual_orient, 
  rec = "Straight or heterosexual = 1; Lesbian, Gay = 2; Bisexual = 3; else = 4"
)

df$gender_coded = as.factor(df$gender_coded)
df$diagnosis_recoded = as.factor(df$diagnosis_recoded)
df$ethnicity_coded = as.factor(df$ethnicity_coded)
df$sexual_orient_coded = as.factor(df$sexual_orient_coded)

# Reversing measures------------------------------------------------------------
# Reversing AEFI (higher = greater EF)
df$AEFI_1r = recode(df$AEFI_1, "1=3; 2=2; 3=1")
df$AEFI_2r = recode(df$AEFI_2, "1=3; 2=2; 3=1")
df$AEFI_3r = recode(df$AEFI_3, "1=3; 2=2; 3=1")
df$AEFI_5r = recode(df$AEFI_5, "1=3; 2=2; 3=1")
df$AEFI_7r = recode(df$AEFI_7, "1=3; 2=2; 3=1")
df$AEFI_8r = recode(df$AEFI_8, "1=3; 2=2; 3=1")
df$AEFI_9r = recode(df$AEFI_9, "1=3; 2=2; 3=1")
df$AEFI_10r = recode(df$AEFI_10, "1=3; 2=2; 3=1")

# Reversing CFS (higher = greater CFS)
df$CFS_2r = recode(df$CFS_2, "1=6; 2=5; 3=4; 4=3; 5=2; 6=1")
df$CFS_3r = recode(df$CFS_3, "1=6; 2=5; 3=4; 4=3; 5=2; 6=1")
df$CFS_5r = recode(df$CFS_5, "1=6; 2=5; 3=4; 4=3; 5=2; 6=1")
df$CFS_10r = recode(df$CFS_10, "1=6; 2=5; 3=4; 4=3; 5=2; 6=1")

# Reversing IRI (higher = greater PT)
df$IRI_1r = recode(df$IRI_1, "0=4; 1=3; 2=2; 3=1; 4=0")
df$IRI_4r = recode(df$IRI_4, "0=4; 1=3; 2=2; 3=1; 4=0")

# Reversing BIS (higher = worse inhib)
df$BIS_1r = recode(df$BIS_1, "1=4; 2=3; 3=2; 4=1")
df$BIS_4r = recode(df$BIS_4, "1=4; 2=3; 3=2; 4=1")
df$BIS_5r = recode(df$BIS_5, "1=4; 2=3; 3=2; 4=1")
df$BIS_6r = recode(df$BIS_6, "1=4; 2=3; 3=2; 4=1")

# Reversing INCOM (higher = greater soc comparison)
df$INCOM_3r = recode(df$INCOM_3, "1=5; 2=4; 3=3; 4=2; 5=1")

# Reversing INQ (higher = greater bel needs)
df$INQ_1r = recode(df$INQ_1, "1=7; 2=6; 3=5; 4=4; 5=3; 6=2; 7=1")
df$INQ_2r = recode(df$INQ_2, "1=7; 2=6; 3=5; 4=4; 5=3; 6=2; 7=1")
df$INQ_4r = recode(df$INQ_4, "1=7; 2=6; 3=5; 4=4; 5=3; 6=2; 7=1")
df$INQ_7r = recode(df$INQ_7, "1=7; 2=6; 3=5; 4=4; 5=3; 6=2; 7=1")
df$INQ_8r = recode(df$INQ_8, "1=7; 2=6; 3=5; 4=4; 5=3; 6=2; 7=1")
df$INQ_9r = recode(df$INQ_9, "1=7; 2=6; 3=5; 4=4; 5=3; 6=2; 7=1")

# Reversing ISMI (higher = greater IS)
df$ISMI_2r = recode(df$ISMI_2, "1=4; 2=3; 3=2; 4=1")
df$ISMI_9r = recode(df$ISMI_9, "1=4; 2=3; 3=2; 4=1")

# Reversing CATQ (higher = greater camouflaging)
df$CATQ_3r = recode(df$CATQ_3, "1=7; 2=6; 3=5; 4=4; 5=3; 6=2; 7=1")
df$CATQ_12r = recode(df$CATQ_12, "1=7; 2=6; 3=5; 4=4; 5=3; 6=2; 7=1")
df$CATQ_19r = recode(df$CATQ_19, "1=7; 2=6; 3=5; 4=4; 5=3; 6=2; 7=1")
df$CATQ_22r = recode(df$CATQ_22, "1=7; 2=6; 3=5; 4=4; 5=3; 6=2; 7=1")
df$CATQ_24r = recode(df$CATQ_24, "1=7; 2=6; 3=5; 4=4; 5=3; 6=2; 7=1")

# Reversing SMS (higher = greater SM)
df$SMS_9r = recode(df$SMS_9, "0=5; 1=4; 2=3; 3=2; 4=1; 5=0")
df$SMS_12r = recode(df$SMS_12, "0=5; 1=4; 2=3; 3=2; 4=1; 5=0")

# Reversing SRF* (higher = greater fatigue)
df$SRF_1r = recode(df$SRF_1, "1=5; 2=4; 3=3; 4=2; 5=1")
df$SRF_2r = recode(df$SRF_2, "1=5; 2=4; 3=3; 4=2; 5=1")
df$SRF_5r = recode(df$SRF_5, "1=5; 2=4; 3=3; 4=2; 5=1")
df$SRF_9r = recode(df$SRF_9, "1=5; 2=4; 3=3; 4=2; 5=1")
df$SRF_10r = recode(df$SRF_10, "1=5; 2=4; 3=3; 4=2; 5=1")
df$SRF_12r = recode(df$SRF_12, "1=5; 2=4; 3=3; 4=2; 5=1")
df$SRF_15r = recode(df$SRF_15, "1=5; 2=4; 3=3; 4=2; 5=1")
df$SRF_18r = recode(df$SRF_18, "1=5; 2=4; 3=3; 4=2; 5=1")

# Reversing KGAI* (higher = greater authenticity), KGAI2 check***
df$KGAI_1r = recode(df$KGAI_1, "1=5; 2=4; 3=3; 4=2; 5=1")
df$KGAI_2r = recode(df$KGAI_2, "1=5; 2=4; 3=3; 4=2; 5=1")
df$KGAI_3r = recode(df$KGAI_3, "1=5; 2=4; 3=3; 4=2; 5=1")
df$KGAI_4r = recode(df$KGAI_4, "1=5; 2=4; 3=3; 4=2; 5=1")
df$KGAI_5r = recode(df$KGAI_5, "1=5; 2=4; 3=3; 4=2; 5=1")
df$KGAI_6r = recode(df$KGAI_6, "1=5; 2=4; 3=3; 4=2; 5=1")
df$KGAI_7r = recode(df$KGAI_7, "1=5; 2=4; 3=3; 4=2; 5=1")
df$KGAI_8r = recode(df$KGAI_8, "1=5; 2=4; 3=3; 4=2; 5=1")
df$KGAI_9r = recode(df$KGAI_9, "1=5; 2=4; 3=3; 4=2; 5=1")
df$KGAI_10r = recode(df$KGAI_10, "1=5; 2=4; 3=3; 4=2; 5=1")
df$KGAI_11r = recode(df$KGAI_11, "1=5; 2=4; 3=3; 4=2; 5=1")

# Reversing SCC (higher = greater clarity)
df$SCC_1r = recode(df$SCC_1, "1=5; 2=4; 3=3; 4=2; 5=1")
df$SCC_2r = recode(df$SCC_2, "1=5; 2=4; 3=3; 4=2; 5=1")
df$SCC_3r = recode(df$SCC_3, "1=5; 2=4; 3=3; 4=2; 5=1")
df$SCC_4r = recode(df$SCC_4, "1=5; 2=4; 3=3; 4=2; 5=1")
df$SCC_5r = recode(df$SCC_5, "1=5; 2=4; 3=3; 4=2; 5=1")
df$SCC_7r = recode(df$SCC_7, "1=5; 2=4; 3=3; 4=2; 5=1")
df$SCC_8r = recode(df$SCC_8, "1=5; 2=4; 3=3; 4=2; 5=1")
df$SCC_9r = recode(df$SCC_9, "1=5; 2=4; 3=3; 4=2; 5=1")
df$SCC_10r = recode(df$SCC_10, "1=5; 2=4; 3=3; 4=2; 5=1")
df$SCC_12r = recode(df$SCC_12, "1=5; 2=4; 3=3; 4=2; 5=1")

# Reversing SATQ* (higher = greater autistic traits)
df$SATQ_1r = recode(df$SATQ_1, "0=3; 1=2; 2=1; 3=0")
df$SATQ_3r = recode(df$SATQ_3, "0=3; 1=2; 2=1; 3=0")
df$SATQ_4r = recode(df$SATQ_4, "0=3; 1=2; 2=1; 3=0")
df$SATQ_5r = recode(df$SATQ_5, "0=3; 1=2; 2=1; 3=0")
df$SATQ_7r = recode(df$SATQ_7, "0=3; 1=2; 2=1; 3=0")
df$SATQ_9r = recode(df$SATQ_9, "0=3; 1=2; 2=1; 3=0")
df$SATQ_11r = recode(df$SATQ_11, "0=3; 1=2; 2=1; 3=0")
df$SATQ_12r = recode(df$SATQ_12, "0=3; 1=2; 2=1; 3=0")
df$SATQ_13r = recode(df$SATQ_13, "0=3; 1=2; 2=1; 3=0")
df$SATQ_14r = recode(df$SATQ_14, "0=3; 1=2; 2=1; 3=0")
df$SATQ_15r = recode(df$SATQ_15, "0=3; 1=2; 2=1; 3=0")
df$SATQ_17r = recode(df$SATQ_17, "0=3; 1=2; 2=1; 3=0")
df$SATQ_19r = recode(df$SATQ_19, "0=3; 1=2; 2=1; 3=0")
df$SATQ_21r = recode(df$SATQ_21, "0=3; 1=2; 2=1; 3=0")
df$SATQ_23r = recode(df$SATQ_23, "0=3; 1=2; 2=1; 3=0")

# Create constructs-------------------------------------------------------------
# AEFI factor
df$AEFI = df$AEFI_1r + df$AEFI_2r + df$AEFI_3r + df$AEFI_4 + df$AEFI_5r +
  df$AEFI_6 + df$AEFI_7r + df$AEFI_8r + df$AEFI_9r + df$AEFI_10r

# ASRS factor
df$ASRS = df$ASRS_1 + df$ASRS_2 + df$ASRS_3 + df$ASRS_4 + df$ASRS_5 + df$ASRS_6

# COMP factor
df$COMP = df$COMP_1 + df$COMP_2 + df$COMP_3 + df$COMP_4 + df$COMP_5 + 
  df$COMP_6 + df$COMP_7 + df$COMP_8 + df$COMP_9 + df$COMP_10 + df$COMP_11 +
  df$COMP_12 + df$COMP_13 + df$COMP_14 + df$COMP_15 + df$COMP_16 + df$COMP_17 +
  df$COMP_18 + df$COMP_19 + df$COMP_20 + df$COMP_21 + df$COMP_22 + df$COMP_23 +
  df$COMP_24 + df$COMP_25

# CFS factor
df$CFS = df$CFS_1 + df$CFS_2r + df$CFS_3r + df$CFS_4 + df$CFS_5r + df$CFS_6 +
  df$CFS_7 + df$CFS_8 + df$CFS_9 + df$CFS_10r + df$CFS_11 + df$CFS_12

# IRI factor
df$IRI = df$IRI_1r + df$IRI_2 + df$IRI_3 + df$IRI_4r + df$IRI_5 + df$IRI_6 +
  df$IRI_7

# BIS factor
df$BIS = df$BIS_1r + df$BIS_2 + df$BIS_3 + df$BIS_4r + df$BIS_5r + df$BIS_6r +
  df$BIS_7 + df$BIS_8

# ERQ factor
df$ERQ = df$ERQ_1 + df$ERQ_2 + df$ERQ_3 + df$ERQ_4 + df$ERQ_5 + df$ERQ_6 +
  df$ERQ_7 + df$ERQ_8 + df$ERQ_9 + df$ERQ_10

# SelfCS factor
df$SelfCS = df$SelfCS_1 + df$SelfCS_2 + df$SelfCS_3 + df$SelfCS_4 +
  df$SelfCS_5 + df$SelfCS_6 + df$SelfCS_7

# INCOM factor
df$INCOM = df$INCOM_1 + df$INCOM_2 + df$INCOM_3r + df$INCOM_4 + df$INCOM_5 +
  df$INCOM_6

# INQ factor
df$INQ = df$INQ_1r + df$INQ_2r + df$INQ_3 + df$INQ_4r + df$INQ_5 + df$INQ_6 +
  df$INQ_7r + df$INQ_8r + df$INQ_9r

# ISMI factor
df$ISMI = df$ISMI_1 + df$ISMI_2r + df$ISMI_3 + df$ISMI_4 + df$ISMI_5 +
  df$ISMI_6 + df$ISMI_7 + df$ISMI_8 + df$ISMI_9r + df$ISMI_10

# RSE factor
df$RSE = df$RSE_1 + df$RSE_2 + df$RSE_3 + df$RSE_4 + df$RSE_5 + df$RSE_6 +
  df$RSE_7 + df$RSE_8

# PSCS factor
df$PSCS = df$PSCS_1 + df$PSCS_2 + df$PSCS_3 + df$PSCS_4

# CATQ factor
df$CATQ = df$CATQ_1 + df$CATQ_2 + df$CATQ_3r + df$CATQ_4 + df$CATQ_5 +
  df$CATQ_6 + df$CATQ_7 + df$CATQ_8 + df$CATQ_9 + df$CATQ_10 + df$CATQ_11 +
  df$CATQ_12r + df$CATQ_13 + df$CATQ_14 + df$CATQ_15 + df$CATQ_16 + df$CATQ_17 +
  df$CATQ_18 + df$CATQ_19r + df$CATQ_20 + df$CATQ_21 + df$CATQ_22r + 
  df$CATQ_23 + df$CATQ_24r + df$CATQ_25

# SMS factor
df$SMS = df$SMS_1 + df$SMS_2 + df$SMS_3 + df$SMS_4 + df$SMS_5 + df$SMS_6 +
  df$SMS_7 + df$SMS_8 + df$SMS_9r + df$SMS_10 + df$SMS_11 + df$SMS_12r +
  df$SMS_13

# SPT factor
df$SPT = df$SPT_1 + df$SPT_2 + df$SPT_3 + df$SPT_4 + df$SPT_5 + df$SPT_6 +
  df$SPT_7 + df$SPT_8 + df$SPT_9 + df$SPT_10 + df$SPT_11 + df$SPT_12 +
  df$SPT_13 + df$SPT_14 + df$SPT_15 + df$SPT_16 + df$SPT_17 + df$SPT_18 +
  df$SPT_19 + df$SPT_20 + df$SPT_21 + df$SPT_22 + df$SPT_23 + df$SPT_24 +
  df$SPT_25 + df$SPT_26 + df$SPT_27 + df$SPT_28 + df$SPT_29 + df$SPT_30 +
  df$SPT_31 + df$SPT_32 + df$SPT_33 + df$SPT_34 + df$SPT_35 + df$SPT_36 +
  df$SPT_37 + df$SPT_38 + df$SPT_39 + df$SPT_40 + df$SPT_41 + df$SPT_42 +
  df$SPT_43 + df$SPT_44 + df$SPT_45 + df$SPT_46 + df$SPT_47 + df$SPT_48 +
  df$SPT_49 + df$SPT_50 + df$SPT_51 + df$SPT_52 + df$SPT_53 + df$SPT_54 +
  df$SPT_55 + df$SPT_56 + df$SPT_57 + df$SPT_58 + df$SPT_59 + df$SPT_60 +
  df$SPT_61 + df$SPT_62 + df$SPT_63

# SCS factor
df$SCS = df$SCS_1 + df$SCS_2 + df$SCS_3 + df$SCS_4 + df$SCS_5 + df$SCS_6 +
  df$SCS_7 + df$SCS_8 + df$SCS_9 + df$SCS_10

# GAD factor
df$GAD = df$GAD_1 + df$GAD_2 + df$GAD_3 + df$GAD_4 + df$GAD_5 + df$GAD_6 +
  df$GAD_7

# PHQ factor
df$PHQ = df$PHQ_1 + df$PHQ_2 + df$PHQ_3 + df$PHQ_4 + df$PHQ_5 + df$PHQ_6 +
  df$PHQ_7 + df$PHQ_8 + df$PHQ_9

# LSAS anxiety factor
df$LSAS_anx = df$LSAS_1_1 + df$LSAS_2_1 + df$LSAS_3_1 + df$LSAS_4_1 +
  df$LSAS_5_1 + df$LSAS_6_1 + df$LSAS_7_1 + df$LSAS_8_1 + df$LSAS_9_1 +
  df$LSAS_10_1 + df$LSAS_11_1 + df$LSAS_12_1 + df$LSAS_13_1 + df$LSAS_14_1 +
  df$LSAS_15_1 + df$LSAS_16_1 + df$LSAS_17_1 + df$LSAS_18_1 + df$LSAS_19_1 +
  df$LSAS_20_1 + df$LSAS_21_1 + df$LSAS_22_1 + df$LSAS_23_1 + df$LSAS_24_1

# LSAS avoidance factor
df$LSAS_avd = df$LSAS_1_2 + df$LSAS_2_2 + df$LSAS_3_2 + df$LSAS_4_2 +
  df$LSAS_5_2 + df$LSAS_6_2 + df$LSAS_7_2 + df$LSAS_8_2 + df$LSAS_9_2 +
  df$LSAS_10_2 + df$LSAS_11_2 + df$LSAS_12_2 + df$LSAS_13_2 + df$LSAS_14_2 +
  df$LSAS_15_2 + df$LSAS_16_2 + df$LSAS_17_2 + df$LSAS_18_2 + df$LSAS_19_2 +
  df$LSAS_20_2 + df$LSAS_21_2 + df$LSAS_22_2 + df$LSAS_23_2 + df$LSAS_24_2

# SWEMWS factor
df$SWEMWS = df$SWEMWS_1 + df$SWEMWS_2 + df$SWEMWS_3 + df$SWEMWS_4 +
  df$SWEMWS_5 + df$SWEMWS_6 + df$SWEMWS_7

# SRF factor
df$SRF = df$SRF_1r + df$SRF_2r + df$SRF_3 + df$SRF_4 + df$SRF_5r + df$SRF_6 +
  df$SRF_7 + df$SRF_8 + df$SRF_9r + df$SRF_10r + df$SRF_11 + df$SRF_12r +
  df$SRF_13 + df$SRF_14 + df$SRF_15r + df$SRF_16 + df$SRF_17 + df$SRF_18r

# KGAI factor
df$KGAI = df$KGAI_1r + df$KGAI_2r + df$KGAI_3r + df$KGAI_4r + df$KGAI_5r +
  df$KGAI_6r + df$KGAI_7r + df$KGAI_8r + df$KGAI_9r + df$KGAI_10r + 
  df$KGAI_11r + df$KGAI_12 + df$KGAI_13 + df$KGAI_14 + df$KGAI_15 + df$KGAI_16 +   
  df$KGAI_17 + df$KGAI_18 + df$KGAI_19 + df$KGAI_20

# SCC factor
df$SCC = df$SCC_1r + df$SCC_2r + df$SCC_3r + df$SCC_4r + df$SCC_5r + 
  df$SCC_6 + df$SCC_7r + df$SCC_8r + df$SCC_9r + df$SCC_10r + df$SCC_11 + 
  df$SCC_12r

# SATQ factor
df$SATQ = df$SATQ_1r + df$SATQ_2 + df$SATQ_3r + df$SATQ_4r + df$SATQ_5r +
  df$SATQ_6 + df$SATQ_7r + df$SATQ_8 + df$SATQ_9r + df$SATQ_10 + df$SATQ_11r +
  df$SATQ_12r + df$SATQ_13r + df$SATQ_14r + df$SATQ_15r + df$SATQ_16 + 
  df$SATQ_17r + df$SATQ_18 + df$SATQ_19r + df$SATQ_20 + df$SATQ_21r + 
  df$SATQ_22 + df$SATQ_23r + df$SATQ_24

# Standardizing--------------------------------------------------
df$CATQ_sd <- scale(df$CATQ)
df$education_years_1_sd <- scale(df$education_years_1)
df$age_sd <- scale(df$age)
df$ASRS_sd <- scale(df$ASRS)
df$SATQ_sd <- scale(df$SATQ)

# HLR: Predictors of CATQ--------------------------------------------------------
## Step/Block 1 (demographics covariates): 
step1 <- lm(
  CATQ_sd ~ ethnicity_coded + gender_coded + sexual_orient_coded + 
    education_years_1_sd + age_sd, data = df
)
summary(step1)

## Step/Block 2 (added diagnosis info)
step2 <- lm(
  CATQ_sd ~ ethnicity_coded + gender_coded + sexual_orient_coded + 
    education_years_1_sd + age_sd + diagnosis_recoded + ASRS_sd + SATQ_sd, 
  data = df
)
summary(step2)
confint(step2, level = 0.95) # find confidence intervals of regression coefs

## Correct for multiple comparisons with BH and extract partial coef
pvals <- summary(step2)$coefficients[,4]
q <- p.adjust(pvals, method = "BH")
q 

partial_r2(step2, "gender_coded2")
partial_r2(step2, "age_sd")
partial_r2(step2, "SATQ_sd")
partial_r2(step2, "ASRS_sd")

partial_f2(step2, "gender_coded2")
partial_f2(step2, "age_sd")
partial_f2(step2, "SATQ_sd")
partial_f2(step2, "ASRS_sd")

## Model comparison
anova(step1, step2) # model 2 better

## Change in R-squared
summary(step2)$r.squared - summary(step1)$r.squared
