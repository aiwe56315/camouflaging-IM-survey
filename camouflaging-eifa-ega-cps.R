--------------------------------------------------------------------------------
# Exploratory Item Factor Analysis and Exploratory Graph Analysis of CAT-Q's 
# Dimensional Structure 
--------------------------------------------------------------------------------
  
library(tidyverse)
library(stats)
library(psych)
library(REdaS)
library(GPArotation)
library(lavaan)
library(haven)
library(semPlot)
library(psy)
library(mirt)
library(EFAtools)
library(ltm)
library(car)
library(writexl)
library(qgraph)
library(bootnet)
library(dplyr)
library(network)
library(igraph)
library(ggraph)
library(tidygraph)
library(visNetwork)
library(networkD3)
library(bootnet)
library(networktools)
library(NetworkComparisonTest)
library(EGAnet)
library(sna)
library(sjmisc)
library(ggpubr)

# Read data
df <- read.csv("survey-final.csv",TRUE,",")

# Cleaning and preparations-----------------------------------------------------
df <- na.omit(df)
df <- filter(df, age < 400)

## Reversing CATQ 
df$CATQ_3r = recode(df$CATQ_3, "1=7; 2=6; 3=5; 4=4; 5=3; 6=2; 7=1")
df$CATQ_12r = recode(df$CATQ_12, "1=7; 2=6; 3=5; 4=4; 5=3; 6=2; 7=1")
df$CATQ_19r = recode(df$CATQ_19, "1=7; 2=6; 3=5; 4=4; 5=3; 6=2; 7=1")
df$CATQ_22r = recode(df$CATQ_22, "1=7; 2=6; 3=5; 4=4; 5=3; 6=2; 7=1")
df$CATQ_24r = recode(df$CATQ_24, "1=7; 2=6; 3=5; 4=4; 5=3; 6=2; 7=1")

## Build data frame 
dfCATQ <- data.frame(
  df$CATQ_1, df$CATQ_2, df$CATQ_3r, df$CATQ_4, df$CATQ_5, df$CATQ_6, df$CATQ_7,
  df$CATQ_8, df$CATQ_9, df$CATQ_10, df$CATQ_11, df$CATQ_12r, df$CATQ_13, 
  df$CATQ_14, df$CATQ_15, df$CATQ_16, df$CATQ_17, df$CATQ_18, df$CATQ_19r, 
  df$CATQ_20, df$CATQ_21, df$CATQ_22r, df$CATQ_23, df$CATQ_24r, df$CATQ_25
)

## Data partition 
set.seed(1234)
ind <- sort(sample(nrow(dfCATQ), nrow(dfCATQ)*.5))
dfCATQ_eifa <- dfCATQ[ind,]
dfCATQ_ega <- dfCATQ[-ind,]

## Check balance of data partition----------------------------------------------
df$CATQ = df$CATQ_1 + df$CATQ_2 + df$CATQ_3r + df$CATQ_4 + df$CATQ_5 +
  df$CATQ_6 + df$CATQ_7 + df$CATQ_8 + df$CATQ_9 + df$CATQ_10 + df$CATQ_11 +
  df$CATQ_12r + df$CATQ_13 + df$CATQ_14 + df$CATQ_15 + df$CATQ_16 + df$CATQ_17 +
  df$CATQ_18 + df$CATQ_19r + df$CATQ_20 + df$CATQ_21 + df$CATQ_22r + 
  df$CATQ_23 + df$CATQ_24r + df$CATQ_25

set.seed(1234) # Sane sample structure and distribution 
ind <- sort(sample(nrow(df), nrow(df)*.5))
df1 <- df[ind,]
df2 <- df[-ind,]

### CAT-Q scores
mean((df1$CATQ)/25)
sd((df1$CATQ)/25)

mean((df2$CATQ)/25)
sd((df2$CATQ)/25)

CATQ_res <- t.test(df1$CATQ, df2$CATQ)
CATQ_res

### Age difference
mean(df1$age)
sd(df1$age)

mean(df2$age)
sd(df2$age)

age_res <- t.test(df1$age, df2$age)
age_res

### Gender/sex proportion difference
df1$gender_coded <- rec(
  df1$gender,
  rec = "Man (also referred to as cisgender man) = 1; Woman (also referred to
  as cisgender woman) = 2; else = 3"
) 
df1$gender_coded <- as.factor(df1$gender_coded)
df1$sex <- as.factor(df1$sex)
summary(df1$gender_coded)
summary(df1$sex)

df2$gender_coded <- rec(
  df2$gender,
  rec = "Man (also referred to as cisgender man) = 1; Woman (also referred to
  as cisgender woman) = 2; else = 3"
) 
df2$gender_coded <- as.factor(df2$gender_coded)
df2$sex <- as.factor(df2$sex)
summary(df2$gender_coded)
summary(df2$sex)

### Ethnicity proportion difference 
df1$ethnicity = as.factor(df1$ethnicity)
summary(df1$ethnicity)

df2$ethnicity = as.factor(df2$ethnicity)
summary(df2$ethnicity)

# Run Chi-Square equality of proportions test (keep format, replace with summary
#stats from proportions)
res <- prop.test(x = c(10, 10), n = c(486, 486))
res

# Scree plot and K1-------------------------------------------------------------
## Set up
pc <- prcomp(dfCATQ_eifa, scale = TRUE)
pc

eigenvalues <- pc$sdev * pc$sdev
eigenvectors <- pc$rotation
summary(pc)

## Visualizations 
screeplot(pc, type = "l", main = "Screeplot of CATQ PCA")
abline(1, 0, col = "red", lty = 2)
scree.plot(dfCATQ, title = "Screeplot of Eigenvalues", simu = "F") # Kaiser rule
biplot(pc, scale = 0) # messy, not too helpful 

# Exploratory Item Factor Analysis (EIFA)---------------------------------------
itemstats(dfCATQ_eifa)

## Unidimensional 
mirt_mod1 <- mirt(
  data = dfCATQ_eifa,
  1,
  rotate = "oblimin",
  method = "MHRM",
  verbose = TRUE
)
mirt_mod1
M2(mirt_mod1, CI = .95)
out <- summary(mirt_mod1, rotate = "oblimin", suppress = .3) 
pander::pander(out)
coef(mirt_mod1)
fs1 <- fscores(mirt_mod1, QMC = TRUE)
fs1

## 2-Factor 
mirt_mod2 <- mirt(
  data = dfCATQ_eifa,
  2,
  rotate = "oblimin",
  method = "MHRM",
  verbose = TRUE
)
mirt_mod2
M2(mirt_mod2, CI = .95)
out <- summary(mirt_mod2, rotate = "oblimin", suppress = .4) 
pander::pander(out)
coef(mirt_mod2)
fs1 <- fscores(mirt_mod2, QMC = TRUE)
fs1

## 3-Factor
mirt_mod3 <- mirt(
  data = dfCATQ_eifa,
  3,
  rotate = "oblimin",
  method = "MHRM",
  verbose = TRUE
)
mirt_mod3
M2(mirt_mod3, CI = .95, QMC = TRUE)
out <- summary(mirt_mod3, rotate = "oblimin") 
pander::pander(out)
coef(mirt_mod3)
fs3 <- fscores(mirt_mod3, QMC = TRUE)
fs3

## 4-Factor
mirt_mod4 <- mirt(
  data = dfCATQ_eifa,
  4,
  rotate = "oblimin",
  method = "MHRM",
  verbose = TRUE
)
mirt_mod4   
M2(mirt_mod4, CI = .95, QMC = TRUE)
out <- summary(mirt_mod4, rotate = "oblimin") 
pander::pander(out)
coef(mirt_mod4)
fs4 <- fscores(mirt_mod4, QMC = TRUE)
fs4

## Model comparisons 
anova(mirt_mod1, mirt_mod2, mirt_mod3, mirt_mod4)

# Internal consistency of 3-factor solution-------------------------------------
## Cronbach's alpha for composite CATQ
cronbach.alpha(dfCATQ, CI = TRUE)

## Cronbach's alpha for Compensation subscale 
dfCATQ_comp <- data.frame(
  df$CATQ_1, df$CATQ_4, df$CATQ_5, df$CATQ_8, df$CATQ_11, df$CATQ_14, 
  df$CATQ_17, df$CATQ_20, df$CATQ_23 
)
cronbach.alpha(dfCATQ_comp, CI = TRUE)

## Cronbach's alpha for Masking subscale
dfCATQ_mask <- data.frame(
  df$CATQ_2, df$CATQ_6, df$CATQ_9, df$CATQ_15, df$CATQ_18, df$CATQ_21
)
cronbach.alpha(dfCATQ_mask, CI = TRUE)

## Cronbach's alpha for Assimilation subscale
dfCATQ_asm <- data.frame(
  df$CATQ_3r, df$CATQ_7, df$CATQ_10, df$CATQ_13, df$CATQ_16, df$CATQ_19r, 
  df$CATQ_22r, df$CATQ_25
)
cronbach.alpha(dfCATQ_asm, CI = TRUE)

## Parametric bootstrap likelihood-ratio (LR) test (whether more restrictive 
## model fits significantly worse than less restrictive model)
mirtCluster()
boot.LR(mirt_mod3, mirt_mod4, R = 100, verbose = TRUE)

mirt::itemplot(mirt_mod3, 3)

# Exploratory Graphical Analysis (EGA)------------------------------------------
## Rename column names for EGA prep
dfCATQ_ega <- dfCATQ_ega %>% 
  rename(CATQ1 = df.CATQ_1,
         CATQ2 = df.CATQ_2,
         CATQ3 = df.CATQ_3r,
         CATQ4 = df.CATQ_4,
         CATQ5 = df.CATQ_5,
         CATQ6 = df.CATQ_6,
         CATQ7 = df.CATQ_7,
         CATQ8 = df.CATQ_8,
         CATQ9 = df.CATQ_9,
         CATQ10 = df.CATQ_10,
         CATQ11 = df.CATQ_11,
         CATQ12 = df.CATQ_12r,
         CATQ13 = df.CATQ_13,
         CATQ14 = df.CATQ_14,
         CATQ15 = df.CATQ_15,
         CATQ16 = df.CATQ_16,
         CATQ17 = df.CATQ_17,
         CATQ18 = df.CATQ_18,
         CATQ19 = df.CATQ_19r,
         CATQ20 = df.CATQ_20,
         CATQ21 = df.CATQ_21,
         CATQ22 = df.CATQ_22r,
         CATQ23 = df.CATQ_23,
         CATQ24 = df.CATQ_24r,
         CATQ25 = df.CATQ_25
  )

## Create network using EBICglasso
names <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
           "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24",
           "25")

itemnames <- c("Copy body language","Monitor body lang/expressions (relaxed)", 
               "Rarely put on act","Developed script","Repeat phrases heard",
               "Adjust body lang/expressions (interested)",
               "Performing rather than being myself","Use learned behaviors",
               "Always think about impression","Need support to socialize",
               "Practice expressions/body lang","Eye contact",
               "Forced to interact","Improve social understanding",
               "Monitor body lang/expressions (interested)","Avoid interacting",
               "Research social skills","Always aware of impression",
               "Free to be myself", "Learn use body/faces from TV",
               "Adjust body lang/expressions (relaxed)", "Conversation flows",
               "Learn social skills from TV/films", 
               "Do not pay attention to face/body","Pretend to be normal")

grp <- list("Compensation" = c(1,4,5,8,11,14,17,20,23), "Masking" = 
              c(2,6,9,12,15,18,21,24), "Assimilation" = c(3,7,10,13,16,19,22,25))

## Create correlation network (basic)
corMat <- cor_auto(dfCATQ_ega)

View(corMat)
g1 <- autograph(corMat, labels = names, layout = "circle", vsize = 6, groups = grp, 
                color = c('#718a8a', '#d7a761', '#b75e50'), border.width = 1.5, 
                nodeNames = itemnames, legend.cex = .35, cut = 0, details = TRUE)
g1

## Regularized model 
GGM <- estimateNetwork(dfCATQ_ega, default = "EBICglasso") # lasso regularized 
GGM_plot <- plot(GGM, labels = names, layout = "spring", vsize = 6, cut = 0,
                 border.width = 1.5, border.color = "black", groups = grp, 
                 color = c('#718a8a', '#d7a761', '#b75e50'), 
                 nodeNames = itemnames, legend.cex = .4)
pdf("GGM.pdf", width = 9, height = 5)

c1 <- centralityPlot(GGM_plot, include = "ExpectedInfluence")
CentralityTable <- centralityTable(GGM)
write.csv(CentralityTable, "CentralityTable.csv")

# Estimate network stability 
set.seed(1234)
bootdata <- bootEGA(
  dfCATQ_ega,
  iter= 1000, 
  type = "resampling",
)

bootdata$summary.table # summarize stability metrics 
bootdata$frequency # frequency of factor solution replication
itemStability(bootdata) # item stability 

b1 <- bootnet(
  GGM, 
  boots = 1000, 
  nCores = 4, 
  statistics = c("strength", "expectedInfluence", "edge")
)

b2 <- bootnet(
  GGM, 
  boots = 1000, 
  nCores = 4, 
  type = "case",
  statistics = c("strength", "expectedInfluence", "edge")
)

save(b1, file = "b1.Rdata")
save(b2, file = "b2.Rdata")
load("b1.Rdata")
load("b2.Rdata")
corStability(b2)

# Alternative Procedure for EGA-------------------------------------------------
# Assign names to nodes
nodes <- c("CATQ1", "CATQ2", "CATQ3r", "CATQ4", "CATQ5", "CATQ6", "CATQ7", 
           "CATQ8", "CATQ9", "CATQ10", "CATQ11", "CATQ12r", "CATQ13", "CATQ14", 
           "CATQ15", "CATQ16", "CATQ17", "CATQ18", "CATQ19", "CATQ20", "CATQ21", 
           "CATQ22r", "CATQ23", "CATQ24r", "CATQ25")

# Estimate network using default methods
network <- estimateNetwork(dfCATQ_ega, default = "EBICglasso")

# Estimate network stability 
set.seed(1234)
b1 <- bootnet(
  network, 
  boots = 1000, 
  nCores = 4, 
  statistics = c("strength", "expectedInfluence", "edge")
)

b2 <- bootnet(
  network, 
  boots = 1000, 
  nCores = 4, 
  type = "case",
  statistics = c("strength", "expectedInfluence", "edge")
)
corStability(b2)

# Plot network
NetPlot <- plot(network, layout = "spring", vsize = 6, border.color = "black", 
                nodeNames = "nodes", color = "lightblue", legend = FALSE)

# Create igraph object, un-directed network, and convert all weights to absolute value
netI <- as.igraph(NetPlot, attributes = TRUE)
class(netI)

net.sym <- as.undirected(netI, mode= "collapse", edge.attr.comb=list(weight="sum", "ignore"))
E(net.sym)$weight <- abs(E(net.sym)$weight)
E(netI)$weight <- abs(E(netI)$weight)

# Community detection using Leiden 
ldc <- cluster_leiden(
  net.sym,
  objective_function = c("modularity"),
  weights = NULL,
  resolution_parameter = 1,
  beta = 0.01,
  initial_membership = NULL,
  n_iterations = 10,
  vertex_weights = NULL
)

# Community detection using FastGreedy (optimizes modularity score through 
# iterative sampling of random edges (kept if improve))
cfg <- cluster_fast_greedy(net.sym)
plot(cfg, net.sym, nodeNames = nodes)
plot_dendrogram(cfg)
length(cfg)

# Detected 3 communities (overlap with CATQ factor structure)
sizes(cfg)
membership(cfg)
modularity(cfg)



# Convergent validity test------------------------------------------------------
# CATQ factor
df$CATQ = df$CATQ_1 + df$CATQ_2 + df$CATQ_3r + df$CATQ_4 + df$CATQ_5 +
  df$CATQ_6 + df$CATQ_7 + df$CATQ_8 + df$CATQ_9 + df$CATQ_10 + df$CATQ_11 +
  df$CATQ_12r + df$CATQ_13 + df$CATQ_14 + df$CATQ_15 + df$CATQ_16 + df$CATQ_17 +
  df$CATQ_18 + df$CATQ_19r + df$CATQ_20 + df$CATQ_21 + df$CATQ_22r + 
  df$CATQ_23 + df$CATQ_24r + df$CATQ_25

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

cor.test(df$CATQ, df$SPT, method = "pearson")
