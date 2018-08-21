# !!! FIRST set your working directory to be wherever 'code_r' etc live !!!

# make sure we always use same pseudo-random numbers 
set.seed(27081975)
# run the IPM definition script
source(file.path(".", "code_r", "1_IPM_definition.R"))

#============================================================================
#
# (I) Basic LTRE analysis
#
#============================================================================

iPar <- mk_intpar(nm=100, l=7, u=20) # set up the iteration parameters
nSimYr <- 100000 # number of years to simulate (takes a few mins w/ 100k years)
envSet <- mk_env(mParSetStore, eStateSet, nYrs=nSimYr) # make a set of environments
system.time(
  simOut <- sim_IPM(envSet$mParSet, envSet$eStateSet, iPar) # simulate forwards
)

# save(envSet, simOut, file = "LTRE_data.rda")
# load(file = "LTRE_data.rda")

# "data" for the analysis
LTREdata <- within(envSet, {
  eNamesSuffix <- colnames(eStateSet)
  eNames <- paste("GrowAJ", eNamesSuffix, sep=".")
  GrowAJ.I0 <- mParSet[,"GrowAJ.I"] + apply(mParSet[,eNames] * eStateSet, 1, sum)
  GrowAJ.I1 <- mParSet[,"GrowAJ.Idiff"] + GrowAJ.I0
  eNames <- paste("GrowJA", eNamesSuffix, sep=".")
  GrowJA.I0 <- mParSet[,"GrowJA.I"] + apply(mParSet[,eNames] * eStateSet, 1, sum)
  GrowJA.I1 <- mParSet[,"GrowJA.Idiff"] + GrowJA.I0
  eNames <- paste("RecrSz", eNamesSuffix, sep=".")
  RecrSz.I <- mParSet[,"RecrSz.I"] + apply(mParSet[,eNames] * eStateSet, 1, sum)
  eNames <- paste("Surv",   eNamesSuffix, sep=".")
  Surv.I   <- mParSet[,"Surv.I"  ] + apply(mParSet[,eNames] * eStateSet, 1, sum)  
  eNames <- paste("Repr",   eNamesSuffix, sep=".")
  Repr.I   <- mParSet[,"Repr.I"  ] + apply(mParSet[,eNames] * eStateSet, 1, sum)  
  eNames <- paste("Recr",   eNamesSuffix, sep=".")
  Recr.I   <- mParSet[,"Recr.I"  ] + apply(mParSet[,eNames] * eStateSet, 1, sum)
  SprT <- eStateSet[, "SprT"]
  WinT <- eStateSet[, "WinT"]
  BrGd <- eStateSet[, "BrGd"]
  # clean up
  rm(eNamesSuffix, eNames)
})

# clean up
LTREdata$mParSet <- NULL
LTREdata$eStateSet <- NULL
LTREdata <- data.frame(as.data.frame(LTREdata), Lambda = simOut$lamt)

# add lagged variables
parNames <- c("GrowAJ.I0", "GrowJA.I0", "GrowAJ.I1", "GrowJA.I1", 
         "RecrSz.I", "Surv.I", "Repr.I", "Recr.I")
laggedData <- head(LTREdata, -1)[, parNames]
names(laggedData) <- paste(names(laggedData), "lag1", sep=".")
LTREdata <- cbind(tail(LTREdata, -1), laggedData)

# 
varLam <- var(log(LTREdata$Lambda))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Helper functions
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

calcLTREcontribs <- function(LTREterms) {
  allCoVar <- cov(LTREterms)
  isLT <- lower.tri(allCoVar, diag = FALSE)
  allCoVar[isLT] <- 2 * allCoVar[isLT]
  isUT <- upper.tri(allCoVar, diag = FALSE)
  allCoVar[isUT] <- NA
  covNames <- expand.grid(v1 = parNames, v2 = parNames, stringsAsFactors = FALSE)
  allCoVar <- data.frame(covNames, contrib = allCoVar[TRUE])
  tbl_df(allCoVar)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## A. Implement original Ellner / Rees approach using linear terms
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vector with names of the parameters
parNames <- c("GrowAJ.I0", "GrowJA.I0", "GrowAJ.I1", "GrowJA.I1", 
              "RecrSz.I", "Surv.I", "Repr.I", "Recr.I")

# create the linear predictor formula
form <- paste("log(Lambda) ~", paste(parNames, collapse=" + "))
# fit the linear model
linLTREmod <- lm(as.formula(form), data=LTREdata)
summary(linLTREmod) # 91.1% of the variance

# model predictions by term
linLTREterms <- predict(linLTREmod, type="terms")
colnames(linLTREterms) <- parNames
# store the associated covariances
allCoVar <- calcLTREcontribs(linLTREterms)
# extract the error variance of the model
residVar <- mean(residuals(linLTREmod)^2)

# sanity check--is the sum of model contribs equal to variance in lambda?
sum(allCoVar$contrib, na.rm = TRUE) + residVar
varLam

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## B. Repeat w/ linear terms *and* lag 1 parameters
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parNames <- c("GrowAJ.I0", "GrowJA.I0", "GrowAJ.I1", "GrowJA.I1", 
              "RecrSz.I", "Surv.I", "Repr.I", "Recr.I")
parNames <- c(parNames, paste(parNames, "lag1", sep="."))

# create the linear predictor formula
form <- paste("log(Lambda) ~", paste(parNames, collapse=" + "))
# fit the linear model
lagLTREmod <- lm(as.formula(form), data = LTREdata)
summary(lagLTREmod) # 97.5% of variance

# model predictions by term
lagLTREterms <- predict(lagLTREmod, type="terms")
colnames(lagLTREterms) <- parNames
# store the associated covariances in a named vector
allCoVar <- calcLTREcontribs(lagLTREterms)
# extract the error variance of the model
residVar <- mean(residuals(lagLTREmod)^2)

# sanity check--is the sum of model contribs equal to variance in lambda?
sum(allCoVar$contrib, na.rm = TRUE) + residVar
varLam

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## C. Repeat w/ smooth terms (i.e. gam) *and* lag 1 parameters
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 3. Repeat with lag 1 parameters, but now use smooth terms (i.e. a GAM)

parNames <- c("Surv.I", "Repr.I", "Recr.I", 
              "GrowAJ.I0", "GrowJA.I0", "GrowAJ.I1", "GrowJA.I1", "RecrSz.I")
parNames <- c(parNames, paste(parNames, "lag1", sep="."))

form <- paste("log(Lambda) ~", 
              paste(sapply(parNames, function(term) paste0("s(",term,")")), collapse=" + "))
gamLagLTREmod <- gam(as.formula(form), data=LTREdata)
summary(gamLagLTREmod) # 99.5% of the variance

# model predictions by term
gamLagLTREterms <- predict(gamLagLTREmod, type="terms")
colnames(gamLagLTREterms) <- parNames
# store the associated covariances in a named vector
allCoVar <- calcLTREcontribs(gamLagLTREterms)
# extract the error variance of the model
residVar <- mean(residuals(gamLagLTREmod)^2)

#============================================================================
#
# (III) Plots summarising the sensitivity surfaces -- using the GAM model
#
#============================================================================

keep <- c("Surv.I", "Repr.I", "Recr.I", "RecrSz.I", 
          "GrowJA.I0", "GrowJA.I1", "GrowAJ.I0", "GrowAJ.I1")
keep <- c(keep, paste0(keep, ".lag1"))
pNames <- names(LTREdata)[names(LTREdata) %in% keep]

xRange <- 2.6
nGrid <- 200

# 1. make a data frame with the sensitivity functions

templateDF <- head(LTREdata, nGrid)
templateDF[,] <- 0
for (pNow in pNames) templateDF[[pNow]] <- mean(LTREdata[[pNow]])

parVal <- seq(-xRange/2, +xRange/2, length = nGrid)

sensFunc <- setNames(vector("list", length = length(pNames)), pNames)
for (pNow in pNames) {
  predDF <- templateDF
  predDF[[pNow]] <- mean(LTREdata[[pNow]]) + parVal
  predLam <- predict(gamLagLTREmod, type = "response", newdata = predDF)
  sensFunc[[pNow]] <- data.frame(predLam = predLam, parVal = parVal, parName = pNow)
}
sensFunc <- do.call("rbind", sensFunc)

# 2. make a data frame holding some rug data

isamp <- sample(nrow(LTREdata), 1000)
LTREdataSample <- LTREdata[isamp,]

templateDF <- LTREdataSample
templateDF[,] <- 0
for (pNow in pNames) templateDF[[pNow]] <- mean(LTREdata[[pNow]])

rugData <- setNames(vector("list", length = length(pNames)), pNames)
for (pNow in pNames) {
  predDF <- templateDF
  predDF[[pNow]] <- LTREdataSample[[pNow]]
  predLam <- predict(gamLagLTREmod, type = "response", newdata = predDF)
  parVal <- LTREdataSample[[pNow]] - mean(LTREdata[[pNow]])
  rugData[[pNow]] <- 
    data.frame(predLam = predLam, parVal = parVal, parName = pNow)
}
rugData <- do.call("rbind", rugData)

# 4. just the non-lagged effects...

pKeep1 <- c(`\nSurvival`                 = "Surv.I", 
            `\nReproduction`             = "Repr.I", 
            `\nRecruitment`              = "Recr.I", 
            `\nRecruit Size`             = "RecrSz.I", 
            `Summer Growth\n Juveniles ` = "GrowJA.I0", 
            `Summer Growth\n Older`      = "GrowJA.I1", 
            `Winter Growth\n Juveniles ` = "GrowAJ.I0", 
            `Winter Growth\n Older`      = "GrowAJ.I1")

pltdat1 <- subset(sensFunc, parName %in% pKeep1)
pltdat1$parName <- factor(pltdat1$parName, levels = pKeep1) 
levels(pltdat1$parName) <- names(pKeep1)

pltdat2 <- subset(rugData, parName %in% pKeep1)
pltdat2$parName <- factor(pltdat2$parName, levels = pKeep1) 
levels(pltdat2$parName) <- names(pKeep1)

pp1 <- ggplot(pltdat1, aes(x = parVal, y = predLam, group = parName)) + 
  geom_line() + 
  scale_x_continuous(breaks = c(-1, 0, +1)) +  
  xlab("Centered Parameter Value") + 
  ylab(expression(paste("Log (", italic(lambda)[t], ")"))) + 
  ggtitle("a) Direct Effects") +
  geom_rug(data = pltdat2, col=rgb(.5, 0, 0, alpha=.05)) + 
  facet_wrap(~ parName, nrow = 2, ncol = 4) + 
  theme_bw() + theme(strip.text.x = element_text(size = 11, color="black"),
  plot.title = element_text(color="black", size=14, face="bold.italic"),
  axis.title.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=14, face="bold"),
  axis.text.x = element_text(colour = "black", size=12),
  axis.text.y = element_text(colour = "black", size=12)
)
pp1

# 5. just the lagged effects...

pKeep2 <- paste0(pKeep1, ".lag1")
names(pKeep2) <- names(pKeep1)

pltdat1 <- subset(sensFunc, parName %in% pKeep2)
pltdat1$parName <- factor(pltdat1$parName, levels = pKeep2) 
levels(pltdat1$parName) <- names(pKeep2)

pltdat2 <- subset(rugData, parName %in% pKeep2)
pltdat2$parName <- factor(pltdat2$parName, levels = pKeep2) 
levels(pltdat2$parName) <- names(pKeep2)

pp2 <- ggplot(pltdat1, aes(x = parVal, y = predLam, group = parName)) + 
  geom_line() + 
  scale_x_continuous(breaks = c(-1, 0, +1)) + 
  xlab("Centered Parameter Value") + 
  ylab(expression(paste("Log (", italic(lambda)[t], ")"))) + 
  ggtitle("b) Delayed Effects") +
  geom_rug(data = pltdat2, col=rgb(.5, 0, 0, alpha=.05)) + 
  facet_wrap(~ parName, nrow = 2, ncol = 4) + 
  theme_bw() + theme(strip.text.x = element_text(size = 11, color="black"),
                     plot.title = element_text(color="black", size=14, face="bold.italic"),
                     axis.title.x = element_text(color="black", size=14, face="bold"),
                     axis.title.y = element_text(color="black", size=14, face="bold"),
                     axis.text.x = element_text(colour = "black"),
                     axis.text.y = element_text(colour = "black")
  )
pp2

# 6. both together 
grid.arrange(pp1, pp2, nrow = 2, ncol = 1)

# save multiplot for Figure 1
fig1 <- ggarrange(pp1, pp2,  ncol = 1, nrow = 2, align = "v") #font.label = list(size = 17, color = "black", face = "bold", family = NULL),
#ggsave("fig_1v2.png", plot =  fig1, width=6, height=8, dpi=350)
# end

#============================================================================
#
# (III) Plots summarising the contributions -- using the GAM model
#
#============================================================================

allCoVarSc <- mutate(allCoVar, contrib = 100 * contrib / varLam)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Helper functions
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

label_groups <- function(x) {
  x <- strsplit(x, split = ".", fixed = TRUE)
  x <- sapply(x, function(x) x[1])
  x <- sub("JA", "", x, fixed = TRUE)
  x <- sub("AJ", "", x, fixed = TRUE)
  return(x)
}

split_groups <- function(x, which) {
  x <- strsplit(x, split = ".", fixed = TRUE)
  x <- sapply(x, function(x) x[which])
  return(x)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plots -- 1. sum over the lagged and direct contributions as well as the 
## different sources of growth variation -> not a great plot
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

allCoVarScGrp <- allCoVarSc %>%
  mutate(v1 = label_groups(v1), 
         v2 = label_groups(v2),
         v3 = ifelse(v1 > v2, 
                     paste(v1, v2, sep = "."), 
                     paste(v2, v1, sep = "."))) %>%
  group_by(v3) %>%
  summarise(contrib = sum(contrib, na.rm = TRUE)) %>% 
  mutate(v1 = split_groups(v3, 1), 
         v2 = split_groups(v3, 2),
         posneg = ifelse(contrib > 0, "+", "-")) %>%
  mutate(v1 = ifelse(v3 == "RecrSz.Grow", "Grow",   v1),
         v2 = ifelse(v3 == "RecrSz.Grow", "RecrSz", v2)) %>%
  mutate(v1 = ifelse(v3 == "RecrSz.Recr", "Recr",   v1),
         v2 = ifelse(v3 == "RecrSz.Recr", "RecrSz", v2))
  
pltParNames <- c("Surv", "Repr", "Recr", "Grow", "RecrSz")

allCoVarScGrpPlot <- 
  mutate(allCoVarScGrp, 
         contrib = ifelse(abs(contrib) < 0.1, NA, contrib))

ggplot(allCoVarScGrpPlot,
       aes(x = v1, y = v2, size = abs(contrib), colour = posneg)) +
  geom_point() +  theme_bw() + 
  scale_x_discrete(limits = (pltParNames)) +
  scale_y_discrete(limits = rev(pltParNames)) +
  scale_colour_manual(values = c("steelblue", "tomato"), guide=FALSE) +
  scale_size_area("Contribution", max_size = 20, 
                  limits = c(0, 80), breaks = c(5,10,20,40,80)) + 
  theme(axis.ticks = element_blank(), 
        axis.title = element_blank(),
        axis.text.x = element_text(colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 14),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plots -- 2. 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df1 <- allCoVarSc %>%
  # positive / negative contribution flag
  mutate(posneg  = ifelse(contrib > 0, "+", "-")) %>%
  # keep only the contributions from direct terms
  filter(!v1 %in% v1[grep(".lag1", v1, fixed = TRUE)] &
           !v2 %in% v2[grep(".lag1", v2, fixed = TRUE)]) %>% 
  # remove the lagged label in individual values
  mutate(v1 = sub(".lag1", "", v1, fixed = TRUE),
         v2 = sub(".lag1", "", v2, fixed = TRUE)) %>%
  # now add a lagged / not lagged indicator
  mutate(which = "Direct")

df2 <- allCoVarSc %>%
  # positive / negative contribution flag
  mutate(posneg  = ifelse(contrib > 0, "+", "-")) %>%
  # keep only the contributions from lagged terms
  filter(v1 %in% v1[grep(".lag1", v1, fixed = TRUE)] &
         v2 %in% v2[grep(".lag1", v2, fixed = TRUE)]) %>% 
  # remove the lagged label in individual values
  mutate(v1 = sub(".lag1", "", v1, fixed = TRUE),
         v2 = sub(".lag1", "", v2, fixed = TRUE)) %>%
  # now add a lagged / not lagged indicator
  mutate(which = "Delayed")

####
#### 1. ** Make the bar plot
####

pltData1 <- 
  bind_rows(df1, df2) %>%
  # remove negligable recruitment contribution 
  filter(v1 != "Recr.I" & v2 != "Recr.I") %>%
  # swap a few labels 
  mutate(v1_reduced = 
           ifelse(v1 %in% v1[grep("Grow", v1, fixed = TRUE)], "Grow.I", v1),
         v2_reduced = 
           ifelse(v2 %in% v2[grep("Grow", v2, fixed = TRUE)], "Grow.I", v2)) %>%
  #
  group_by(which, v1_reduced, v2_reduced) %>%
  summarise(contrib = sum(contrib, na.rm = TRUE)) %>% ungroup() %>%
  # remove the NAs
  na.omit %>%
  #
  arrange(desc(abs(contrib)))

# quick look at the combined contributions
# View(pltData1)

# filter combined contributions that are <.5% of var(lambda)
pltData1 <- 
  filter(pltData1, abs(contrib) > 1) %>%
  mutate(term = paste0(which, ": ", v1_reduced, " * ", v2_reduced))

x_labels <- 
  c( `Direct: Surv.I * Surv.I` = "Survival\nVariance",
     `Direct: Repr.I * Repr.I` = "Reproduction\nVariance",
    `Delayed: Repr.I * Repr.I` = "Reproduction\nVariance",
     `Direct: Repr.I * Surv.I` = "Reproduction-Survival\nCovariance",
    `Delayed: Grow.I * Repr.I` = "Growth-Reproduction\nCovariance",
    `Delayed: Grow.I * Grow.I` = "Growth\nVariance")

p1 <- 
  ggplot(pltData1, 
         aes(x = term, y = contrib, fill = which)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual("Source", values = c("thistle4", "wheat3")) + 
  scale_x_discrete(limits=names(x_labels), labels=x_labels) +
  geom_abline(intercept = 0, slope = 0, size = 0.2) + 
  coord_flip() +
  ylab(expression(paste("Contribution to variance of log(", italic(lambda)[t], ")"))) +
  theme_bw() + ggtitle("a)") +
  theme(axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(color="black", size = 14, face="bold"),
        legend.title=element_text(size=12, face="bold"), 
        legend.text=element_text(size=11),
        plot.margin = unit(c(1,0,1,1.5), "cm"))
p1

####
#### 2) ** Make a points plot showing both direct and delayed effects
####

pltData2 <- 
  df2 %>%
  # filter contributions that are > .1% of var(lambda)
  mutate(contrib = ifelse(abs(contrib) < 0.1, NA, contrib)) %>%
  # remove negligable recruitment contribution
  filter(v1 != "Recr.I"   & v2 != "Recr.I" &
         v1 != "Surv.I"   & v2 != "Surv.I" &
         v1 != "RecrSz.I" & v2 != "RecrSz.I") %>%
  # flag to determine which labels to swap
  mutate(
    is_swap = 
      (v1 == "GrowAJ.I1" & v2 == "GrowJA.I0") |
      (v1 == "GrowAJ.I1" & v2 == "GrowAJ.I0") |
      (v1 == "GrowJA.I1" & v2 == "GrowJA.I0")
  ) %>%
  # swap a few labels 
  mutate(
    v1_new = ifelse(is_swap, v2, v1),
    v2_new = ifelse(is_swap, v1, v2)
  ) %>% select(-v1, -v2) %>% rename(v1 = v1_new, v2 = v2_new) %>%
  # remove the NAs
  na.omit

pltParNames <- c("Repr.I", "GrowJA.I0", "GrowJA.I1", "GrowAJ.I0", "GrowAJ.I1")

x_labels <- 
  c(`GrowJA.I0` = "Summer Growth\nJuveniles",
    `GrowJA.I1` = "Summer Growth\nOlder",
    `GrowAJ.I0` = "Winter Growth\nJuveniles",
    `GrowAJ.I1` = "Winter Growth\nOlder",
    `Repr.I`    = "Reproduction")

x_labels <- rev(x_labels)

p2 <- 
  ggplot(pltData2, 
       aes(x = v1, y = v2, size = abs(contrib), colour = posneg)) +
  geom_point() +  
  scale_x_discrete(limits=rev(names(x_labels)), labels=rev(x_labels)) +
  scale_y_discrete(limits=names(x_labels), labels=x_labels) +
  scale_colour_manual("Direction", values = c("steelblue", "tomato")) +
  scale_size_area("Contribution (%)", max_size = 18, 
                  limits = c(0, 10), breaks = c(1,2,4,8)) + 
  theme_bw() + ggtitle("b)") +
  theme(axis.ticks = element_blank(), 
        axis.title = element_blank(),
        axis.text.x = element_text(size = 12, angle = 50, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, hjust = 1, color = "black"),
        legend.title=element_text(size=12, face="bold"), 
        legend.text=element_text(size=11),
        plot.margin = unit(c(1,0,0,1.5), "cm"))
p2


## save multiplot for Figure 2. 
# 1 column 2 rows
width = 6.5
height = (1 + 9/16) * width
# ggsave("fig_2v2b.png", plot = grid.draw(cbind(rbind(ggplotGrob(p1), ggplotGrob(p2), size="last"), size='last')), width= width, height= height, dpi=350)
# ggsave("fig_2v2a.png", plot = grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size="last")),
#        width= width, height= height, dpi=350)

# 2 columns 1 row
fig2 <- grid.arrange(p1, p2, nrow = 1, ncol = 2, respect=TRUE) #widths = c(1.2, 1),
width = 12
height = 9/16 * width
#ggsave("fig_2v2.png", plot =  fig2, width= width, height= height, dpi=350)

# end

####
#### 3) ** Make a points plot showing both direct and delayed effects 
####       (not in paper)
####

pltData3 <- 
  bind_rows(df1, df2) %>%
  # filter contributions that are > 1% of var(lambda)
  mutate(contrib = ifelse(abs(contrib) < 1, NA, contrib)) %>%
  # remove negligable recruitment contribution
  filter(v1 != "Recr.I" & v2 != "Recr.I") %>%
  # remove the NAs
  na.omit

pltParNames <- c("Surv.I", "Repr.I",
                 "GrowJA.I0", "GrowJA.I1", "GrowAJ.I0", "GrowAJ.I1", "RecrSz.I")

# 
ggplot(pltData3, 
       aes(x = v1, y = v2, size = abs(contrib), colour = posneg)) +
  geom_point() +  
  scale_x_discrete(limits = rev(pltParNames)) +
  scale_y_discrete(limits = (pltParNames)) +
  scale_colour_manual(values = c("steelblue", "tomato"), guide=FALSE) +
  scale_size_area("Contribution", max_size = 8, 
                  limits = c(0, 80), breaks = c(5,10,20,40, 80)) + 
  facet_wrap(~ which, nrow = 1) +
  theme_bw() + 
  theme(axis.ticks = element_blank(), 
        strip.text.x = element_text(size = 8, color="black"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 55, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"))


#============================================================================
#
# (IV) Environmental LTRE analysis 
#
#============================================================================

envLTREdata <- within(envSet, {
  eNamesSuffix <- colnames(eStateSet)
  
  eNames <- paste("GrowAJ", eNamesSuffix, sep=".")
  GrowAJ <- cbind(I0 = mParSet[, "GrowAJ.I"], t(mParSet[1, eNames] * t(eStateSet)))
  GrowAJ <- cbind(I1 = mParSet[,"GrowAJ.Idiff"] + GrowAJ[,"I0"], GrowAJ)
  
  eNames <- paste("GrowJA", eNamesSuffix, sep=".")
  GrowJA <- cbind(I0 = mParSet[, "GrowJA.I"], t(mParSet[1, eNames] * t(eStateSet)))
  GrowJA <- cbind(I1 = mParSet[,"GrowJA.Idiff"] + GrowJA[,"I0"], GrowJA)
  
  eNames <- paste("RecrSz", eNamesSuffix, sep=".")
  RecrSz <- cbind(I = mParSet[, "RecrSz.I"], t(mParSet[1, eNames] * t(eStateSet)))
  
  eNames <- paste("Surv", eNamesSuffix, sep=".")
  Surv <- cbind(I = mParSet[, "Surv.I"], t(mParSet[1, eNames] * t(eStateSet)))
  
  eNames <- paste("Repr", eNamesSuffix, sep=".")
  Repr <- cbind(I = mParSet[, "Repr.I"], t(mParSet[1, eNames] * t(eStateSet)))
  
  eNames <- paste("Recr", eNamesSuffix, sep=".")
  Recr <- cbind(I = mParSet[, "Recr.I"], t(mParSet[1, eNames] * t(eStateSet)))
  
  rm(eNamesSuffix, eNames)
})

# clean up
envLTREdata$mParSet <- NULL; envLTREdata$eStateSet <- NULL
envLTREdata <- suppressWarnings(as.data.frame(envLTREdata))
envLTREdata <- data.frame(envLTREdata, Lambda = simOut$lamt)

# add lagged variables
laggedData <- head(subset(envLTREdata, select = -Lambda), -1)
names(laggedData) <- paste(names(laggedData), "lag1", sep=".")
envLTREdata <- cbind(tail(envLTREdata, -1), laggedData)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compute all the contributions
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

k_vals  <- c(0, 1)
ij_vals <- c("Surv", "Repr", "GrowJA.0", "GrowJA.1", "GrowAJ.0", "GrowAJ.1")
pq_vals <- c("SprT", "WinT", "BrGd")

beta_name <- function(k, ij) {
  nm <- ifelse(grepl("[.]", ij), sub("[.]", ".I", ij), paste0(ij, ".I"))
  nm <- ifelse(k, paste0(nm, ".lag", k), nm)
  return(nm)
}

term_name <- function(k, ij, pq) {
  nm <- sub("[.][[:alnum:]]", "", ij)
  nm <- paste0(nm, ".", pq)
  nm <- ifelse(k, paste0(nm, ".lag", k), nm)
  return(nm)
}

calc_cov <- function(term1, term2) {
  cov(envLTREdata[[term1]], envLTREdata[[term2]])  
}

##
## 1. refit the linear-lagged LTRE model w/ only the terms that matter
##

parNames <- c("Surv.I", "Repr.I", 
              "GrowAJ.I0.lag1", "GrowJA.I0.lag1", "GrowAJ.I1.lag1", "GrowJA.I1.lag1",
              "RecrSz.I.lag1", "Repr.I.lag1")

form <- paste("log(Lambda) ~", paste(parNames, collapse=" + "))
envLTREmod <- lm(as.formula(form), data = LTREdata)
betas <- coef(envLTREmod)[-1]

##
## 2. Calculate the contributions we can attribute to the env covars
## 

calc_cov <- function(term1, term2) {
  cov(envLTREdata[[term1]], envLTREdata[[term2]])  
}

all_envir_terms <- 
  expand.grid(k = k_vals, 
              i = ij_vals, j = ij_vals,
              p = pq_vals, q = pq_vals, 
              stringsAsFactors = FALSE) %>%
  mutate(
    # construct the term names to reference in 'envLTREdata'
    term1 = term_name(k, i, p), term2 = term_name(k, j, q),
    # extract the corresponding covariance for the terms
    covr_contrib = mapply(calc_cov, term1, term2),
    # work out the sensitivities to scale things by
    beta_contrib = betas[beta_name(k, i)] * betas[beta_name(k, j)],
    # compute the total contribution relative to total variance
    total_contrib = 100 * covr_contrib * beta_contrib / var(log(envLTREdata$Lambda))
  ) %>% na.omit()

# quick look over the output
arrange(all_envir_terms, desc(abs(total_contrib))) %>% View

##
## 3. Calculate the remaining unexplained contributions
## 

names(envLTREdata)

all_other_terms <- 
  expand.grid(k = k_vals, 
              i = ij_vals, j = ij_vals,
              stringsAsFactors = FALSE) %>%
  mutate(
    # construct the term names to reference in 'envLTREdata'
    term1 = beta_name(k, i), term2 = beta_name(k, j),
    # extract the corresponding covariance for the terms
    covr_contrib = mapply(calc_cov, term1, term2),
    # work out the sensitivities to scale things by
    beta_contrib = betas[term1] * betas[term2],
    # compute the total contribution relative to total variance
    total_contrib = 100 * covr_contrib * beta_contrib / var(log(envLTREdata$Lambda))
  ) %>% na.omit()

# quick look over the output
arrange(all_other_terms, desc(abs(total_contrib))) %>% View

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Make the summary plots
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##
## 1. Clean up the two sets of contributions
##

summary_envir_terms <-
  all_envir_terms %>%
  mutate(    
    # remove the `.lag1` from term labels
    term1 = sub(".lag1", "", term1, fixed = TRUE),
    term2 = sub(".lag1", "", term2, fixed = TRUE),
    # make a unique term using alphabetical order
    term = ifelse(term1 < term2, 
                  paste0(term1, "-", term2), paste0(term2, "-", term1))
  ) %>% 
  # now sum over the two parts of the covar contribs
  rename(is_lag = k) %>%
  group_by(term, is_lag) %>%
  summarise(total_contrib = sum(total_contrib)) %>% 
  mutate(
    # split the label back up again 
    term1 = strsplit(term, "[.-]")[[1]][1],
    term2 = strsplit(term, "[.-]")[[1]][3],
    env1  = strsplit(term, "[.-]")[[1]][2],
    env2  = strsplit(term, "[.-]")[[1]][4]
  ) %>% 
  arrange(desc(abs(total_contrib)))

## 

summary_other_terms <-
  all_other_terms %>%
  mutate(    
    # remove the `.lag1` from term labels
    term1 = sub(".lag1", "", term1, fixed = TRUE),
    term2 = sub(".lag1", "", term2, fixed = TRUE),
    # make a unique term using alphabetical order
    term = ifelse(term1 < term2, 
                  paste0(term1, "-", term2), paste0(term2, "-", term1))
  ) %>% 
  # now sum over the two parts of the covar contribs
  rename(is_lag = k) %>%
  group_by(term, is_lag) %>%
  summarise(total_contrib = sum(total_contrib)) %>% 
  mutate(
    # split the label back up again 
    term1 = strsplit(term, "[.-]")[[1]][1],
    term2 = strsplit(term, "[.-]")[[1]][3]
  ) %>% 
  arrange(desc(abs(total_contrib)))

## quick look

View(summary_envir_terms)
View(summary_other_terms)

##
## 2. vectors to label ggplot aesthetics
## 

x_labels <- 
  c(`Surv-Surv-Lag0` = "Survival\n(Direct)",
    `Repr-Repr-Lag0` = "Reproduction\n(Direct)",
    `Repr-Repr-Lag1` = "Reproduction\n(Delayed)",
    `Repr-Surv-Lag0` = "Reproduction-Survival\n(Direct)",
    `Grow-Repr-Lag1` = "Growth-Reproduction\n(Delayed)",
    `Grow-Grow-Lag1` = "Growth\n(Delayed)")

fill_labels_env <- 
  c(`BrGd-BrGd` = "Var(SF)",
    `BrGd-SprT` = bquote(paste("Cov(SF, ", T[spring], ")")),
    `BrGd-WinT` = bquote(paste("Cov(SF, ", T[winter], ")")),
    `SprT-SprT` = bquote(paste("Var(", T[spring], ")")),
    `SprT-WinT` = bquote(paste("Cov(", T[spring], ", ",T[winter],")")),
    `WinT-WinT` = bquote(paste("Var(", T[winter], ")")))

fill_labels_src <- 
  c(`envir`   = "Environment", `other` = "Unexplained")

##
## 3. Prepare dataset to plot environmental contributions
##

plt_data <- list()

plt_data[[1]] <- 
  summary_envir_terms %>%
  mutate(
    term1 = sub("JA|AJ", "", term1),
    term2 = sub("JA|AJ", "", term2)
  ) %>%
  group_by(is_lag, term1, term2) %>%
  summarise(total_contrib = sum(total_contrib)) %>%
  mutate(source = "envir")

plt_data[[2]] <- 
  summary_other_terms %>%
  mutate(
    term1 = sub("JA|AJ", "", term1),
    term2 = sub("JA|AJ", "", term2)
  ) %>%
  group_by(is_lag, term1, term2) %>%
  summarise(total_contrib = sum(total_contrib)) %>%
  mutate(source = "other")

plt_data <- 
  bind_rows(plt_data) %>%
  mutate(term = paste0(term1, "-", term2, "-Lag", is_lag))

##
## Plot environment summary
##
ppp1 <- 
  ggplot(plt_data, 
       aes(x = term, y = total_contrib, fill = source)) +
  geom_bar(stat = "identity", width = 0.7, position = "dodge") +
  scale_x_discrete(limits=names(x_labels), labels=x_labels) +
  scale_fill_discrete(limits=names(fill_labels_src), labels=fill_labels_src) +
  geom_abline(intercept = 0, slope = 0, size = 0.2) + 
  coord_flip() +
  ylab(expression(paste("Contribution to variance of log(", italic(lambda)[t], ")"))) +
  labs(fill = "Source of variation") +
  ggtitle("a)") + 
  theme_light() + 
  theme(axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 12, color = " black"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, color = " black"),
        legend.title=element_text(size=12, face="bold"), 
        legend.text=element_text(size=11),
        axis.title.x = element_text(size=14,face="bold"))

ppp1

##
## Prepare data to plot environmental sources   
##

plt_data <- 
  summary_envir_terms %>% 
  ungroup() %>%
  mutate(
    # need to make sure same order is used
    env3 = ifelse(env1 > env2, env1, env2),
    env1 = ifelse(env1 > env2, env2, env1),
    env2 = ifelse(env3 > env2, env3, env2),
    #
    term1 = sub("JA|AJ", "", term1),
    term2 = sub("JA|AJ", "", term2),
    #
    term = paste0(term1, "-", term2, "-Lag", is_lag)
  ) %>% select(-env3) %>%
  group_by(is_lag, term, env1, env2) %>%
  summarise(total_contrib = sum(total_contrib)) %>%
  arrange(desc(abs(total_contrib))) 

sum(abs(plt_data$total_contrib))
sum(plt_data$total_contrib)

sum(filter(plt_data, env1 == "WinT" | env2 == "WinT")$total_contrib)

## 
## Plot environemntal sources
##

ppp2 <- 
  ggplot(plt_data, 
       aes(x = term, y = total_contrib, fill = paste0(env1,"-",env2))) +
  geom_bar(stat = "identity", width = 0.7, position = "dodge") +
  scale_x_discrete(limits=names(x_labels), labels=x_labels) +
  scale_fill_discrete(limits=names(fill_labels_env), labels=fill_labels_env) +
  geom_abline(intercept = 0, slope = 0, size = 0.2) + 
  coord_flip() +
  ylab(expression(paste("Contribution to variance of log(", italic(lambda)[t], ")"))) +
  labs(fill = "Source of variation") +
  ggtitle("b)") + 
  theme_light() + 
  theme(axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.y = element_blank(),
        legend.title=element_text(size=12, face="bold"), 
        legend.text=element_text(size=11),
        axis.title.x = element_text(size=14,face="bold"))

ppp2

## both together
grid.arrange(ppp1, ppp2, nrow = 2, ncol = 1)

# save multiplot for Figure 3
fig3 <- ggarrange(ppp1, ppp2, font.label = list(size = 17, color = "black", face = "bold", family = NULL), ncol = 1, nrow = 2, align = "v")
# ggsave("fig_3.png", plot =  fig3, width=8, height=10, dpi=350)
# end







