### Load Packages
lapply(c("dplyr", "RMark", "plyr", "janitor", "chron", "ggplot2", "patchwork"), require, character.only = T)


### Load AIC results
uni.aic <- read.csv("CapMort - Univariate w Interaction AIC.csv")
all.aic <- read.csv("CapMort - AllComboAIC.csv")
threshold.aic <- read.csv("EWT.capmort.threshold.results.csv")
# Generate List of Covariates
cov.names <- colnames(tcm.data.raw[,7:ncol(tcm.data.raw)])
uni.mods <- paste0(cov.names, " * LN")

### Extract coefficient estimates from best models. Compare to univariate/original approach
#For mutate_at quick fix
scale.fix <- function(x){hold <- scale(x)
return(hold[,1])}

# Load original data
tcm.data.raw <- read.csv(file = "TurkCapMort_EH.csv") %>%
  mutate_at(c('Study.Area'), .funs = function(x){ gsub(" ", "", x)}) %>%
  mutate_at(c('Study.Area', 'Location', 'Sex', 'SexTT', 'TurkAge', 'Trans.Type',
              'Hematoma', 'Year', 'SampType', 'MG', 'LPDV', 'REV', 'Group'), as.factor) %>%
  select(-Location, -SexTT, -Pat.Tag, -Weight) %>%
  mutate_at(c('YearDay', 'HandTime', 'Prcp.Day', 'Prcp.Week',
              'Tavg.Day', 'Tavg.Week', 'BodyCon'), scale.fix)

# Prep RMark Object for re-running models
capmort.process = process.data(tcm.data.raw,model="Nest",nocc=30,
                               groups=c('Study.Area', 'Sex', 'TurkAge', 'Trans.Type',
                                        'Hematoma', 'Year', 'SampType', 'MG', 'LPDV', 'REV', 'Group'))
capmort.ddl = make.design.data(capmort.process)
#Log-transform of date effect#
# First create 'group covariate' for log-transformation of time
Daypost<- data.frame(seq(1,30,1))
colnames(Daypost)<- c("time")
Daypost$LN<-log(Daypost$time)
# Add LN date data to design data
capmort.ddl$S=merge_design.covariates(capmort.ddl$S,Daypost,bygroup=FALSE, bytime=TRUE)


#Find top model for supported univariates
supported.covs <- uni.aic %>%
  filter(AICc < .$AICc[which(npar == 2)]) %>%
  mutate(model = gsub(" \\* LN\\)", "", model)) %>%
  mutate(model = gsub("S\\(\\~", "", model)) %>%
  pull(model)

top.for.cov <- all.aic[sapply(supported.covs, FUN = function(x){min(which(grepl(x, all.aic$S)))}),] %>%
  mutate(Cov = supported.covs) %>%
  select(Cov, S, AICc)







##########################################################################
### Model Weights by parameter
all.aic$S

cov.model <- sapply(gsub(" \\* LN", "", uni.mods), FUN = function(x){
  ifelse(grepl(x, ifelse(grepl("YearDay", all.aic$S), "YDay", all.aic$S)), 1, 0)
})

aic.weights.by.cov <- cbind(cov.model, all.aic) %>%
  mutate(DeltaAICc = AICc - first(AICc)) %>%
  mutate(weight = exp(-0.5*DeltaAICc)) %>%
  mutate(weight = weight/sum(weight)) %>%
  mutate_at(gsub(" \\* LN", "", uni.mods), ~ . * weight) %>%
  select(gsub(" \\* LN", "", uni.mods)) %>%
  colSums() %>% sort()
