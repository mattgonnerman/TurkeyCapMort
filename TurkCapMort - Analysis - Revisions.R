#Load Packages
lapply(c("dplyr", "RMark", "plyr", "janitor", "chron", "ggplot2", 'stringr'), require, character.only = T)

scale.fix <- function(x){hold <- scale(x)
return(hold[,1])}


#############################################################################################
### A Priori Models
tcm.data.raw <- read.csv(file = "TurkCapMort_EH.csv") %>%
  mutate_at(c('Study.Area'), .funs = function(x){ gsub(" ", "", x)}) %>%
  mutate_at(c('Study.Area', 'Location', 'Sex', 'SexTT', 'TurkAge', 'Trans.Type',
              'Hematoma', 'Year', 'SampType', 'MG', 'LPDV', 'REV', 'Group'), as.factor) %>%
  select(-Location, -SexTT, -Pat.Tag, -Weight) %>%
  mutate_at(c('YearDay', 'HandTime', 'Prcp.Day', 'Prcp.Week',
              'Tavg.Day', 'Tavg.Week', 'BodyCon'), scale.fix)
  

### Prepare for RMark
capmort.process = process.data(tcm.data.raw,model="Nest",nocc=30,
                               groups=c('Study.Area', 'Sex', 'TurkAge', 'Trans.Type',
                                        'Hematoma', 'Year', 'SampType', 'MG', 'LPDV', 'REV', 'Group'))
capmort.ddl = make.design.data(capmort.process)

###Log-transform of date effect###
## First create 'group covariate' for log-transformation of time
Daypost<- data.frame(seq(1,30,1))
colnames(Daypost)<- c("time")
Daypost$LN<-log(Daypost$time)

## Add LN date data to design data
capmort.ddl$S=merge_design.covariates(capmort.ddl$S,Daypost,bygroup=FALSE, bytime=TRUE)

#Generate List of Models to run
cov.names <- colnames(tcm.data.raw[,7:(ncol(tcm.data.raw)-2)])
allcomb.1.mods <- c("1", "LN", paste0("LN + ", cov.names))
allcomb.2.mods <- paste0("LN + ", combn(cov.names, 2, FUN = paste, collapse = ' + '))
allcomb.3.mods <- paste0("LN + ", combn(cov.names, 3, FUN = paste, collapse = ' + '))
allcomb.mods <- c(allcomb.1.mods, allcomb.2.mods, allcomb.3.mods)
models.df <- data.frame(ID = 1:length(allcomb.mods), Model = allcomb.mods)

setwd("./Outputs/")
run.all.combo <- function(){
  for(x in 1:nrow(models.df)){
    assign(paste0("S.Model", x),
           list(formula = as.formula(paste0("~ ", models.df$Model[x]))))
  }
  
  capmort.model.list = create.model.list("Nest")
  
  capmort.results = mark.wrapper(capmort.model.list, data = capmort.process, ddl = capmort.ddl,
                                 output=FALSE,  invisible = TRUE)
  
  # Return model table and list of models
  return(capmort.results)
}

apriori.results <- run.all.combo()
apriori.results

setwd("./../")
all.aic <- apriori.results$model.table
write.csv(all.aic, file = "CapMort - AllComboAIC.csv", row.names=FALSE)

uni.aic <- as.data.frame(all.aic) %>%
  filter(str_count(S, "\\+") < 2) %>% arrange(AICc) %>%
  mutate(DeltaAICc = AICc - first(AICc)) %>% 
  select(model, npar, AICc, DeltaAICc, Deviance) %>%
  mutate(Weight = exp(-0.5*DeltaAICc)) %>%
  mutate(Weight = Weight/sum(Weight)) 

write.csv(uni.aic, file = "CapMort - UNIAIC.csv", row.names=FALSE)
  

#####################################################################################################################################
### Ad Hoc Models

## TIME PERIOD MODELS
setwd("./Outputs/")
S.Year = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ Year * LN, link="loglog")))
S.Feb19 = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ BefFeb19 * LN, link="loglog")))
S.Mar19 = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ BefMar19 * LN, link="loglog")))
setwd("./../")

S.Year$results$beta

## TRANSMITTER TYPE vs SEX
tcm.data.Fonly <- tcm.data.raw %>% filter(Sex == "F")
tcm.data.Neckonly <- tcm.data.raw %>% filter(Trans.Type == "Neck")
capmort.process.Fonly = process.data(tcm.data.Fonly, model="Nest",nocc=30,
                               groups=c("Sex", "Trans.Type"))
capmort.process.Neckonly = process.data(tcm.data.Neckonly, model="Nest",nocc=30,
                                     groups=c("Sex", "Trans.Type"))
capmort.ddl.Fonly = make.design.data(capmort.process.Fonly)
capmort.ddl.Neckonly = make.design.data(capmort.process.Neckonly)
capmort.ddl.Fonly$S=merge_design.covariates(capmort.ddl.Fonly$S,Daypost,bygroup=FALSE, bytime=TRUE)
capmort.ddl.Neckonly$S=merge_design.covariates(capmort.ddl.Neckonly$S,Daypost,bygroup=FALSE, bytime=TRUE)
setwd("./Outputs/")
S.TT.FOnly = mark(capmort.process.Fonly, capmort.ddl.Fonly, model.parameters=list(S=list(formula =~ Trans.Type * LN, link="loglog")))
S.Sex.NeckOnly = mark(capmort.process.Neckonly, capmort.ddl.Neckonly, model.parameters=list(S=list(formula =~ Sex * LN, link="loglog")))
setwd("./../")

## HANDLING TIME vs GROUP/SOLO
tcm.data.Group <- tcm.data.raw %>% filter(Group == 1)
tcm.data.Solo <- tcm.data.raw %>% filter(Group == 0)
capmort.process.Group = process.data(tcm.data.Group, model="Nest",nocc=30)
capmort.process.Solo = process.data(tcm.raw.Solo, model="Nest",nocc=30)
capmort.ddl.Group = make.design.data(capmort.process.Group)
capmort.ddl.Solo = make.design.data(capmort.process.Solo)
capmort.ddl.Group$S=merge_design.covariates(capmort.ddl.Group$S,Daypost,bygroup=FALSE, bytime=TRUE)
capmort.ddl.Solo$S=merge_design.covariates(capmort.ddl.Solo$S,Daypost,bygroup=FALSE, bytime=TRUE)

S.HandTime.Group = mark(capmort.process.Group, capmort.ddl.Group, model.parameters=list(S=list(formula =~ HandTime * LN, link="loglog")))
S.HandTimeSolo = mark(capmort.process.Solo, capmort.ddl.Solo, model.parameters=list(S=list(formula =~ HandTime * LN, link="loglog")))

#####################################################################################################################################
###Identify threshold date for capture mortality
#Use best performing model
thresh.modform <- apriori.results$model.table$S[1]
setwd("./Outputs/")
# Add Threshold model for individual days
lapply(2:29, function(x){
  assign(paste0("capmort", x, ".ddl"),
         add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,x,30), name=paste0("Days",x-1), right=F, replace=T),
         envir = .GlobalEnv)
  })

#Run Models
rm(list=ls(pattern="^S."))
EWT.capmort.threshold.models=function(){
  lapply(1:28, function(x){
    assign(paste0('S.Days',x),
           mark(capmort.process, get(paste0("capmort", x + 1, ".ddl")), 
                model.parameters=list(
                  S=list(
                    formula = paste0(thresh.modform, "+ Days", x), link="loglog"))),
           envir = .GlobalEnv)
  })
  
  return()
} 
EWT.capmort.threshold.models()
EWT.capmort.threshold.results <- collect.models()
setwd("./../")


write.csv(EWT.capmort.threshold.results$model.table, file = "ThresholdAIC.csv", row.names=FALSE)
