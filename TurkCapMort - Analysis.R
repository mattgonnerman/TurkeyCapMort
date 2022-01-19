#Load Packages
lapply(c("dplyr", "RMark", "plyr", "janitor", "chron", "ggplot2"), require, character.only = T)

#############################################################################################
### A Priori Models
EWT.EH.cov <- read.csv(file = "CapMortEncounterHistory_withcovariatesBC.csv") %>%
  mutate(Study.Area = as.factor(Study.Area),
         Location = as.factor(Location),
         Sex = as.factor(Sex),
         SexTT = as.factor(SexTT),
         TurkAge = as.factor(TurkAge),
         Trans.Type = as.factor(Trans.Type),
         Hematoma = as.factor(Hematoma),
         Pat.Tag = as.factor(Pat.Tag),
         Year = as.factor(Year))

capmort.process = process.data(EWT.EH.cov,model="Nest",nocc=30,
                               groups=c("Study.Area", "Location", "Sex", "TurkAge", "Trans.Type", "SexTT", "Hematoma", "Pat.Tag", "Year"))
capmort.ddl = make.design.data(capmort.process)

###Log-transform of date effect###
## First create 'group covariate' for log-transformation of time
Daypost<- data.frame(seq(1,30,1))
colnames(Daypost)<- c("time")
Daypost$LN<-log(Daypost$time)

## Add LN date data to design data
capmort.ddl$S=merge_design.covariates(capmort.ddl$S,Daypost,bygroup=FALSE, bytime=TRUE)

#Create Model List
EWT.capmort.models=function()
{
  #Time Models
  S.dot = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ 1, link="loglog")))
  S.LN = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ LN, link="loglog")))
  S.time = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ Time, link="loglog")))
  S.time2 = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ Time + I(Time^2), link="loglog")))
  
  #Capture and Handling Models
  S.turkage = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ TurkAge * LN, link="loglog")))
  S.sex = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ Sex * LN, link="loglog")))
  S.transtype = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ Trans.Type * LN, link="loglog")))
  S.TTSex = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ SexTT * LN, link="loglog")))
  S.studyarea = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ Study.Area * LN, link="loglog")))
  S.Location = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ Location * LN, link="loglog"))) #This will likely be explained better by weather covariates
  S.LPDV = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ LPDV * LN, link="loglog")))
  S.MG = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ MG * LN, link="loglog")))
  S.REV = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ REV * LN, link="loglog")))
  S.Hematoma = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ Hematoma * LN, link="loglog")))
  S.DayYear = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ YearDay * LN, link="loglog")))
  S.HandlingTime = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ HandTime * LN, link="loglog")))
  S.BodyCondition = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ BodyCon * LN, link="loglog")))
  S.OrdDay = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ YearDay * LN, link="loglog")))
  
  #Weather Models
  S.mt.CapDate = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ mtCDy * LN, link="loglog")))
  S.avgmt.Week = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ avgmtWk * LN, link="loglog")))
  S.prec.CapDate = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ pcpCDy * LN, link="loglog")))
  S.avgprec.Week = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ avgpcpWk * LN, link="loglog")))
  
  return(collect.models() )
}

EWT.capmort.results = EWT.capmort.models()
EWT.capmort.results

write.csv(EWT.capmort.results$model.table, file = "EWT.capmort.results.csv", row.names=FALSE)

#Save Individual Model Outputs
#Export summaries as single file
sink("CapMort_APrioriModels_BetaEst.csv") #At Home
cat("exp(-exp(-x))")
cat('\n')
cat("For interpretting the beta use reverse of loglog")
cat('\n')
cat('\n')
cat("Model Estimates")
cat('\n')
cat('\n')
cat("Handling Time")
cat('\n')
write.csv(EWT.capmort.results$S.HandlingTime$results$beta)
cat('______________________________________')
cat('\n')
cat('\n')
cat("Transmitter Type")
cat('\n')
write.csv(EWT.capmort.results$S.transtype$results$beta)
cat('______________________________________')
cat('\n')
cat('\n')
cat("Turkey Age")
cat('\n')
write.csv(EWT.capmort.results$S.turkage$results$beta)
cat('______________________________________')
cat('\n')
cat('\n')
cat("REV")
cat('\n')
write.csv(EWT.capmort.results$S.REV$results$beta)
cat('______________________________________')
cat('\n')
cat('\n')
cat("Sex")
cat('\n')
write.csv(EWT.capmort.results$S.Sex$results$beta)
cat('______________________________________')
cat('\n')
cat('\n')
sink()

#####################################################################################################################################
### Ad Hoc Models

## YEAR
S.Year = mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ Year * LN, link="loglog")))
S.Year$results$beta

## TRANSMITTER TYPE
EWT.EH.cov.TT <- EWT.EH.cov %>% filter(Sex == "F")

capmort.process.TT = process.data(EWT.EH.cov.TT,model="Nest",nocc=30,
                               groups=c("Study.Area", "Location", "Sex", "TurkAge", "Trans.Type", "Hematoma", "Pat.Tag", "Year"))
capmort.ddl.TT = make.design.data(capmort.process.TT)

###Log-transform of date effect###
## First create 'group covariate' for log-transformation of time
Daypost<- data.frame(seq(1,30,1))
colnames(Daypost)<- c("time")
Daypost$LN<-log(Daypost$time)

## Add LN date data to design data
capmort.ddl.TT$S=merge_design.covariates(capmort.ddl.TT$S,Daypost,bygroup=FALSE, bytime=TRUE)

S.TT.solo = mark(capmort.process.TT, capmort.ddl.TT, model.parameters=list(S=list(formula =~ Trans.Type * LN, link="loglog")))
S.TT.solo$results$beta

## SEX
EWT.EH.cov.Sex <- EWT.EH.cov %>% filter(Trans.Type == "Neck")

capmort.process.Sex = process.data(EWT.EH.cov.Sex,model="Nest",nocc=30,
                                  groups=c("Study.Area", "Location", "Sex", "TurkAge", "Trans.Type", "Hematoma", "Pat.Tag", "Year"))
capmort.ddl.Sex = make.design.data(capmort.process.Sex)

###Log-transform of date effect###
## First create 'group covariate' for log-transformation of time
Daypost<- data.frame(seq(1,30,1))
colnames(Daypost)<- c("time")
Daypost$LN<-log(Daypost$time)

## Add LN date data to design data
capmort.ddl.Sex$S=merge_design.covariates(capmort.ddl.Sex$S,Daypost,bygroup=FALSE, bytime=TRUE)

S.Sex.solo = mark(capmort.process.Sex, capmort.ddl.Sex, model.parameters=list(S=list(formula =~ Sex * LN, link="loglog")))
S.Sex.solo$results$beta

sink("CapMort_AdHocModels_CoefEst.csv") #At Home
cat("exp(-exp(-x))")
cat('\n')
cat("For interpretting the beta use reverse of loglog")
cat('\n')
cat('\n')
cat("Model Estimates")
cat('\n')
cat('\n')
cat("Year")
cat('\n')
write.csv(S.Year$results$beta)
cat('______________________________________')
cat('\n')
cat('\n')
cat("Transmitter Type")
cat('\n')
write.csv(S.TT.solo$results$beta)
cat('______________________________________')
cat('\n')
cat('\n')
cat("Turkey Age")
cat('\n')
write.csv(S.Sex.solo$results$beta)
sink()

#####################################################################################################################################
#Identify threshold date for capture mortality

### Add Threshold model for individual days
capmort2.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,2,30), name="Days1", right=F, replace=T)
capmort3.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,3,30), name="Days2", right=F, replace=T)
capmort4.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,4,30), name="Days3", right=F, replace=T)
capmort5.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,5,30), name="Days4", right=F, replace=T)
capmort6.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,6,30), name="Days5", right=F, replace=T)
capmort7.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,7,30), name="Days6", right=F, replace=T)
capmort8.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,8,30), name="Days7", right=F, replace=T)
capmort9.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,9,30), name="Days8", right=F, replace=T)
capmort10.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,10,30), name="Days9", right=F, replace=T)
capmort11.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,11,30), name="Days10", right=F, replace=T)
capmort12.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,12,30), name="Days11", right=F, replace=T)
capmort13.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,13,30), name="Days12", right=F, replace=T)
capmort14.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,14,30), name="Days13", right=F, replace=T)
capmort15.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,15,30), name="Days14", right=F, replace=T)
capmort16.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,16,30), name="Days15", right=F, replace=T)
capmort17.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,17,30), name="Days16", right=F, replace=T)
capmort18.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,18,30), name="Days17", right=F, replace=T)
capmort19.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,19,30), name="Days18", right=F, replace=T)
capmort20.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,20,30), name="Days19", right=F, replace=T)
capmort21.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,21,30), name="Days20", right=F, replace=T)
capmort22.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,22,30), name="Days21", right=F, replace=T)
capmort23.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,23,30), name="Days22", right=F, replace=T)
capmort24.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,24,30), name="Days23", right=F, replace=T)
capmort25.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,25,30), name="Days24", right=F, replace=T)
capmort26.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,26,30), name="Days25", right=F, replace=T)
capmort27.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,27,30), name="Days26", right=F, replace=T)
capmort28.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,28,30), name="Days27", right=F, replace=T)
capmort29.ddl=add.design.data(capmort.process, capmort.ddl, parameter="S", type="time",bins=c(1,29,30), name="Days28", right=F, replace=T)

EWT.capmort.threshold.models=function()
{
  S.Days1A = mark(capmort.process, capmort2.ddl, model.parameters=list(S=list(formula =~ Days1 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days2A = mark(capmort.process, capmort3.ddl, model.parameters=list(S=list(formula =~ Days2 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days3A = mark(capmort.process, capmort4.ddl, model.parameters=list(S=list(formula =~ Days3 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days4A = mark(capmort.process, capmort5.ddl, model.parameters=list(S=list(formula =~ Days4 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days5A = mark(capmort.process, capmort6.ddl, model.parameters=list(S=list(formula =~ Days5 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days6A = mark(capmort.process, capmort7.ddl, model.parameters=list(S=list(formula =~ Days6 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days7A = mark(capmort.process, capmort8.ddl, model.parameters=list(S=list(formula =~ Days7 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days8A = mark(capmort.process, capmort9.ddl, model.parameters=list(S=list(formula =~ Days8 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days9A = mark(capmort.process, capmort10.ddl, model.parameters=list(S=list(formula =~ Days9 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days10A = mark(capmort.process, capmort11.ddl, model.parameters=list(S=list(formula =~ Days10 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days11A = mark(capmort.process, capmort12.ddl, model.parameters=list(S=list(formula =~ Days11 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days12A = mark(capmort.process, capmort13.ddl, model.parameters=list(S=list(formula =~ Days12 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days13A = mark(capmort.process, capmort14.ddl, model.parameters=list(S=list(formula =~ Days13 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days14A = mark(capmort.process, capmort15.ddl, model.parameters=list(S=list(formula =~ Days14 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days15A = mark(capmort.process, capmort16.ddl, model.parameters=list(S=list(formula =~ Days15 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days16A = mark(capmort.process, capmort17.ddl, model.parameters=list(S=list(formula =~ Days16 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days17A = mark(capmort.process, capmort18.ddl, model.parameters=list(S=list(formula =~ Days17 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days18A = mark(capmort.process, capmort19.ddl, model.parameters=list(S=list(formula =~ Days18 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days19A = mark(capmort.process, capmort20.ddl, model.parameters=list(S=list(formula =~ Days19 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days20A = mark(capmort.process, capmort21.ddl, model.parameters=list(S=list(formula =~ Days20 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days21A = mark(capmort.process, capmort22.ddl, model.parameters=list(S=list(formula =~ Days21 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days22A = mark(capmort.process, capmort23.ddl, model.parameters=list(S=list(formula =~ Days22 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days23A = mark(capmort.process, capmort24.ddl, model.parameters=list(S=list(formula =~ Days23 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days24A = mark(capmort.process, capmort25.ddl, model.parameters=list(S=list(formula =~ Days24 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days25A = mark(capmort.process, capmort26.ddl, model.parameters=list(S=list(formula =~ Days25 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days26A = mark(capmort.process, capmort27.ddl, model.parameters=list(S=list(formula =~ Days26 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days27A = mark(capmort.process, capmort28.ddl, model.parameters=list(S=list(formula =~ Days27 + Trans.Type + TurkAge + HandTime, link="loglog")))
  S.Days28A = mark(capmort.process, capmort29.ddl, model.parameters=list(S=list(formula =~ Days28 + Trans.Type + TurkAge + HandTime, link="loglog")))
  
  return(collect.models() )
}

EWT.capmort.threshold.results = EWT.capmort.threshold.models()
EWT.capmort.threshold.results

write.csv(EWT.capmort.threshold.results$model.table, file = "EWT.capmort.threshold.results.csv", row.names=FALSE)
