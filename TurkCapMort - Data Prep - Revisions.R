#Load Packages
lapply(c("dplyr", "RMark", "plyr", "janitor", "chron", "ggplot2"), require, character.only = T)

########################################
#### Creating DSR Encounter History ####
########################################
telem.raw <- read.csv("Telemetry_Data - Telemetry.csv")
trap.raw <- read.csv("Trapping - Data.csv")
pathogen.raw <- read.csv("Pathogen_Status.csv")

#Clean up dataframes by removing unnecessary columns
telem1 <- telem.raw %>% dplyr::select(AlumBand, TransFreq, Date, Fate) 
trap <- trap.raw %>% dplyr::select(AlumBand, Date, TransFreq) %>%
  filter(!is.na(TransFreq)) %>% #remove birds without transmitters
  mutate(Fate = "L") #add fate column for first encounter

telem <- rbind(trap, telem1) %>%
  mutate(Date = as.POSIXlt(Date, format="%m/%d/%Y")) %>%
  mutate(JulDay = Date$yday)
for (i in 1:length(telem$Date)){
  if(telem$Date$year[i] >120){
    telem$JulDay[i] <- telem$Date$yday[i] + (((telem$Date$year[i] - 118)*365)+1)
  }else{telem$JulDay[i] <- telem$Date$yday[i] + ((telem$Date$year[i] - 118)*365)
  }
}


## Subset alive data and mortality data
EWT.Alive<- subset(telem, Fate=="L")
EWT.Dead<- subset(telem, Fate=="D")

#Return First day found for each bird
First.Found <- aggregate(x=EWT.Alive$JulDay, by=list(EWT.Alive$AlumBand), FUN="min")
colnames(First.Found)<- c("ID", "FirstFound")

##Returns last alive record for each bird
Last.Alive<- aggregate(x=EWT.Alive$JulDay, by=list(EWT.Alive$AlumBand), FUN="max")
colnames(Last.Alive)<- c("ID", "LastAlive")

##This may come in handy defining the live/dead (success/failure) column
Last.Alive$Status.alive<- 0

##Gives earliest time the bird is recorded dead.  Obviously not all birds
##die but we will dead with that in the next step
First.Dead<- aggregate(x=EWT.Dead$JulDay, by=list(EWT.Dead$AlumBand), FUN="min")
colnames(First.Dead)<- c("ID", "FirstDead")

First.Dead$Status.dead<- 1

##Merge back together first the two live frames
EWT.Merge.Alive<- merge(First.Found, Last.Alive, by="ID")

## And merge the combined live with the dead.  There is probably a way to do a 3-way merge 
## but this works.  The all.x=TRUE keeps the alive data with no corresponding dead data.
EWT.Merge.All<- merge(EWT.Merge.Alive, First.Dead, by="ID", all.x=TRUE)

##Last Checked
EWT.Merge.All$LastChecked<- pmax(EWT.Merge.All$LastAlive, EWT.Merge.All$FirstDead, na.rm=TRUE)

##Fate column, which will be the pmax() of the status.alive and status.dead columns
EWT.Merge.All$Fate<- pmax(EWT.Merge.All$Status.alive, EWT.Merge.All$Status.dead, na.rm=TRUE)

#We want the Fate to be 0 if the bird survived past 45 days (can change number of days.)
#Also want to set everything in terms of first day the bird was found
for(i in 1:length(EWT.Merge.All$ID)){
  #Set last day found alive in terms of first found
  EWT.Merge.All$LastAlive[i] <- EWT.Merge.All$LastAlive[i] - EWT.Merge.All$FirstFound[i] + 1
  #Set last checked in terms of first found
  EWT.Merge.All$LastChecked[i] <- EWT.Merge.All$LastChecked[i] - EWT.Merge.All$FirstFound[i] + 1
  #Set first day found dead in terms of first found
  if(!is.na(EWT.Merge.All$FirstDead[i])){
    EWT.Merge.All$FirstDead[i] <- EWT.Merge.All$FirstDead[i] - EWT.Merge.All$FirstFound[i] + 1
  }
  #Set first found to 1
  EWT.Merge.All$FirstFound[i] <- 1
  #Set Fate to 0 if bird bird survived more than 45 days
  if(EWT.Merge.All$LastAlive[i] >= 30){
    EWT.Merge.All$Fate[i] <- 0
  }else{
    EWT.Merge.All$Fate[i] <- 1
  }
  
  #We need LastAlive to be 1 less than LastChecked for Dead Birds
  #We need LastAlive = LastChecked for Alive Birds
  if(EWT.Merge.All$Fate[i] == 0){
    EWT.Merge.All$LastChecked[i] <- 30
    EWT.Merge.All$LastAlive[i] <- 30
  }
}

EWT.Merge.censored1 <- EWT.Merge.All %>% 
  filter(is.na(FirstDead) & Fate == 1)
cens_turks1 <- unique(EWT.Merge.censored1$ID)

EWT.Merge.censored2 <- EWT.Merge.All %>% 
  filter(LastChecked > 30)
cens_turks2 <- unique(EWT.Merge.censored2$ID)

cens_turks <- unique(c(cens_turks1, cens_turks2))

'%notin%' <- Negate('%in%') #Create a "not in" function

EWT.EH <- EWT.Merge.All %>%
  dplyr::select(Bird.ID = ID, FirstFound, LastPresent = LastAlive, LastChecked, Fate) %>%
  filter(Bird.ID %notin% cens_turks) %>% #Remove birds that were lost before 45 days
  mutate(Freq = 1)

write.csv(EWT.EH, file = "CapMortEncounterHistory_nocovariates.csv", row.names=FALSE)

### Prepare and Merge Covariates to EH
require(chron)

trap.cov <- trap.raw %>%
  filter(!is.na(TransFreq)) %>%
  dplyr::select(Bird.ID = AlumBand, YearDay = Date, Study.Area, Location, Trans.Type, Sex, TurkAge = Age, 
                Weight = Weight..lbs., Cap.Time = Time.Fired, Release.Time = Release.Time, Hematoma, Pat.Tag = Pat.Tag) %>%
  mutate(YearDay = as.POSIXlt(YearDay, format="%m/%d/%Y")) %>%
  mutate(YearDay = YearDay$yday) %>% #Return the julian day of that year
  mutate(YearDay = ifelse(YearDay > 200, 1, YearDay)) %>% #If caught in decemberm change Julday to 1 to avoid issues
  mutate(Cap.Time = as.POSIXct(Cap.Time, format = "%H:%M")) %>%
  mutate(Release.Time = as.POSIXct(Release.Time, format = "%H:%M")) %>%
  mutate(HandTime = as.numeric(Release.Time - Cap.Time)) %>%
  mutate(HandTime = ifelse(is.na(HandTime), mean(na.omit(HandTime)), HandTime)) %>%
  dplyr::select(-Release.Time, -Cap.Time) %>%
  mutate(Weight = as.numeric(Weight)) %>%
  mutate(Pat.Tag = ifelse(Pat.Tag == "", "N", "Y")) %>%
  mutate(Trans.Type = ifelse(Trans.Type == "VHF_Neck", "Neck", ifelse(Trans.Type == "VHF_Neckmod" | Trans.Type == "VHF_Back" | Trans.Type == "GPS_Back", "Back", NA)))
#There is an error for the NAs introduce, no worries there I think

EWT.EH.cov <- merge(EWT.EH, trap.cov, by = "Bird.ID", all.x = T) %>%
  mutate(Hematoma = ifelse(is.na(Hematoma), 0, Hematoma))

## Pathogen Status ##
#Pathogen_Data <- gs_title("Pathogen_Status")
#Pathogen_sheet <- gs_read(ss = Pathogen_Data)
#pathogen.raw <- as.data.frame(Pathogen_sheet)

infect.cov <- pathogen.raw %>%
  dplyr::select(Bird.ID = BirdID, MG = MG_result , LPDV = LPDV_result , REV = REV_result ) %>%
  mutate(MG = as.numeric(MG)) %>%
  mutate(LPDV = as.numeric(LPDV)) %>%
  mutate( REV = as.numeric(REV))

EWT.EH.cov1 <- merge(EWT.EH.cov, infect.cov, by="Bird.ID", all.x = TRUE) #all.x keeps rows with no LPDV results


#NA for infection status needs to be changed to the average of the sampled individuals.
mean(na.omit(infect.cov$LPDV))
mean(na.omit(infect.cov$MG))
mean(na.omit(infect.cov$REV))

for (i in 1:length(EWT.EH.cov1$LPDV)){
  if(is.na(EWT.EH.cov1$LPDV[i])){
    EWT.EH.cov1$LPDV[i] <- mean(na.omit(infect.cov$LPDV))
  }
}

for (i in 1:length(EWT.EH.cov1$MG)){
  if(is.na(EWT.EH.cov1$MG[i])){
    EWT.EH.cov1$MG[i] <- mean(na.omit(infect.cov$MG))
  }
}

for (i in 1:length(EWT.EH.cov1$REV)){
  EWT.EH.cov1$REV[i] <- ifelse(is.na(EWT.EH.cov1$REV[i]), mean(na.omit(infect.cov$REV)), EWT.EH.cov1$REV[i])
}

### Weather

### Precipitation and Temperature at the nearest site
require(rnoaa)
require(lubridate)
options(noaakey = "tEJZffukUjubhCZapAulcVngXiYjeyPD")

#Load study site locations and get xy and dates for each capture
cap.info <- trap.raw %>%
  filter(!is.na(Trans.Type) & Trans.Type != "") %>%
  select(Bird.ID = AlumBand, Location, Date)

site.info <- read.csv("CaptureSites.csv") %>%
  select(Location = Location.Name, latitude = Latitude, longitude = Longitude) %>%
  merge(., cap.info, by = "Location", all = T) %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>%
  filter(!is.na(Bird.ID))

#Find the nearest weather station
site.date <- site.info %>% select(Date, latitude, longitude) %>% distinct() %>% 
  dplyr::mutate(site.id = row_number()) %>%
  mutate(Date2 = Date + 6) %>%
  select(site.id, latitude, longitude, Date.min = Date, Date.max = Date2)

#Precipitation measured in tenths of a mm
stations.prcp <- ghcnd_stations() %>% filter(first_year < 2018, last_year > 2020) %>%
  filter(element %in% c("PRCP"))
prcp.fun <- function(df, station = stations.prcp){
  #Find nearest station to capture site
  met.stat = meteo_nearby_stations(lat_lon_df = data.frame(id = df[1], latitude = df[2], longitude = df[3]),
                                   limit = 10, var = "prcp", year_min = 2018, year_max = 2020, station_data = station)
  met.stat.df = do.call("rbind", met.stat)
  #Pull weather monitor
  met.pull = meteo_pull_monitors(monitors = met.stat.df$id, var = "PRCP",
                      date_min = df[4], date_max = df[5])
  #Check all data is available
  met.check = meteo_coverage(met.pull, obs_start_date = df[4], obs_end_date = df[5])
  stat.pull = met.check$summary$id[min(which(met.check$summary$total_obs == 7))]
  #Keep only nearest full dataset
  met.subset = met.pull %>% filter(id == stat.pull)
  #If none of the nearest checked have full data, rerun
  if(is.na(mean(met.subset$prcp))){
    prcp.fun(df, station = station %>% filter(!id %in% met.stat.df$id))
  }else{
    
  dist = met.stat.df %>% filter(id == stat.pull)

  data.frame(id = met.subset$id[1],
             distance = dist$distance[1],
             Date = met.subset$date[7],
             Prcp.Day = met.subset$prcp[7], #tenths of mm
             Prcp.Week = mean(met.subset$prcp))
  }
}
prcp.list <- apply(site.date, 1, prcp.fun)
prcp.df <- do.call("rbind", prcp.list) %>%
  select(Prcp.Day, Prcp.Week)

#Temperature measured in tenths of degree Celcius
stations.tavg <- ghcnd_stations() %>% filter(first_year < 2018, last_year > 2020) %>%
  filter(element %in% c("TAVG"))
tavg.fun <- function(df, station = stations.tavg){
  #Find nearest station to capture site
  met.stat = meteo_nearby_stations(lat_lon_df = data.frame(id = df[1], latitude = df[2], longitude = df[3]),
                                   limit = 10, var = "tavg", year_min = 2018, year_max = 2020, station_data = station)
  met.stat.df = do.call("rbind", met.stat)
  #Pull weather monitor
  met.pull = meteo_pull_monitors(monitors = met.stat.df$id, var = "TAVG",
                                 date_min = df[4], date_max = df[5])
  #Check all data is available
  met.check = meteo_coverage(met.pull, obs_start_date = df[4], obs_end_date = df[5])
  stat.pull = met.check$summary$id[min(which(met.check$summary$total_obs == 7))]
  #Keep only nearest full dataset
  met.subset = met.pull %>% filter(id == stat.pull)
  #If none of the nearest checked have full data, rerun
  if(is.na(mean(met.subset$tavg))){
    tavg.fun(df, station = station %>% filter(!id %in% met.stat.df$id))
  }else{
    
    dist = met.stat.df %>% filter(id == stat.pull)
    
    data.frame(id = met.subset$id[1],
               distance = dist$distance[1],
               Date = met.subset$date[7],
               tavg.Day = met.subset$tavg[7], #tenths of mm
               tavg.Week = mean(met.subset$tavg))
  }
}
tavg.list <- apply(site.date, 1, tavg.fun)
tavg.df <- do.call("rbind", tavg.list) %>%
  select(Tavg.Day = tavg.Day, Tavg.Week = tavg.Week)

weather.cov <- as.data.frame(cbind(site.date, prcp.df, tavg.df)) %>%
  dplyr::rename(Date = Date.min) %>%
  merge(., site.info, by = c("latitude", "longitude", "Date")) %>%
  select(Bird.ID, Prcp.Day, Prcp.Week, Tavg.Day, Tavg.Week)
  
EWT.EH.cov2 <- merge(EWT.EH.cov1, weather.cov, by = "Bird.ID", all.x = T)


### Group Versus Solo Release
require(chron)

group.cov1 <- trap.raw %>% filter(!is.na(TransFreq)) %>% 
  select(Bird.ID = AlumBand, Location, Date, Release.Time) %>%
  mutate(DateTime = as.POSIXct(paste(Date, Release.Time, sep = " "),
                               format = "%m/%d/%Y %H:%M"),
         Censor = ifelse(is.na(Release.Time) | Release.Time == "", 1, 0))
group.cov <- group.cov1 %>%
  group_by(Location, DateTime) %>%
  dplyr::summarize(Total = n()) %>%
  mutate(Group = ifelse(Total > 1, 1, 0)) %>%
  merge(test, ., by = c("Location", "DateTime")) %>%
  mutate(Group = ifelse(Censor == 1, NA, Group)) %>%
  select(Bird.ID, Group)

EWT.EH.cov3 <- merge(EWT.EH.cov2, group.cov, by = "Bird.ID", all.x = T)


### Blood Collection Method
bld.smpl <- read.csv("CaptMort_SampType.csv") %>%
  dplyr::rename(Bird.ID = ID) %>% mutate(Fill = 1) %>%
  tidyr::pivot_wider(names_from = "SampType", values_from = "Fill") %>%
  merge(., EWT.EH.cov2 %>% select(Bird.ID), by = "Bird.ID", all.y = T) %>%
  select(Bird.ID, MG = MG_result, QB, BC) %>%
  mutate(SampType = ifelse(QB == 1, 1, NA)) %>%
  mutate(SampType = ifelse(is.na(SampType) & is.na(MG) & BC == 1, 2, SampType)) %>%
  mutate(SampType = ifelse(is.na(SampType) & !is.na(MG) & BC == 1, 3, SampType)) %>%
  mutate(SampType = ifelse(is.na(SampType), 0, SampType)) %>%
  select(Bird.ID, SampType)
  
EWT.EH.cov3 <- merge(EWT.EH.cov2, bld.smpl, by = "Bird.ID", all.x = T)


################################################################################################
### Body Condition
require(lubridate)

bodycondition.raw <- trap.raw %>% 
  filter(Recapture != "Y") %>%
  dplyr::select(Bird.ID = AlumBand, BodyMass = Weight..lbs., Tarsus, CapDate = Date, Age) %>%
  mutate(CapDate = as.Date(CapDate, format = "%m/%d/%Y")) %>%
  mutate(CapDate_J = yday(CapDate)) %>%
  mutate(CapDate_J = ifelse(CapDate_J > 200, (CapDate_J - 365), CapDate_J)) %>%
  mutate(BodyMass = as.numeric(BodyMass)) %>%
  mutate(Tarsus = as.numeric(Tarsus)) %>%
  filter(!is.na(BodyMass)) %>%
  filter(!is.na(Tarsus)) %>%
  mutate(CapYear = year(CapDate)) %>%
  filter(!is.na(Bird.ID))

bodycondition.lm <- lm(BodyMass ~ Tarsus + CapDate + Age, data = bodycondition.raw)
bodycondition.raw$BC_Score <- resid(bodycondition.lm)

bodycondition.cov <- bodycondition.raw %>% dplyr::select(Bird.ID, BodyCondition = BC_Score)

EWT.EH.cov4 <- merge(EWT.EH.cov3, bodycondition.cov, by = "Bird.ID")

### ADDING ADDITIONAL INFO ON YEAR
year.cov <- trap.raw %>%
  filter(AlumBand %in% EWT.EH.cov3$Bird.ID) %>%
  dplyr::select(Bird.ID = AlumBand, Date) %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>%
  mutate(Year = year(Date)) %>%
  dplyr::select(-Date)

EWT.EH.cov5 <- merge(EWT.EH.cov4, year.cov, by = "Bird.ID")

EWT.EH.covFinal <- EWT.EH.cov5 %>% #Slightly reduced sample size to include body condition (loss of 6 birds)
  dplyr::rename(BodyCon = BodyCondition) %>%
  mutate(SexTT = as.factor(ifelse(Sex == "F" & Trans.Type == "Back", "F.B",
                                  ifelse(Sex == "F" & Trans.Type == "Neck", "F.N",
                                         ifelse(Sex == "M" & Trans.Type == "Neck", "M.N", "M.B")))))

write.csv(EWT.EH.covFinal, file = "TurkCapMort_EH.csv", row.names=FALSE)
