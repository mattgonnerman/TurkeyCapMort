#Load Packages
lapply(c("dplyr", "RMark", "plyr", "janitor", "chron", "ggplot2"), require, character.only = T)

#Set working directory
setwd("E:/Maine Drive/Analysis/Capture Mortality") #at home


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

#Weather Data
#install.packages("prism")
require(prism)
require(rgdal)
require(raster)
require(sp)
#Create list of unique trapping events
#Subset the point/polygon for each town trapping event took place in
#For each unique event, need temps/precip for day of and following 7 days
#Create new DF with ID, Date of Capture, and Town
#Will use this to pull weather data from PRISM
weather.cov <- trap.raw %>%
  filter(!is.na(Trans.Type)) %>%
  dplyr::select(Bird.ID = AlumBand, Town, CapDate = Date) %>%
  mutate(CapDate = as.Date(CapDate, format = "%m/%d/%Y")) %>%
  mutate(Event = paste(Town, CapDate, sep = "_"))

#Get the location info for the town
#This can be skipped once we have the capture locations (I think)
town.polygons.all <- readOGR("./Maine_Towns", "Maine_Boundaries_Town_Polygon")

#Maybe want to get info by event to save time running code
events <- unique(weather.cov$Event)

#Load PRISM Data
#http://ropensci.github.io/prism/
options(prism.path = "~/prismtmp")
#This will load all PRISM data from the date of first capture through 7 days after the last capture
get_prism_dailys(type = "tmin", minDate = as.Date("2018-01-01"), maxDate = as.Date("2018-04-15"), keepZip = F)
get_prism_dailys(type = "ppt", minDate = as.Date("2018-01-01"), maxDate = as.Date("2018-04-15"), keepZip = F)
get_prism_dailys(type = "tmin", minDate = as.Date("2018-12-01"), maxDate = as.Date("2019-04-15"), keepZip = F)
get_prism_dailys(type = "ppt", minDate = as.Date("2018-12-01"), maxDate = as.Date("2019-04-15"), keepZip = F)
get_prism_dailys(type = "tmin", minDate = as.Date("2019-12-01"), maxDate = as.Date("2020-04-15"), keepZip = F)
get_prism_dailys(type = "ppt", minDate = as.Date("2019-12-01"), maxDate = as.Date("2020-04-15"), keepZip = F)

#create list of PRISM rasters to access
prismdata.list <- ls_prism_data()

weather.event <- as.data.frame(events) %>%
  dplyr::select(Event = events)
for(i in 1:length(events)){
  #townname <- sub("\\_.*", "", events[i]) #Event Town
  eventdate <- gsub("-", x = sub(".*_", "", events[i]), "")
  eventasdate <- as.Date(sub(".*_", "", events[i]))
  
  #Get a polygon for the town boundaries
  townname <- sub("\\_.*", "", events[i])
  town.polygon <- subset(town.polygons.all, town.polygons.all@data$TOWN == townname)
  town.polygon.sp <- spTransform(town.polygon, CRS("+proj=longlat +datum=WGS84"))
  
  #Get the Average min Temp for day of Capture                       
  prismlist.pos <- which(grepl(gsub("-", x = eventasdate, ""), ls_prism_data()[,1]))
  weather.raster <- raster(ls_prism_data(absPath=T)[prismlist.pos[2],2]) #Gets mintemp info as raster
  weather.raster.crop <- mask(weather.raster, town.polygon.sp) #Crops raster to town boundary
  weather.points <- rasterToPoints(weather.raster.crop) #Converts raster to points
  avgmintemp1 <- mean(weather.points[,3]) #Averages the mintemp for that day across points
  
  #Get the Average Precip for day of Capture                       
  weather.raster <- raster(ls_prism_data(absPath=T)[prismlist.pos[1],2]) #Gets mintemp info as raster
  weather.raster.crop <- mask(weather.raster, town.polygon.sp) #Crops raster to town boundary
  weather.points <- rasterToPoints(weather.raster.crop) #Converts raster to points
  avgprecip1 <- mean(weather.points[,3]) #Averages the mintemp for that day across points
  
  #Get the Average min Temp for 1 day post Capture                       
  prismlist.pos <- which(grepl(gsub("-", x = (eventasdate + 1) , ""), ls_prism_data()[,1]))
  weather.raster <- raster(ls_prism_data(absPath=T)[prismlist.pos[2],2]) #Gets mintemp info as raster
  weather.raster.crop <- mask(weather.raster, town.polygon.sp) #Crops raster to town boundary
  weather.points <- rasterToPoints(weather.raster.crop) #Converts raster to points
  avgmintemp2 <- mean(weather.points[,3]) #Averages the mintemp for that day across points
  
  #Get the Average Precip for 1 day post Capture                       
  weather.raster <- raster(ls_prism_data(absPath=T)[prismlist.pos[1],2]) #Gets mintemp info as raster
  weather.raster.crop <- mask(weather.raster, town.polygon.sp) #Crops raster to town boundary
  weather.points <- rasterToPoints(weather.raster.crop) #Converts raster to points
  avgprecip2 <- mean(weather.points[,3]) #Averages the mintemp for that day across points
  
  #Get the Average min Temp for 2 days post Capture                       
  prismlist.pos <- which(grepl(gsub("-", x = (eventasdate + 2) , ""), ls_prism_data()[,1]))
  weather.raster <- raster(ls_prism_data(absPath=T)[prismlist.pos[2],2]) #Gets mintemp info as raster
  weather.raster.crop <- mask(weather.raster, town.polygon.sp) #Crops raster to town boundary
  weather.points <- rasterToPoints(weather.raster.crop) #Converts raster to points
  avgmintemp3 <- mean(weather.points[,3]) #Averages the mintemp for that day across points
  
  #Get the Average Precip for 2 days post Capture                       
  weather.raster <- raster(ls_prism_data(absPath=T)[prismlist.pos[1],2]) #Gets mintemp info as raster
  weather.raster.crop <- mask(weather.raster, town.polygon.sp) #Crops raster to town boundary
  weather.points <- rasterToPoints(weather.raster.crop) #Converts raster to points
  avgprecip3 <- mean(weather.points[,3]) #Averages the mintemp for that day across points
  
  #Get the Average min Temp for 3 days post Capture                       
  prismlist.pos <- which(grepl(gsub("-", x = (eventasdate + 3) , ""), ls_prism_data()[,1]))
  weather.raster <- raster(ls_prism_data(absPath=T)[prismlist.pos[2],2]) #Gets mintemp info as raster
  weather.raster.crop <- mask(weather.raster, town.polygon.sp) #Crops raster to town boundary
  weather.points <- rasterToPoints(weather.raster.crop) #Converts raster to points
  avgmintemp4 <- mean(weather.points[,3]) #Averages the mintemp for that day across points
  
  #Get the Average Precip for 3 days post Capture                       
  weather.raster <- raster(ls_prism_data(absPath=T)[prismlist.pos[1],2]) #Gets mintemp info as raster
  weather.raster.crop <- mask(weather.raster, town.polygon.sp) #Crops raster to town boundary
  weather.points <- rasterToPoints(weather.raster.crop) #Converts raster to points
  avgprecip4 <- mean(weather.points[,3]) #Averages the mintemp for that day across points
  
  #Get the Average min Temp for 4 days post Capture                       
  prismlist.pos <- which(grepl(gsub("-", x = (eventasdate + 4) , ""), ls_prism_data()[,1]))
  weather.raster <- raster(ls_prism_data(absPath=T)[prismlist.pos[2],2]) #Gets mintemp info as raster
  weather.raster.crop <- mask(weather.raster, town.polygon.sp) #Crops raster to town boundary
  weather.points <- rasterToPoints(weather.raster.crop) #Converts raster to points
  avgmintemp5 <- mean(weather.points[,3]) #Averages the mintemp for that day across points
  
  #Get the Average Precip for 4 days post Capture                       
  weather.raster <- raster(ls_prism_data(absPath=T)[prismlist.pos[1],2]) #Gets mintemp info as raster
  weather.raster.crop <- mask(weather.raster, town.polygon.sp) #Crops raster to town boundary
  weather.points <- rasterToPoints(weather.raster.crop) #Converts raster to points
  avgprecip5 <- mean(weather.points[,3]) #Averages the mintemp for that day across points
  
  #Get the Average min Temp for 5 days post Capture                       
  prismlist.pos <- which(grepl(gsub("-", x = (eventasdate + 5) , ""), ls_prism_data()[,1]))
  weather.raster <- raster(ls_prism_data(absPath=T)[prismlist.pos[2],2]) #Gets mintemp info as raster
  weather.raster.crop <- mask(weather.raster, town.polygon.sp) #Crops raster to town boundary
  weather.points <- rasterToPoints(weather.raster.crop) #Converts raster to points
  avgmintemp6 <- mean(weather.points[,3]) #Averages the mintemp for that day across points
  
  #Get the Average Precip for 5 days post Capture                       
  weather.raster <- raster(ls_prism_data(absPath=T)[prismlist.pos[1],2]) #Gets mintemp info as raster
  weather.raster.crop <- mask(weather.raster, town.polygon.sp) #Crops raster to town boundary
  weather.points <- rasterToPoints(weather.raster.crop) #Converts raster to points
  avgprecip6 <- mean(weather.points[,3]) #Averages the mintemp for that day across points
  
  #Get the Average min Temp for 6 days post Capture                       
  prismlist.pos <- which(grepl(gsub("-", x = (eventasdate + 6) , ""), ls_prism_data()[,1]))
  weather.raster <- raster(ls_prism_data(absPath=T)[prismlist.pos[2],2]) #Gets mintemp info as raster
  weather.raster.crop <- mask(weather.raster, town.polygon.sp) #Crops raster to town boundary
  weather.points <- rasterToPoints(weather.raster.crop) #Converts raster to points
  avgmintemp7 <- mean(weather.points[,3]) #Averages the mintemp for that day across points
  
  #Get the Average Precip for 6 days post Capture                       
  weather.raster <- raster(ls_prism_data(absPath=T)[prismlist.pos[1],2]) #Gets mintemp info as raster
  weather.raster.crop <- mask(weather.raster, town.polygon.sp) #Crops raster to town boundary
  weather.points <- rasterToPoints(weather.raster.crop) #Converts raster to points
  avgprecip7 <- mean(weather.points[,3]) #Averages the mintemp for that day across points
  
  weather.event$mtCDy[i] <- avgmintemp1
  #weather.event$mt3dy[i] <- min(c(avgmintemp1, avgmintemp2, avgmintemp3))
  #weather.event$mtWk[i] <- min(c(avgmintemp1, avgmintemp2, avgmintemp3, avgmintemp4, avgmintemp5, avgmintemp6, avgmintemp7))
  #weather.event$avgmt3dy[i] <- mean(c(avgmintemp1, avgmintemp2, avgmintemp3))
  weather.event$avgmtWk[i] <- mean(c(avgmintemp1, avgmintemp2, avgmintemp3, avgmintemp4, avgmintemp5, avgmintemp6, avgmintemp7))
  weather.event$pcpCDy[i] <- avgprecip1
  #weather.event$tpcp3dy[i] <- sum(c(avgprecip1, avgprecip2, avgprecip3))
  #weather.event$tpcpWk[i] <- sum(c(avgprecip1, avgprecip2, avgprecip3, avgprecip4, avgprecip5, avgprecip6, avgprecip7))
  #weather.event$avgpcp3dy[i] <- mean(c(avgprecip1, avgprecip2, avgprecip3))
  weather.event$avgpcpWk[i] <- mean(c(avgprecip1, avgprecip2, avgprecip3, avgprecip4, avgprecip5, avgprecip6, avgprecip7))
}

weather.cov1 <- merge(weather.cov, weather.event, by = "Event", all.x = T) %>%
  dplyr::select(-Event, -Town, -CapDate)

EWT.EH.cov2 <- merge(EWT.EH.cov1, weather.cov1, by = "Bird.ID", all.x = T)
write.csv(EWT.EH.cov2, file = "CapMortEncounterHistory_withcovariates.csv", row.names=FALSE)

################################################################################################
### Body Condition
require(lubridate)

bodycondition.raw <- trap.raw %>% 
  filter(Recapture != "Y") %>%
  dplyr::select(BirdID = AlumBand, BodyMass = Weight..lbs., Tarsus, CapDate = Date, Age) %>%
  mutate(CapDate = as.Date(CapDate, format = "%m/%d/%Y")) %>%
  mutate(CapDate_J = yday(CapDate)) %>%
  mutate(CapDate_J = ifelse(CapDate_J > 200, (CapDate_J - 365), CapDate_J)) %>%
  mutate(BodyMass = as.numeric(BodyMass)) %>%
  mutate(Tarsus = as.numeric(Tarsus)) %>%
  filter(!is.na(BodyMass)) %>%
  filter(!is.na(Tarsus)) %>%
  mutate(CapYear = year(CapDate)) %>%
  filter(!is.na(BirdID))

bodycondition.lm <- lm(BodyMass ~ Tarsus + CapDate + Age, data = bodycondition.raw)
bodycondition.raw$BC_Score <- resid(bodycondition.lm)

bodycondition.cov <- bodycondition.raw %>% dplyr::select(Bird.ID = BirdID, BodyCondition = BC_Score)

EWT.EH.cov3 <- merge(EWT.EH.cov2, bodycondition.cov, by = "Bird.ID")

### ADDING ADDITIONAL INFO ON YEAR
year.cov <- trap.raw %>%
  filter(AlumBand %in% EWT.EH.cov3$Bird.ID) %>%
  dplyr::select(Bird.ID = AlumBand, Date) %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>%
  mutate(Year = year(Date)) %>%
  dplyr::select(-Date)

EWT.EH.cov4 <- merge(EWT.EH.cov3, year.cov, by = "Bird.ID")

EWT.EH.covFinal <- EWT.EH.cov4 %>% #Slightly reduced sample size to include body condition (loss of 6 birds)
  mutate(Study.Area = as.factor(Study.Area),
         Location = as.factor(Location),
         Sex = as.factor(Sex),
         TurkAge = as.factor(TurkAge),
         Trans.Type = as.factor(Trans.Type),
         Hematoma = as.factor(Hematoma),
         Pat.Tag = as.factor(Pat.Tag),
         Year = as.factor(Year)) %>%
  dplyr::rename(BodyCon = BodyCondition) %>%
  dplyr::mutate(Weight = as.numeric(scale(Weight, center = T, scale = T)),
                HandTime = as.numeric(scale(HandTime, center = T, scale = T)),
                mtCDy = as.numeric(scale(mtCDy, center = T, scale = T)),
                avgmtWk = as.numeric(scale(avgmtWk, center = T, scale = T)),
                pcpCDy = as.numeric(scale(pcpCDy, center = T, scale = T)),
                avgpcpWk = as.numeric(scale(avgpcpWk, center = T, scale = T)),
                BodyCon = as.numeric(scale(BodyCon, center = T, scale = T)),
                YearDay = as.numeric(scale(YearDay, center = T, scale = T))) %>%
  ### Want to better imputate missing values for LPDV, etc.
  mutate(MG = ifelse(MG == max(MG) | MG == min(MG), MG, NA),
         LPDV = ifelse(LPDV == max(LPDV) | LPDV == min(LPDV), LPDV, NA),
         REV = ifelse(REV == max(REV) | REV == min(REV), REV, NA)) %>%
  mutate(MG = as.numeric(scale(MG, center = T, scale = T)),
         LPDV = as.numeric(scale(LPDV, center = T, scale = T)),
         REV = as.numeric(scale(REV, center = T, scale = T))) %>%
  mutate(MG = ifelse(is.na(MG), 0, MG),
         LPDV = ifelse(is.na(LPDV), 0, LPDV),
         REV = ifelse(is.na(REV), 0, REV)) %>%
  mutate(SexTT = as.factor(ifelse(Sex == "F" & Trans.Type == "Back", "F.B",
                                  ifelse(Sex == "F" & Trans.Type == "Neck", "F.N",
                                         ifelse(Sex == "M" & Trans.Type == "Neck", "M.N", "M.B")))))

write.csv(EWT.EH.covFinal, file = "CapMortEncounterHistory_withcovariatesBC.csv", row.names=FALSE)
