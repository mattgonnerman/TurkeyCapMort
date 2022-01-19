#Load Packages
lapply(c("dplyr", "RMark", "plyr", "janitor", "chron", "ggplot2"), require, character.only = T)

#Set working directory
setwd("E:/Maine Drive/Analysis/Capture Mortality") #at home

#Run Analysis if necessary
source("./TurkCapMort - Analysis.R")

#######################################################################################################################
### PLOT FOR TOP MODEL = YEAR
year.survivals <- EWT.capmort.results$S.Year$results$real
year.survivals$DaysPostCap <- rep(1:29)
year.survivals$Year <- c(rep("2018", 29), rep("2019", 29), rep("2020", 29))
year.survivals$Year <- as.factor(year.survivals$Year)

year.est <- year.survivals %>%
  dplyr::select(Year, DaysPostCap, DailyS = estimate, DailyLCL = lcl, DailyUCL = ucl) %>%
  `rownames<-`( NULL ) %>%
  group_by(Year) %>%
  mutate(CumulS = cumprod(DailyS))

year.est.18 <- year.survivals %>%
  filter(Year == "2018") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))
year.est.19 <- year.survivals %>%
  filter(Year == "2019") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))
year.est.20 <- year.survivals %>%
  filter(Year == "2020") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))
year.est <- rbind(year.est.18, year.est.19, year.est.20)

year.plot.daily <- ggplot(data = year.est, aes(x = DaysPostCap, y=estimate, group = Year)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = Year), alpha = .3) +
  geom_line(aes(linetype = Year, color = Year), size = 1.3) +
  labs(x = "Days Post Capture", y = "Daily Survival")+
  theme_classic(base_size = 18) +
  theme(legend.key.width=unit(3,"line")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0))) +
  scale_color_manual(name = "Year",
                     values = c("#FB5904", "#04A6FB", "#5904FB")) +
  scale_fill_manual(name = "Year",
                    values = c("#FB5904", "#04A6FB", "#5904FB")) +
  scale_linetype_manual(name = "Year",
                        values = c("solid", "dashed", "longdash")) +
  theme(legend.position = "none")

year.plot.cum <- ggplot(data = year.est, aes(x = DaysPostCap, y=CumProd, group = Year)) +
  geom_ribbon(aes(ymin = CumProdLCL, ymax = CumProdUCL, fill = Year), alpha = .3) +
  geom_line(aes(linetype = Year, color = Year), size = 1.3) +
  labs(x = "Days Post Capture", y = "Cumulative Survival")+
  theme_classic(base_size = 18) +
  theme(legend.key.width=unit(3,"line")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0))) +
  scale_color_manual(name = "Year",
                     values = c("#FB5904", "#04A6FB", "#5904FB")) +
  scale_fill_manual(name = "Year",
                    values = c("#FB5904", "#04A6FB", "#5904FB")) +
  scale_linetype_manual(name = "Year",
                        values = c("solid", "dashed", "longdash")) +
  theme(legend.position = "none")

year.legend <- ggpubr::get_legend(year.plot.cum, position = "bottom")
  
require(patchwork)
gridlay <- "
AB
CC
"
year.fig <- year.plot.daily + year.plot.cum + year.legend +
  plot_layout(design = gridlay, heights = c(1, .2)) +
  plot_annotation(tag_levels = list(c("A", "B", "")))

ggsave(year.fig, filename = "Figure - Year Daily and Cumulative Survival.jpeg",
       width = 15.5, height = 8, dpi = 600)

  
################################################################################################################################
#Create Single gridded plot for top models
###Plot Survival over time for Trans Type###
transtype.survivals <- EWT.capmort.results$S.transtype$results$real
transtype.survivals$DaysPostCap <- rep(1:29)
transtype.survivals$TransType <- c(rep("Backpack", 29), rep("Necklace", 29))
transtype.survivals$TransType <- as.factor(transtype.survivals$TransType)

transtype.est <- transtype.survivals %>%
  dplyr::select(TransType, DaysPostCap, DailyS = estimate, DailyLCL = lcl, DailyUCL = ucl) %>%
  `rownames<-`( NULL ) %>%
  group_by(TransType) %>%
  mutate(CumulS = cumprod(DailyS))

transtype.est.back <- transtype.survivals %>%
  filter(TransType == "Backpack") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))
transtype.est.neck <- transtype.survivals %>%
  filter(TransType == "Necklace") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))

transtype.est <- rbind(transtype.est.back, transtype.est.neck)

transtype.plot.daily <- ggplot(data = transtype.est, aes(x = DaysPostCap, y=estimate, group = TransType)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = TransType), alpha = .3) +
  geom_line(aes(linetype = TransType, color = TransType), size = 1.3) +
  labs(x = "", y = "Daily Survival")+
  theme_classic(base_size = 18) +
  theme(legend.key.width=unit(3,"line")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0))) +
  scale_color_manual(name = "Transmitter\nType",
                     values = c("#FB5904", "#04A6FB")) +
  scale_fill_manual(name = "Transmitter\nType",
                     values = c("#FB5904", "#04A6FB")) +
  scale_linetype_manual(name = "Transmitter\nType",
                     values = c("solid", "dashed")) +
  theme(legend.position = "none")

transtype.plot.cum <- ggplot(data = transtype.est, aes(x = DaysPostCap, y=CumProd, group = TransType)) +
  geom_ribbon(aes(ymin = CumProdLCL, ymax = CumProdUCL, fill = TransType), alpha = .3) +
  geom_line(aes(linetype = TransType, color = TransType), size = 1.3) +
  labs(x = "", y = "Cumulative Survival")+
  theme_classic(base_size = 18) +
  theme(legend.key.width=unit(3,"line")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0))) +
  scale_color_manual(name = "Transmitter\nType",
                     values = c("#FB5904", "#04A6FB")) +
  scale_fill_manual(name = "Transmitter\nType",
                    values = c("#FB5904", "#04A6FB")) +
  scale_linetype_manual(name = "Transmitter\nType",
                        values = c("solid", "dashed")) +
  theme(legend.position = "none")

transtype.legend <- ggpubr::get_legend(transtype.plot.cum, position = "right")


##############################################################################################
###Plot Survival over time for TurkAge###
turkage.survivals <- EWT.capmort.results$S.turkage$results$real
turkage.survivals$DaysPostCap <- rep(1:29)
turkage.survivals$TurkAge <- c(rep("Adult", 29), rep("Juvenile", 29))
turkage.survivals$TurkAge <- as.factor(turkage.survivals$TurkAge)

turkage.est <- turkage.survivals %>%
  dplyr::select(TurkAge, DaysPostCap, DailyS = estimate, DailyLCL = lcl, DailyUCL = ucl) %>%
  `rownames<-`( NULL ) %>%
  group_by(TurkAge) %>%
  mutate(CumulS = cumprod(DailyS))

turkage.est.A <- turkage.survivals %>%
  filter(TurkAge == "Adult") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))
turkage.est.J <- turkage.survivals %>%
  filter(TurkAge == "Juvenile") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))

turkage.est <- rbind(turkage.est.A, turkage.est.J)

turkage.plot.daily <- ggplot(data = turkage.est, aes(x = DaysPostCap, y=estimate, group = TurkAge)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = TurkAge), alpha = .3) +
  geom_line(aes(linetype = TurkAge, color = TurkAge), size = 1.3) +
  labs(x = "", y = "Daily Survival")+
  theme_classic(base_size = 18) +
  theme(legend.key.width=unit(3,"line")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0))) +
  scale_color_manual(name = "Turkey\nAge",
                     values = c("#FB5904", "#04A6FB")) +
  scale_fill_manual(name = "Turkey\nAge",
                    values = c("#FB5904", "#04A6FB")) +
  scale_linetype_manual(name = "Turkey\nAge",
                        values = c("solid", "dashed")) +
  theme(legend.position = "none")

turkage.plot.cum <- ggplot(data = turkage.est, aes(x = DaysPostCap, y=CumProd, group = TurkAge)) +
  geom_ribbon(aes(ymin = CumProdLCL, ymax = CumProdUCL, fill = TurkAge), alpha = .3) +
  geom_line(aes(linetype = TurkAge, color = TurkAge), size = 1.3) +
  labs(x = "", y = "Cumulative Survival")+
  theme_classic(base_size = 18) +
  theme(legend.key.width=unit(3,"line")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0))) +
  scale_color_manual(name = "Turkey\nAge",
                     values = c("#FB5904", "#04A6FB")) +
  scale_fill_manual(name = "Turkey\nAge",
                    values = c("#FB5904", "#04A6FB")) +
  scale_linetype_manual(name = "Turkey\nAge",
                        values = c("solid", "dashed")) +
  theme(legend.position = "none")

turkage.legend <- ggpubr::get_legend(turkage.plot.cum, position = "right")


##############################################################################################
###Plot Survival over time for Handling Times###
###Plot survival for different handling times
handtime.m <- mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ HandTime * LN, link="loglog")))
handtime <- covariate.predictions(handtime.m, 
                                  data = data.frame(HandTime = rep(c(-2,0,2), each = 29)),
                                  indices = c(1:29))
handtime.survivals <- handtime$estimates %>%
  dplyr::select(HandTime = covdata, DaysPostCap = model.index, estimate, lcl, ucl) %>%
  distinct() %>%
  mutate(HandTime = factor(HandTime, levels = c(-2,0,2), labels = c("Low", "Mean", "High")))

handtime.est.L <- handtime.survivals %>%
  filter(HandTime == "Low") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))
handtime.est.M <- handtime.survivals %>%
  filter(HandTime == "Mean") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))
handtime.est.H <- handtime.survivals %>%
  filter(HandTime == "High") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))

handtime.est <- rbind(handtime.est.L, handtime.est.M, handtime.est.H)

handtime.plot.daily <- ggplot(data = handtime.est, aes(x = DaysPostCap, y=estimate, group = HandTime)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = HandTime), alpha = .3) +
  geom_line(aes(linetype = HandTime, color = HandTime), size = 1.3) +
  labs(x = "Days Post Capture", y = "Daily Survival")+
  theme_classic(base_size = 18) +
  theme(legend.key.width=unit(3,"line")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0))) +
  scale_color_manual(name = "Handling\nTime",
                     values = c("#FB5904", "#04A6FB", "#5904FB")) +
  scale_fill_manual(name = "Handling\nTime",
                    values = c("#FB5904", "#04A6FB", "#5904FB")) +
  scale_linetype_manual(name = "Handling\nTime",
                        values = c("solid", "dashed", "longdash")) +
  theme(legend.position = "none")

handtime.plot.cum <- ggplot(data = handtime.est, aes(x = DaysPostCap, y=CumProd, group = HandTime)) +
  geom_ribbon(aes(ymin = CumProdLCL, ymax = CumProdUCL, fill = HandTime), alpha = .3) +
  geom_line(aes(linetype = HandTime, color = HandTime), size = 1.3) +
  labs(x = "Days Post Capture", y = "Cumulative Survival")+
  theme_classic(base_size = 18) +
  theme(legend.key.width=unit(3,"line")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0))) +
  scale_color_manual(name = "Handling\nTime",
                     values = c("#FB5904", "#04A6FB", "#5904FB")) +
  scale_fill_manual(name = "Handling\nTime",
                    values = c("#FB5904", "#04A6FB", "#5904FB")) +
  scale_linetype_manual(name = "Handling\nTime",
                        values = c("solid", "dashed", "longdash")) +
  theme(legend.position = "none")

handtime.legend <- ggpubr::get_legend(handtime.plot.cum, position = "right")

### COMBINED FIGURE
require(patchwork)
gridlay <- "
ABC
DEF
GHI
"

combo.fig <- transtype.plot.daily + transtype.plot.cum + transtype.legend +
  turkage.plot.daily + turkage.plot.cum + turkage.legend +
  handtime.plot.daily + handtime.plot.cum + handtime.legend +
  plot_layout(design = gridlay, widths = c(1, 1, .2, 1, 1, .2, 1, 1, .2)) +
  plot_annotation(tag_levels = list(c("A", "B", "", "C", "D", "", "E", "F", "")))

ggsave(combo.fig, filename = "Figure - Daily and Cumulative Survival.jpeg",
       width = 17.5, height = 15, dpi = 600)


### SUB NULL MODELS
##############################################################################################
#REV
REV.m <- mark(capmort.process, capmort.ddl, model.parameters=list(S=list(formula =~ REV * LN,
                                                                         link="loglog")))
REV <- covariate.predictions(REV.m,
                                  data = data.frame(REV = rep(c(min(EWT.EH.cov2$REV),
                                                                max(EWT.EH.cov2$REV)),
                                                              each = 29)),
                                  indices = c(1:29))
REV.survivals <- REV$estimates %>%
  dplyr::select(REV = covdata, DaysPostCap = model.index, estimate, lcl, ucl) %>%
  distinct() %>%
  mutate(REV = factor(REV, levels = c(min(EWT.EH.cov2$REV),max(EWT.EH.cov2$REV)),
                      labels = c("Negative", "Positive")))

REV.est <- REV.survivals %>%
  dplyr::select(REV, DaysPostCap, DailyS = estimate, DailyLCL = lcl, DailyUCL = ucl) %>%
  `rownames<-`( NULL ) %>%
  group_by(REV) %>%
  mutate(CumulS = cumprod(DailyS))

REV.est.neg <- REV.survivals %>%
  filter(REV == "Negative") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))
REV.est.pos <- REV.survivals %>%
  filter(REV == "Positive") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))
REV.est <- rbind(REV.est.neg, REV.est.pos)

REV.plot.daily <- ggplot(data = REV.est, aes(x = DaysPostCap, y=estimate, group = REV)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = REV), alpha = .3) +
  geom_line(aes(linetype = REV, color = REV), size = 1.3) +
  labs(x = "Days Post Capture", y = "Daily Survival")+
  theme_classic(base_size = 18) +
  theme(legend.key.width=unit(3,"line")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0))) +
  scale_color_manual(name = "REV",
                     values = c("#FB5904", "#04A6FB")) +
  scale_fill_manual(name = "REV",
                    values = c("#FB5904", "#04A6FB")) +
  scale_linetype_manual(name = "REV",
                        values = c("solid", "dashed")) +
  theme(legend.position = "none")

REV.plot.cum <- ggplot(data = REV.est, aes(x = DaysPostCap, y=CumProd, group = REV)) +
  geom_ribbon(aes(ymin = CumProdLCL, ymax = CumProdUCL, fill = REV), alpha = .3) +
  geom_line(aes(linetype = REV, color = REV), size = 1.3) +
  labs(x = "Days Post Capture", y = "Cumulative Survival")+
  theme_classic(base_size = 18) +
  theme(legend.key.width=unit(3,"line")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0))) +
  scale_color_manual(name = "REV",
                     values = c("#FB5904", "#04A6FB")) +
  scale_fill_manual(name = "REV",
                    values = c("#FB5904", "#04A6FB")) +
  scale_linetype_manual(name = "REV",
                        values = c("solid", "dashed")) +
  theme(legend.position = "none")

REV.legend <- ggpubr::get_legend(REV.plot.cum, position = "bottom")

require(patchwork)
gridlay <- "
AB
CC
"
year.fig <- REV.plot.daily + REV.plot.cum + REV.legend +
  plot_layout(design = gridlay, heights = c(1, .2)) +
  plot_annotation(tag_levels = list(c("A", "B", "")))

ggsave(year.fig, filename = "Figure - REV Daily and Cumulative Survival.jpeg",
       width = 15.5, height = 8, dpi = 600)


# 
# 
# #Mean Temp Capture Date
# turkmt.survivals <- EWT.capmort.results$S.mt.CapDate
# turkmt.data <- data.frame(mtCDy = c(-20.885, -10.880, -1.348))
# turkmt_real <- covariate.predictions(model = turkmt.survivals, data = turkmt.data, indices = c(1:29))[[1]] %>%
#   select(covdata, estimate, se) %>%
#   mutate(mtCDy = as.factor(covdata))
# turkmt_real$DaysPostCap <- rep(1:29)
# turkMT.plot <- ggplot(data = turkmt_real %>% filter(DaysPostCap < 30), aes(x = DaysPostCap, y=estimate)) +
#   geom_line(aes(linetype = mtCDy), size = 1.3) +
#   labs(x = "Days Post Capture", y = "Daily Probability of Survival", linetype = "Mean Temp\nCapture Day")+
#   theme_classic(base_size = 18) +
#   theme(legend.key.width=unit(3,"line")) +
#   theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0))) +
#   theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0)))
# turkMT.plot
# ggsave("E:/Maine Drive/Projects/Annual Report 2020/Figures/Cap Mort MeanTemp CD Plot.jpeg", width = 8, height = 6, units = "in")
# 
# 
# turkmtw.survivals <- EWT.capmort.results$S.avgmt.Week
# turkmtw.data <- data.frame(avgmtWk = c(-17.549, -9.869, -3.504))
# turkmtw_real <- covariate.predictions(model = turkmtw.survivals, data = turkmtw.data, indices = c(1:29))[[1]] %>%
#   select(covdata, estimate, se) %>%
#   mutate(avgmtWk = as.factor(covdata))
# turkmtw_real$DaysPostCap <- rep(1:29)
# turkMTW.plot <- ggplot(data = turkmtw_real %>% filter(DaysPostCap < 30), aes(x = DaysPostCap, y=estimate)) +
#   geom_line(aes(linetype = avgmtWk), size = 1.3) +
#   labs(x = "Days Post Capture", y = "Daily Probability of Survival", linetype = "Mean Temp\nWeek")+
#   theme_classic(base_size = 18) +
#   theme(legend.key.width=unit(3,"line")) +
#   theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0))) +
#   theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0)))
# turkMTW.plot
# ggsave("E:/Maine Drive/Projects/Annual Report 2020/Figures/Cap Mort MeanTemp Week Plot.jpeg", width = 8, height = 6, units = "in")






#####################################################################################################################################
#Identify threshold date for capture mortality
#Create Table with Model Deviances and Add Column for Threshold Day
capmort.threshold.table <- data.frame(model.table(EWT.capmort.threshold.results))
for(i in 1:length(capmort.threshold.table$S)){
  capmort.threshold.table$Day[i] <- as.numeric(gsub("[^0-9.-]", "", capmort.threshold.table[i,1]))
}

#Create a plot that shows how model deviance changes as a function of threshold value
deviance.plot <- ggplot(capmort.threshold.table, aes(x=Day, y=Deviance)) +
  geom_line(size = 1.3)+ 
  ylim((1+max(capmort.threshold.table$Deviance)),(min(capmort.threshold.table$Deviance)-1)) +
  xlim(0,30)+
  ylab ("Model Deviance") +
  xlab ("Threshold Day") +
  theme_classic(base_size = 18) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0))) +
  geom_vline(xintercept = 4, linetype = "dashed", size = 1.3, color = "red")
#annotation_custom(b.label)
deviance.plot
ggsave("Cap Mort Deviance Plot.jpeg", width = 6, height = 6, units = "in", dpi = 1200)
