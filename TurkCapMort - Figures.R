#Load Packages
lapply(c("dplyr", "RMark", "plyr", "janitor", "chron", "ggplot2", "patchwork"), require, character.only = T)

#Run Analysis if necessary
source("TurkCapMort - Analysis.R")

#######################################################################################################################
### PLOT FOR Transmitter * Sex Graph
TTSex.survivals <- EWT.capmort.results$S.TTSex$results$real
TTSex.survivals$DaysPostCap <- rep(1:29)
TTSex.survivals$TTSex <- c(rep("F.Back", 29), rep("F.Neck", 29), rep("M.Back", 29), rep("M.Neck", 29))
TTSex.survivals$TTSex <- as.factor(TTSex.survivals$TTSex)

TTSex.est <- TTSex.survivals %>%
  dplyr::select(TTSex, DaysPostCap, DailyS = estimate, DailyLCL = lcl, DailyUCL = ucl) %>%
  `rownames<-`( NULL ) %>%
  group_by(TTSex) %>%
  mutate(CumulS = cumprod(DailyS))

TTSex.est.FB <- TTSex.survivals %>%
  filter(TTSex == "F.Back") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))
TTSex.est.FN<- TTSex.survivals %>%
  filter(TTSex == "F.Neck") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))
TTSex.est.MB <- TTSex.survivals %>%
  filter(TTSex == "M.Back") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))
TTSex.est.MN <- TTSex.survivals %>%
  filter(TTSex == "M.Neck") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))
TTSex.est <- rbind(TTSex.est.FB, TTSex.est.FN, TTSex.est.MN)

TTSex.plot.daily <- ggplot(data = TTSex.est, aes(x = DaysPostCap, y=estimate, group = TTSex)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = TTSex), alpha = .3) +
  geom_line(aes(linetype = TTSex, color = TTSex), size = 1.3) +
  labs(x = "Days Post Capture", y = "Daily Survival")+
  theme_classic(base_size = 18) +
  theme(legend.key.width=unit(3,"line")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0))) +
  scale_color_manual(name = "Sex*Transmitter",
                     values = c("#FB5904", "#04A6FB", "#5904FB")) +
  scale_fill_manual(name = "Sex*Transmitter",
                    values = c("#FB5904", "#04A6FB", "#5904FB")) +
  scale_linetype_manual(name = "Sex*Transmitter",
                        values = c("solid", "dashed", "longdash")) +
  theme(legend.position = "none")

TTSex.plot.cum <- ggplot(data = TTSex.est, aes(x = DaysPostCap, y=CumProd, group = TTSex)) +
  geom_ribbon(aes(ymin = CumProdLCL, ymax = CumProdUCL, fill = TTSex), alpha = .3) +
  geom_line(aes(linetype = TTSex, color = TTSex), size = 1.3) +
  labs(x = "Days Post Capture", y = "Cumulative Survival")+
  theme_classic(base_size = 18) +
  theme(legend.key.width=unit(3,"line")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0))) +
  scale_color_manual(name = "Sex*Transmitter",
                     values = c("#FB5904", "#04A6FB", "#5904FB")) +
  scale_fill_manual(name = "Sex*Transmitter",
                    values = c("#FB5904", "#04A6FB", "#5904FB")) +
  scale_linetype_manual(name = "Sex*Transmitter",
                        values = c("solid", "dashed", "longdash")) +
  theme(legend.position = "none")

TTSex.legend <- ggpubr::get_legend( TTSex.plot.cum, position = "bottom")

gridlay <- "
AB
CC
"
TTSex.fig <- TTSex.plot.daily + TTSex.plot.cum + TTSex.legend +
  plot_layout(design = gridlay, heights = c(1, .2)) +
  plot_annotation(tag_levels = list(c("A", "B", "")))

ggsave(TTSex.fig, filename = "./Figures/Figure - TTSex Daily and Cumulative Survival.jpeg",
       width = 15.5, height = 8, dpi = 600)

  
################################################################################################################################
#Create Single gridded plot for top models
###Trans Type
#Data
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

#Daily Survival
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

#Cumulative Survival
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


###TurkAge
#Data
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

#Daily Survival
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

#Cumulative Survival
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


###Handling Times
#Data
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

#Daily Survival
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

#Cumulative Survival
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

combo.fig <- handtime.plot.daily + handtime.plot.cum + handtime.legend + 
  transtype.plot.daily + transtype.plot.cum + transtype.legend +
  turkage.plot.daily + turkage.plot.cum + turkage.legend +
  plot_layout(design = gridlay, widths = c(1, 1, .2, 1, 1, .2, 1, 1, .2)) +
  plot_annotation(tag_levels = list(c("A", "B", "", "C", "D", "", "E", "F", "")))

ggsave(combo.fig, filename = "./Figures/Figure - Daily and Cumulative Survival.jpeg",
       width = 17.5, height = 15, dpi = 600)

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
REV.fig <- REV.plot.daily + REV.plot.cum + REV.legend +
  plot_layout(design = gridlay, heights = c(1, .2)) +
  plot_annotation(tag_levels = list(c("A", "B", "")))

ggsave(REV.fig, filename = "Figure - REV Daily and Cumulative Survival.jpeg",
       width = 15.5, height = 8, dpi = 600)

###Sex
year.survivals <- EWT.capmort.results$S.sex$results$real
year.survivals$DaysPostCap <- rep(1:29)
year.survivals$Year <- c(rep("2018", 29), rep("2019", 29), rep("2020", 29))
year.survivals$Year <- as.factor(year.survivals$Year)

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
REV.fig <- REV.plot.daily + REV.plot.cum + REV.legend +
  plot_layout(design = gridlay, heights = c(1, .2)) +
  plot_annotation(tag_levels = list(c("A", "B", "")))

ggsave(REV.fig, filename = "Figure - REV Daily and Cumulative Survival.jpeg",
       width = 15.5, height = 8, dpi = 600)

#######################################################################################################################
### Ad Hoc Year
year.survivals <- S.Year$results$real
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

gridlay <- "
AB
CC
"
year.fig <- year.plot.daily + year.plot.cum + year.legend +
  plot_layout(design = gridlay, heights = c(1, .2)) +
  plot_annotation(tag_levels = list(c("A", "B", "")))

ggsave(year.fig, filename = "./Figures/Figure - Year Daily and Cumulative Survival.jpeg",
       width = 15.5, height = 8, dpi = 600)

#####################################################################################################################################
### A Priori versus Ad Hoc Coefficient Comparisons

#Compare Beta Coefficients
TT.ah <- as.data.frame(S.TT.solo$results$beta[2,]) %>%
  mutate(ID = "ad hoc") %>%
  mutate(Model = "TransType") %>%
  mutate(Beta = "Descriptor")
TT.ap <- as.data.frame(EWT.capmort.results$S.transtype$results$beta[2,]) %>%
  mutate(ID = "a priori") %>%
  mutate(Model = "TransType") %>%
  mutate(Beta = "Descriptor")
Sex.ah <- as.data.frame(S.Sex.solo$results$beta[2,]) %>%
  mutate(ID = "ad hoc") %>%
  mutate(Model = "Sex") %>%
  mutate(Beta = "Descriptor")
Sex.ap <- as.data.frame(EWT.capmort.results$S.sex$results$beta[2,]) %>%
  mutate(ID = "a priori") %>%
  mutate(Model = "Sex") %>%
  mutate(Beta = "Descriptor")

TT.ah.int <- as.data.frame(S.TT.solo$results$beta[4,]) %>%
  mutate(ID = "ad hoc") %>%
  mutate(Model = "TransType") %>%
  mutate(Beta = "Interaction")
TT.ap.int <- as.data.frame(EWT.capmort.results$S.transtype$results$beta[4,]) %>%
  mutate(ID = "a priori") %>%
  mutate(Model = "TransType") %>%
  mutate(Beta = "Interaction")
Sex.ah.int <- as.data.frame(S.Sex.solo$results$beta[4,]) %>%
  mutate(ID = "ad hoc") %>%
  mutate(Model = "Sex") %>%
  mutate(Beta = "Interaction")
Sex.ap.int <- as.data.frame(EWT.capmort.results$S.sex$results$beta[4,]) %>%
  mutate(ID = "a priori") %>%
  mutate(Model = "Sex") %>%
  mutate(Beta = "Interaction")

AH.AP.comp <- rbind(TT.ah,TT.ap,Sex.ah,Sex.ap,
                    TT.ah.int,TT.ap.int,Sex.ah.int,Sex.ap.int) %>%
  mutate(ID = factor(ID, levels = c("ad hoc", "a priori")),
         Model = factor(Model, levels = c("TransType", "Sex")),
         ShapeColor = factor(paste(ID, Model, sep = "."), 
                             levels = c("a priori.Sex","ad hoc.Sex",
                                        "a priori.TransType", "ad hoc.TransType"),
                             labels = c("Sex (a priori)", "Sex (ad hoc)",
                                        "Transmitter (a priori)", "Transmitter (ad hoc)")))


betacomp <- ggplot(data = AH.AP.comp, aes(y = Beta, x = estimate, shape = ShapeColor, color = ShapeColor)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = 2, size = 1.5) +
  geom_point(size = 8,
             position = position_dodge(width = .4)) +
  geom_errorbar(aes(xmin = lcl, xmax = ucl),
                width = 0, size = 2,
                position = position_dodge(width = .4)) +
  theme_classic(base_size = 18) + 
  xlab("Coefficient Estimate") +
  ylab("") +
  labs(color = "Model Type") + 
  theme(legend.title.align=0.5) + 
  scale_colour_manual(name = "Model",
                      values = c("#7B3DC2", "#C23D41", "#84C23D", "#3DC2BE")) +
  scale_shape_manual(name = "Model",
                     values = c(15, 19, 15, 19)) +
  theme(legend.position = c(.85,.8))


###Real Survival Estimates Plotted
#Ad hoc Transmitter Type
AH.transtype.survivals <- S.TT.solo$results$real
AH.transtype.survivals$DaysPostCap <- rep(1:29)
AH.transtype.survivals$TransType <- c(rep("Backpack", 29), rep("Necklace", 29))
AH.transtype.survivals$TransType <- as.factor(AH.transtype.survivals$TransType)

AH.transtype.est <- AH.transtype.survivals %>%
  dplyr::select(TransType, DaysPostCap, DailyS = estimate, DailyLCL = lcl, DailyUCL = ucl) %>%
  `rownames<-`( NULL ) %>%
  group_by(TransType) %>%
  mutate(CumulS = cumprod(DailyS))

AH.transtype.est.back <- AH.transtype.survivals %>%
  filter(TransType == "Backpack") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))
AH.transtype.est.neck <- AH.transtype.survivals %>%
  filter(TransType == "Necklace") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))

AH.transtype.est <- rbind(AH.transtype.est.back, AH.transtype.est.neck)

AH.transtype.plot.daily <- ggplot(data = AH.transtype.est, aes(x = DaysPostCap, y=estimate, group = TransType)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = TransType), alpha = .3) +
  geom_line(aes(linetype = TransType, color = TransType), size = 1.3) +
  labs(x = "Days Post Release", y = "Daily Survival")+
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
  theme(axis.title = element_blank(),
        legend.position = c(.8,.4))

#A Priori Trans Type
AP.transtype.survivals <- EWT.capmort.results$S.transtype$results$real
AP.transtype.survivals$DaysPostCap <- rep(1:29)
AP.transtype.survivals$TransType <- c(rep("Backpack", 29), rep("Necklace", 29))
AP.transtype.survivals$TransType <- as.factor(AP.transtype.survivals$TransType)

AP.transtype.est <- AP.transtype.survivals %>%
  dplyr::select(TransType, DaysPostCap, DailyS = estimate, DailyLCL = lcl, DailyUCL = ucl) %>%
  `rownames<-`( NULL ) %>%
  group_by(TransType) %>%
  mutate(CumulS = cumprod(DailyS))

AP.transtype.est.back <- AP.transtype.survivals %>%
  filter(TransType == "Backpack") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))
AP.transtype.est.neck <- AP.transtype.survivals %>%
  filter(TransType == "Necklace") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))

AP.transtype.est <- rbind(AP.transtype.est.back, AP.transtype.est.neck)

AP.transtype.plot.daily <- ggplot(data = AP.transtype.est, aes(x = DaysPostCap, y=estimate, group = TransType)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = TransType), alpha = .3) +
  geom_line(aes(linetype = TransType, color = TransType), size = 1.3) +
  labs(x = "Days Post Release", y = "Daily Survival")+
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
  theme(legend.position = c(.8,.4))


#Ad hoc sex
AH.sex.survivals <- S.Sex.solo$results$real
AH.sex.survivals$DaysPostCap <- rep(1:29)
AH.sex.survivals$Sex <- c(rep("Female", 29), rep("Male", 29))
AH.sex.survivals$Sex <- as.factor(AH.sex.survivals$Sex)

AH.sex.est <- AH.sex.survivals %>%
  dplyr::select(Sex, DaysPostCap, DailyS = estimate, DailyLCL = lcl, DailyUCL = ucl) %>%
  `rownames<-`( NULL ) %>%
  group_by(Sex) %>%
  mutate(CumulS = cumprod(DailyS))

AH.sex.est.F <- AH.sex.survivals %>%
  filter(Sex == "Female") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))
AH.sex.est.M <- AH.sex.survivals %>%
  filter(Sex == "Male") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))

AH.sex.est <- rbind(AH.sex.est.F, AH.sex.est.M)

AH.sex.plot.daily <- ggplot(data = AH.sex.est, aes(x = DaysPostCap, y=estimate, group = Sex)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = Sex), alpha = .3) +
  geom_line(aes(linetype = Sex, color = Sex), size = 1.3) +
  labs(x = "", y = "", title = "Ad Hoc")+
  theme_classic(base_size = 18) +
  theme(legend.key.width=unit(3,"line")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0))) +
  scale_color_manual(name = "Sex",
                     values = c("#FB5904", "#04A6FB")) +
  scale_fill_manual(name = "Sex",
                    values = c("#FB5904", "#04A6FB")) +
  scale_linetype_manual(name = "Sex",
                        values = c("solid", "dashed")) +
  theme(legend.position = c(.8,.4),
        axis.title = element_blank(),
        plot.title = element_text(size = 22, face = "bold", hjust = .5))


#A priori Sex
AP.sex.survivals <- EWT.capmort.results$S.sex$results$real
AP.sex.survivals$DaysPostCap <- rep(1:29)
AP.sex.survivals$Sex <- c(rep("Female", 29), rep("Male", 29))
AP.sex.survivals$Sex <- as.factor(AP.sex.survivals$Sex)


AP.sex.est <- AP.sex.survivals %>%
  dplyr::select(Sex, DaysPostCap, DailyS = estimate, DailyLCL = lcl, DailyUCL = ucl) %>%
  `rownames<-`( NULL ) %>%
  group_by(Sex) %>%
  mutate(CumulS = cumprod(DailyS))

AP.sex.est.F <- AP.sex.survivals %>%
  filter(Sex == "Female") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))
AP.sex.est.M <- AP.sex.survivals %>%
  filter(Sex == "Male") %>%
  mutate(CumProd = ifelse(DaysPostCap == 1, estimate, estimate*dplyr::lag(cumprod(estimate))),
         CumProdLCL = ifelse(DaysPostCap == 1, lcl, lcl*dplyr::lag(cumprod(lcl))),
         CumProdUCL = ifelse(DaysPostCap == 1, ucl, ucl*dplyr::lag(cumprod(ucl))))

AP.sex.est <- rbind(AP.sex.est.F, AP.sex.est.M)

AP.sex.plot.daily <- ggplot(data = AP.sex.est, aes(x = DaysPostCap, y=estimate, group = Sex)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = Sex), alpha = .3) +
  geom_line(aes(linetype = Sex, color = Sex), size = 1.3) +
  labs(x = "", y = "Daily Survival", title = "A Priori")+
  theme_classic(base_size = 18) +
  theme(legend.key.width=unit(3,"line")) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0))) +
  scale_color_manual(name = "Sex",
                     values = c("#FB5904", "#04A6FB")) +
  scale_fill_manual(name = "Sex",
                    values = c("#FB5904", "#04A6FB")) +
  scale_linetype_manual(name = "Sex",
                        values = c("solid", "dashed")) +
  theme(legend.position = c(.8,.4),
        plot.title = element_text(size = 22, face = "bold", hjust = .5))

TT.legend <- ggpubr::get_legend(AH.transtype.plot.daily, position = "right")

Sex.legend <- ggpubr::get_legend(AP.sex.plot.daily, position = "right")

gridlay <- "
AA
BC
DE
"

AP.AH.comp.fig <- betacomp + AP.sex.plot.daily + AH.sex.plot.daily +
  AP.transtype.plot.daily + AH.transtype.plot.daily +
  plot_layout(design = gridlay) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag.position = c(.03,.95))

ggsave(AP.AH.comp.fig, filename = "./Figures/Figure - A Priori Ad Hoc Comp.jpeg",
       width = 14.5, height = 18, dpi = 600)
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
