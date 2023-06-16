### Find top model for supported univariates
uni.results <- uni.int.aic %>%
  select(S, k = npar, AICc, DeltaAICc, w = weight, Deviance)
write.csv(uni.results, file = "./Results/Univariate AIC Table.csv", row.names=FALSE)

top95.results <- all.aic %>%
  select(S, k = npar, AICc, DeltaAICc, w = weight, Deviance) %>%
  mutate(CumW =cumsum(w)) %>%
  filter(CumW < .95)

top.for.cov <- all.aic[sapply(cov.names, FUN = function(x){min(which(grepl(x, all.aic$S)))}),] %>%
  mutate(Cov = supported.covs) %>%
  select(Cov, S, AICc)


##########################################################################
### Model Weights by parameter
all.aic <- read.csv("./Results/CapMort - AllComboAIC.csv")

cov.model <- sapply(cov.names, FUN = function(x){
  ifelse(grepl(x, all.aic$S), 1, 0)
})

aic.weights.by.cov <- cbind(cov.model, all.aic) %>%
  mutate_at(cov.names, ~ . * weight) %>%
  select(all_of(cov.names)) %>%
  colSums() %>% sort() %>%
  data.frame(CumWeight = ., Covariate = names(.),)

top.models <- apply(cov.model, 2, FUN = function(x){min(which(x == 1))}) %>% sort()
aic.top <- data.frame(TopModel = top.models,
                      Covariate = names(top.models))

allmodels.results <- merge(aic.weights.by.cov, aic.top, by = "Covariate") %>%
  arrange(desc(CumWeight))
write.csv(allmodels.results, "./Results/AllModels - Top and CumW.csv", row.names=F)

##########################################################################
### Top Model Outputs
S.Top = mark(capmort.process, 
             capmort.ddl, 
             model.parameters=list(S=list(formula =~ HandTime + Backpack + Adult + LN, link="loglog")))

S.Top$results$beta


handtime <- covariate.predictions(S.Top, 
                                  data = data.frame(HandTime = rep(c(min(tcm.data.raw$HandTime),
                                                                     0,
                                                                     max(tcm.data.raw$HandTime)), 
                                                                   each = 29),
                                                    Backpack = max(tcm.data.raw$Backpack),
                                                    Adult = max(tcm.data.raw$Adult)),
                                  indices = c(1:29))
handtime.survivals <- handtime$estimates %>%
  dplyr::select(HandTime, DaysPostCap = model.index, estimate, lcl, ucl) %>%
  distinct() %>%
  mutate(HandTime = factor(HandTime, levels = c(min(tcm.data.raw$HandTime),
                                                0,
                                                max(tcm.data.raw$HandTime)),
                           labels = c("Low", "Mean", "High")))

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