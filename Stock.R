## Model soil stocks

## Load packages and data ####
library(tidyverse); theme_set(theme_bw())
library(nlme)
library(emmeans)
library(plyr)

load("data/soil_stock.rda") # Calculated from https://doi.org/10.5281/zenodo.3631861 with bulk density
soil_stock <- droplevels(subset(soil_stock, !myco %in% 'mixed')) # Keep only AM- and EcM- dominated plots

# Add stocks for all horizons per plot
stock_plot <- ddply(soil_stock, .(plot, myco, block), numcolwise(sum))

# divide by 1000 to give kg C m-2 soil to 0.2 m depth
stock_plot$totalC <- with(stock_plot, totalC.vol / 1000)

# calculate C:N, C:Po
stock_plot$CN <- with(stock_plot, totalC.vol / totalN.vol)
stock_plot$CPo <- with(stock_plot, totalC.vol / orgP.vol)
stock_plot$PoPi <- with(stock_plot, orgP.vol / inorgP.vol)
stock_plot$PoP <- with(stock_plot, orgP.vol / totalP.vol)

## Total C ####
stock.totalC <- lme(totalC ~ myco, random = ~ 1 | block, weights = varPower(), data = stock_plot)
anova(stock.totalC)

## Comparisons
stock.totalC.tuckey <- emmeans(stock.totalC, pairwise ~ myco, adjust = "tukey")
stock.totalC.mult <- CLD(stock.totalC.tuckey,alpha=0.05,Letters=letters, adjust="tukey")
stock.totalC.mult$.group <- gsub('\\s+', '', stock.totalC.mult$.group)

stock.totalC.plot <- ggplot(stock.totalC.mult, aes(x = myco, y = emmean, fill = myco)) +
  geom_bar(stat = "identity", colour="black", width =.6) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = .1, size = .5) +
  geom_text(aes(y = emmean+SE, label = .group, vjust = -.4, hjust	= .65)) + 
  scale_x_discrete(name = "Forest type", labels= c("AM", "EcM")) +
  scale_fill_manual(values = c("#b2df8a", "#33a02c")) +
  ylim(c(0,15)) +
  ylab(expression(Soil~C~(kg~C~m^-2~soil))) +
  theme(panel.grid.major.x = element_blank(), legend.position = "none")
stock.totalC.plot

## C:N ####
stock.CN.mod <- lme(CN ~ myco, random = ~ 1 | block, data = stock_plot)
anova(stock.CN.mod)

## Comparisons
tuckey.CN <- emmeans(stock.CN.mod, pairwise ~ myco, adjust = "tukey")
mult.CN <- CLD(tuckey.CN,alpha=0.05,Letters=letters, adjust="tukey")
mult.CN$.group <- gsub('\\s+', '', mult.CN$.group)

stock.CN.myco.plot <- ggplot(mult.CN, aes(x = myco, y = emmean, fill = myco)) +
  geom_bar(stat = "identity", colour="black", width =.6) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = .1, size = .5) +
  geom_text(aes(y = emmean+SE, label = .group, vjust = -.4, hjust	= .65)) + 
  scale_x_discrete(name = "Forest type", labels= c("AM", "EcM")) +
  scale_fill_manual(values = c("#b2df8a", "#33a02c")) +
  ylim(c(0,22)) +
  ylab("Soil C:N (stocks)") +
  theme(panel.grid.major.x = element_blank(), legend.position = "none")
stock.CN.myco.plot

## CPo ####
stock.CPo.mod <- lme(CPo ~ myco, random = ~ 1 | block, data = stock_plot)
anova(stock.CPo.mod)

## Comparisons
tuckey.CPo <- emmeans(stock.CPo.mod, pairwise ~ myco, adjust = "tukey")
mult.CPo <- CLD(tuckey.CPo,alpha=0.05,Letters=letters, adjust="tukey")
mult.CPo$.group <- gsub('\\s+', '', mult.CPo$.group)

stock.CPo.myco.plot <- ggplot(mult.CPo, aes(x = myco, y = emmean, fill = myco)) +
  geom_bar(stat = "identity", colour="black", width =.6) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = .1, size = .5) +
  scale_x_discrete(name = "Forest type", labels= c("AM", "EcM")) +
  scale_fill_manual(values = c("#b2df8a", "#33a02c")) +
  ylim(c(0,750)) +
  ylab("Soil C:Po (stocks)") +
  theme(panel.grid.major.x = element_blank(), legend.position = "none")
stock.CPo.myco.plot

