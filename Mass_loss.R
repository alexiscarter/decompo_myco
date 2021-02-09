## Modeling mass loss

## Load packages and data ####
library(ggpubr)
library(tidyverse); theme_set(theme_bw())
library(nlme)
library(emmeans)

load('data/data_bags.rdata')

## Select bags at time 2 only
t2 <- subset(t1t2, time == 2)

## Average per plot
t1t2_ave<- t1t2 %>% 
  select(-bag.id) %>%
  group_by(soil.prov, soil.type.prov, block, inoc.prov, mesh, horizon, time) %>% 
  summarise_all(funs(mean), na.rm = TRUE) # averaged per replicated bags

t2_ave<- t2 %>% 
  select(-bag.id, -time) %>%
  group_by(soil.prov, soil.type.prov, block, inoc.prov, mesh, horizon) %>% 
  summarise_all(funs(mean), na.rm = TRUE) # averaged per replicated bags

## Full model ####
mod.wt <- lme(wt_loss ~ soil.type.prov*horizon*time + mesh*inoc.prov*horizon*time, random = ~ 1|block, weights = varPower(), data = t1t2_ave)
mod.wt.full.anova <- anova(mod.wt)
plot(mod.wt)

## Pairwise comparisons
## Horizon x Residence x Exclusion x Time
WT.horiz.inoc.mesh <- emmeans(mod.wt, pairwise ~ horizon*inoc.prov*mesh*time, adjust = "tukey")

## Residence x Time
WT.inoc.time <- emmeans(mod.wt, pairwise ~ inoc.prov*time, adjust = "tukey") #only time:inoc.prov significant 
WT.mult.inoc.time <- CLD(WT.inoc.time, alpha=0.05, Letters=letters, adjust="tukey", reversed = FALSE)

## Residence x Time x Exclusion
WT.inoc.mesh.time <- emmeans(mod.wt, pairwise ~ mesh*inoc.prov*time, adjust = "tukey")
WT.mult.inoc.mesh.time <- CLD(WT.inoc.mesh.time, alpha=0.05, Letters=letters, adjust="tukey", reversed = FALSE)
WT.mult.inoc.mesh.time$.group <- gsub('\\s+', '', WT.mult.inoc.mesh.time$.group)
WTloss.plot.inoc.mesh.time <- ggplot(WT.mult.inoc.mesh.time, aes(x = inoc.prov, y = emmean, fill = mesh)) + 
  geom_bar(stat = "identity", colour="black", position=position_dodge(width = .7), width = .6) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.7), width = .2, size = .5) +
  geom_text(aes(y = emmean+SE, label = .group, vjust = -.4), position=position_dodge(width=0.7), show.legend = FALSE) +
  scale_x_discrete(labels= c("AM", "EcM")) +
  ylim(c(0,43)) +
  facet_grid(~time, labeller = as_labeller(c('1'= 'After one year', '2'= 'After two years'))) +
  scale_fill_manual(name="Mycorrhizal\nhyphae\nexcluded", labels=c("Yes", "No"), values = c("#d7191c", "#f7f7f7")) +
  labs(x="Forest of residence", y= "Mass loss (%)") +
  theme(panel.grid.major.x = element_blank())

## Provenance x Horizons x Time
WT.type.horiz.time <- emmeans(mod.wt, pairwise ~ soil.type.prov*horizon*time, adjust = "tukey") #only time:inoc.prov significant 
WT.mult.type.horiz.time <- CLD(WT.type.horiz.time, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
WT.mult.type.horiz.time$.group <- gsub('\\s+', '', WT.mult.type.horiz.time$.group)
WT.plot.mult.type.horiz.time <- ggplot(WT.mult.type.horiz.time, aes(x = horizon, y = emmean, fill = soil.type.prov)) + 
  geom_bar(stat = "identity", colour="black", position=position_dodge(width = .7), width = .6) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.7), width = .2, size = .5) +
  geom_text(aes(y = emmean+SE, label = .group, vjust = -.4), position=position_dodge(width=0.7), show.legend = FALSE) +
  ylim(c(0,50)) +
  facet_grid(~time, labeller = as_labeller(c('1'= 'After one year', '2'= 'After two years'))) +
  scale_fill_manual(name="Soil\nprovenance", labels=c("AM", "EcM"), values = c("#b2df8a", "#33a02c")) +
  labs(x="Horizons", y= "Mass loss (%)") +
  theme(panel.grid.major.x = element_blank())

## HFA figure
WT <- emmeans(mod.wt, pairwise ~ inoc.prov*soil.type.prov*mesh*time, adjust = "tukey") #only time:inoc.prov significant 
WT.mult <- CLD(WT, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
WT.mult$.group <- gsub('\\s+', '', WT.mult$.group)
WT.mult$time <- gsub('1', 'one', WT.mult$time)
WT.mult$time <- gsub('2', 'two', WT.mult$time)

WT.plot.mult.type.horiz.inoc <- ggplot(WT.mult, aes(x = soil.type.prov, y = emmean, fill = inoc.prov)) + 
  geom_bar(stat = "identity", colour="black", position=position_dodge(width = .7), width = .6) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.7), width = .2, size = .5) +
  geom_text(aes(y = emmean+SE, label = .group, vjust = -.4), position=position_dodge(width=0.7), show.legend = FALSE) +
  ylim(c(0,43)) +
  facet_grid(mesh~time, labeller = as_labeller(c('one'= 'After one year', 'two'= 'After two years', '1'= '1', '44'='44')))+#, 'apetit'= expression(paste(delta^{15},'N', ' (‰)')), 'bgrand'= 'mesh 44'))) +
  scale_fill_manual(name="Forest of\nincubation", labels=c("AM", "EcM"), values = c("#b2df8a", "#33a02c")) +
  scale_x_discrete( labels=c("AM", "EcM")) +
  labs(x="Soil provenance", y= "Mass loss (%)") +
  theme(panel.spacing=unit(0, "lines"), strip.text.y = element_text(angle = 0), panel.grid.major.x = element_blank(), text = element_text(family = 'Times'), strip.background = element_rect(fill="white"))

## Model for ECM forest ####
ECM <- subset(t2_ave, inoc.prov == "ECM")
mod.wt.ecm <- lme(wt_loss ~ soil.type.prov*horizon + mesh*horizon, random = ~ 1|block, weights = varPower(), data = ECM)
mod.wt.ecm.anova <- anova(mod.wt.ecm)
plot(mod.wt.ecm)

## Exclusion x Horizon
WT.mesh.horiz.ecm <- emmeans(mod.wt.ecm, pairwise ~ mesh*horizon, adjust = "tukey")
WT.mult.mesh.horiz.ecm <- CLD(WT.mesh.horiz.ecm, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
WT.mult.mesh.horiz.ecm$.group <- gsub('\\s+', '', WT.mult.mesh.horiz.ecm$.group)
WTloss.plot.mesh.horiz.ecm <- ggplot(WT.mult.mesh.horiz.ecm, aes(x = horizon, y = emmean, fill = mesh)) + 
  geom_bar(stat = "identity", colour="black", position=position_dodge(width = .7), width = .6) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.7), width = .2, size = .5) +
  geom_text(aes(y = emmean+SE, label = .group, hjust = .5, vjust = -.4), position=position_dodge(width=0.7), show.legend = FALSE) +
  ylim(c(0,54)) +
  #scale_fill_manual(name="Mycorrhizal\nhyphae\nexcluded", labels=c("Yes", "No"), values = c("#d7191c", "#f7f7f7")) +
  scale_fill_manual(name="Pore size\nmesh", labels=c(expression("1 "~mu~"m"), expression("44 "~mu~"m")), values = c("darkgrey", "#f7f7f7")) +
  labs(x="Horizons", y= "Mass loss in EcM forest (%)")

## Exclusion
WT.mesh.ecm <- emmeans(mod.wt.ecm, pairwise ~ mesh, adjust = "tukey")
WT.mult.mesh.ecm <- CLD(WT.mesh.ecm, alpha=0.05, Letters=letters, adjust="tukey", reversed = FALSE)
WTloss.plot.mesh.ecm <- ggplot(WT.mult.mesh.ecm, aes(x = mesh, y = emmean, fill = mesh)) + 
  geom_bar(stat = "identity", colour="black", position=position_dodge(width = .7), width = .6) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.7), width = .15, size = .5) +
  geom_text(aes(y = emmean+SE, label = .group, hjust = .6, vjust = -.4), position=position_dodge(width=0.7), show.legend = FALSE) +
  ylim(c(0,49)) +
  scale_fill_manual(name="Mycorrhizal\nhyphae\nexcluded", labels=c("Yes", "No"), values = c("#d7191c", "#f7f7f7")) +
  theme(axis.text.x=element_blank()) +
  labs(x="", y= "Mass loss in EcM forest (%)")

## Model for AM forest ####
AM <- subset(t2_ave, inoc.prov == "AM")
mod.wt.am <- lme(wt_loss ~ soil.type.prov*horizon + mesh*horizon, random = ~ 1|block, weights = varPower(), data = AM)
mod.wt.am.anova <- anova(mod.wt.am)
plot(mod.wt.am)

## Exclusion x Horizon
WT.mesh.horiz.am <- emmeans(mod.wt.am, pairwise ~ mesh*horizon, adjust = "tukey")
WT.mult.mesh.horiz.am <- CLD(WT.mesh.horiz.am, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
WT.mult.mesh.horiz.am$.group <- gsub('\\s+', '', WT.mult.mesh.horiz.am$.group)
WTloss.plot.mesh.horiz.am <- ggplot(WT.mult.mesh.horiz.am, aes(x = horizon, y = emmean, fill = mesh)) + 
  geom_bar(stat = "identity", colour="black", position=position_dodge(width = .7), width = .6) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.7), width = .2, size = .5) +
  geom_text(aes(y = emmean+SE, label = .group, hjust = .5, vjust = -.4), position=position_dodge(width=0.7), show.legend = FALSE) +
  ylim(c(0,54)) +
  #scale_fill_manual(name="Mycorrhizal\nhyphae\nexcluded", labels=c("Yes", "No"), values = c("#d7191c", "#f7f7f7")) +
  scale_fill_manual(name="Pore size\nmesh", labels=c(expression("1 "~mu~"m"), expression("44 "~mu~"m")), values = c("darkgrey", "#f7f7f7")) +
  labs(x="Horizons", y= "Mass loss in AM forest (%)")

## Exclusion
WT.mesh.am <- emmeans(mod.wt.am, pairwise ~ mesh, adjust = "tukey")
WT.mult.mesh.am <- CLD(WT.mesh.am, alpha=0.05, Letters=letters, adjust="tukey", reversed = FALSE)
WTloss.plot.mesh.am <- ggplot(WT.mult.mesh.am, aes(x = mesh, y = emmean, fill = mesh)) + 
  geom_bar(stat = "identity", colour="black", position=position_dodge(width = .7), width = .6) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.7), width = .15, size = .5) +
  geom_text(aes(y = emmean+SE, label = .group, hjust = .6, vjust = -.4), position=position_dodge(width=0.7), show.legend = FALSE) +
  ylim(c(0,55)) +
  scale_fill_manual(name="Mycorrhizal\nhyphae\nexcluded", labels=c("Yes", "No"), values = c("#d7191c", "#f7f7f7")) +
  theme(axis.text.x=element_blank(), panel.grid.major.x = element_blank()) +
  labs(x="", y= "Mass loss in AM forest (%)")

# Arrange AM+EcM in one plot
wtloss.amecm <- ggarrange(WTloss.plot.mesh.horiz.am, WTloss.plot.mesh.horiz.ecm,
                          labels = c("a", "b"),
                          align = "h",
                          common.legend = TRUE, legend = "right")
