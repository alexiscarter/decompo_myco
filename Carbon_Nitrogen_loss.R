## Model Carbon and Nitrogen loss

## Load packages and data ####
library(ggpubr)
library(tidyverse); theme_set(theme_bw())
library(nlme)
library(emmeans)
library(grid)

load('data/data_bags.rdata')

## Select bags at time 2 only
t2 <- subset(t1t2, time == 2)

## Average per plot (2 replicates or only 1 if problem with the other replicate)
t2_ave<- t2 %>% 
  select(-bag.id) %>% 
  group_by(soil.prov, soil.type.prov, block, inoc.prov, mesh, horizon, time) %>% summarise_all(funs(mean), na.rm = TRUE)

t1t2_ave<- t1t2 %>% 
  select(-bag.id) %>% 
  group_by(soil.prov, soil.type.prov, block, inoc.prov, mesh, horizon, time) %>% summarise_all(funs(mean), na.rm = TRUE)

## Carbon loss ####
C.mod <- lme(gC_loss ~ soil.type.prov*horizon*time + mesh*inoc.prov*horizon*time, random = ~ 1|block, 
             data = t1t2_ave,
             weights = varPower())
anova(C.mod)
plot(C.mod)

## Pairwise comparisons
C.mesh.type.horiz <- emmeans(C.mod, pairwise ~ mesh*inoc.prov*horizon, adjust = "tukey")
C.mult.mesh.type.horiz <- CLD(C.mesh.type.horiz, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
C.mult.mesh.type.horiz$.group <- gsub('\\s+', '', C.mult.mesh.type.horiz$.group)
## Plot
gCloss.plot.mesh.type.horiz <- ggplot(C.mult.mesh.type.horiz, aes(x = horizon, y = emmean, fill = mesh)) + 
  geom_bar(stat = "identity", colour="black", position=position_dodge(width = .7), width = .6) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.7), width = .2, size = .5) +
  geom_text(aes(y = emmean+SE, label = .group, vjust = -.4), position=position_dodge(width=0.7), show.legend = FALSE) +
  ylim(c(0,50)) +
  facet_grid(~inoc.prov, labeller = as_labeller(c('AM'= 'Residence in AM forest', 'ECM'= 'Residence in EcM forest'))) +
  scale_fill_manual(name="Mycorrhizal\nhyphae\nexcluded", labels=c("Yes", "No"), values = c("#d7191c", "#f7f7f7")) +
  labs(x="Horizons", y= "Carbon loss (%)") +
  theme(panel.grid.major.x = element_blank())

## Nitrogen loss ####
N.mod <- lme(gN_loss ~ soil.type.prov*horizon + mesh*inoc.prov*horizon, random = ~ 1|block, data = t2_ave,
             weights = varExp())
anova(N.mod)
plot(N.mod)

## Comparisons
N.inoc <- emmeans(N.mod, pairwise ~ inoc.prov, adjust = "tukey")
N.mesh <- emmeans(N.mod, pairwise ~ mesh, adjust = "tukey")

N.type.horiz <- emmeans(N.mod, pairwise ~ horizon*mesh*inoc.prov, adjust = "tukey")
N.mult.type.horiz <- CLD(N.type.horiz, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
N.mult.type.horiz$.group <- gsub('\\s+', '', N.mult.type.horiz$.group)

gNloss.plot.type.horiz <- ggplot(N.mult.type.horiz, aes(x = horizon, y = emmean, fill = mesh)) + 
  geom_bar(stat = "identity", colour="black", position=position_dodge(width = .7), width = .6) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.7), width = .2, size = .5) +
  geom_text(aes(y = emmean+SE, label = .group, vjust = -.4), position=position_dodge(width=0.7), show.legend = FALSE) +
  facet_grid(~inoc.prov, labeller = as_labeller(c('AM'= 'Residence in AM forest', 'ECM'= 'Residence in EcM forest'))) +
  scale_fill_manual(name="Mycorrhizal\nhyphae\nexcluded", labels=c("Yes", "No"), values = c("#d7191c", "#f7f7f7")) +
  labs(x="Horizons", y= "N loss (%)") +
  theme(panel.grid.major.x = element_blank())
# Color strip
g <- ggplot_gtable(ggplot_build(gNloss.plot.type.horiz))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("#b2df8a","#33a02c")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g, recording=TRUE)

