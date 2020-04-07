## Model ratios change

## Load packages and data ####
library(ggpubr)
library(tidyverse); theme_set(theme_bw())
library(nlme)
library(emmeans)
library(grid)

load('data/data_bags.rdata')

## Average by plots
t1t2_ave<- t1t2 %>% 
  select(-bag.id) %>% 
  group_by(soil.prov, soil.type.prov, block, inoc.prov, mesh, horizon, time) %>% summarise_all(funs(mean), na.rm = TRUE)

## C:N ratio ####
t1t2_ave <- t1t2_ave %>% 
  mutate(CN = gC_t0/gN_t0, CN_t = gC_t/gN_t) %>% 
  mutate(CN_change = CN/CN_t)

mod.CN <- lme(CN_change ~ soil.type.prov*horizon*time*mesh + mesh*inoc.prov*horizon*time, random = ~ 1|block,
              weights =  varIdent(form = ~ 1 | horizon), data = t1t2_ave)
anova(mod.CN)
plot(mod.CN)

## Comparisons
CN.horiz.prov.mesh <- emmeans(mod.CN, pairwise ~ horizon*soil.type.prov*mesh, adjust = "tukey")
CN.mult.horiz.prov.mesh <- CLD(CN.horiz.prov.mesh, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
CN.mult.horiz.prov.mesh$.group <- gsub('\\s+', '', CN.mult.horiz.prov.mesh$.group)
CN.plot.mult.horiz.prov.mesh <- ggplot(CN.mult.horiz.prov.mesh, aes(x = horizon, y = emmean-1, fill = mesh)) + 
  geom_bar(stat = "identity", colour="black", position=position_dodge(width = .7), width = .6) +
  facet_grid(~soil.type.prov, labeller = labeller(soil.type.prov = c('AS'='AM soil provenance', 'FG'='EcM soil provenance'))) +
  geom_errorbar(aes(ymin = emmean-1-SE, ymax = emmean-1+SE), position=position_dodge(width=0.7), width = .2, size = .5) +
  geom_text(aes(y = emmean-1-SE, label = .group, vjust = 1.25), position=position_dodge(width=0.7), show.legend = FALSE) +
  scale_fill_manual(name="Mycorrhizal\nhyphae\nexcluded", labels=c("Yes", "No"), values = c("red", "white")) +
  labs(x="Horizons", y= "C:N ratio change") +
  scale_y_continuous(labels = function(y) y + 1)
# Color strip
g <- ggplot_gtable(ggplot_build(CN.plot.mult.horiz.prov.mesh))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("#b2df8a","#33a02c")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g, recording=TRUE)

## Lignin:N ratio ####
LN_ave <- t1t2_ave %>%  
  mutate(ligN = gLignin_t0/gN_t0, ligN_t = gLignin_t/gN_t) %>% 
  mutate(ligN_change = ligN/ligN_t) %>% 
  subset(horizon == "Litter")

mod.LN <- lme(ligN_change ~ soil.type.prov*mesh*inoc.prov, random = ~ 1|block,
              weights =  NULL, data = LN_ave)
anova(mod.LN)
plot(mod.LN)

## Comparisons
LN.horiz.prov.mesh <- emmeans(mod.LN, pairwise ~ soil.type.prov*mesh, adjust = "tukey")
LN.mult.horiz.prov.mesh <- CLD(LN.horiz.prov.mesh, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
LN.mult.horiz.prov.mesh$.group <- gsub('\\s+', '', LN.mult.horiz.prov.mesh$.group)
LN.plot.mult.horiz.prov.mesh <- ggplot(LN.mult.horiz.prov.mesh, aes(x = soil.type.prov, y = emmean-1, fill = mesh)) + 
  geom_bar(stat = "identity", colour="black", position=position_dodge(width = .7), width = .6) +
  geom_errorbar(aes(ymin = emmean-1-SE, ymax = emmean-1+SE), position=position_dodge(width=0.7), width = .2, size = .5) +
  geom_text(aes(y = emmean-1-SE, label = .group, vjust = 1.25), position=position_dodge(width=0.7), show.legend = FALSE) +
  scale_fill_manual(name="Mycorrhizal\nhyphae\nexcluded", labels=c("Yes", "No"), values = c("red", "white")) +
  labs(x="Soil provenance ", y= "Litter lignin:N ratio change") +
  scale_x_discrete(labels= c("AM", "EcM")) +
  scale_y_continuous(labels = function(y) y + 1, limits = c(-0.1,0.23))
# Color strip
g <- ggplot_gtable(ggplot_build(LN.plot.mult.horiz.prov.mesh))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("#b2df8a","#33a02c")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g, recording=TRUE)

## d15N ####
d15N_ave <- t1t2_ave %>%  
  subset(time == 2) %>% 
  subset(horizon == "Fragmented") %>% 
  mutate(dN = dN15.t0-d15N) %>% 
  group_by(soil.prov, soil.type.prov, block, inoc.prov, mesh, horizon, time) %>% summarise_all(funs(mean), na.rm = TRUE)

d15N.mod <- lme(d15N ~ soil.type.prov * mesh * inoc.prov, random = ~ 1|block, data = d15N_ave,
                weights = varPower())
anova(d15N.mod)
plot(d15N.mod)

## Comparisons
d15N.inoc <- emmeans(d15N.mod, pairwise ~ inoc.prov, adjust = "tukey")

d15N.mod.contrast <- emmeans(d15N.mod, pairwise ~ soil.type.prov*mesh*inoc.prov, adjust = "tukey")
d15N.mod.mult <- CLD(d15N.mod.contrast, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
d15N.mod.mult$.group <- gsub('\\s+', '', d15N.mod.mult$.group)

d15N.mod.mult.plot  <- ggplot(d15N.mod.mult, aes(x = soil.type.prov, y = emmean, fill = mesh)) + 
  geom_bar(stat = "identity", colour="black", position=position_dodge(width = .7), width = .6) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.7), width = .2, size = .5) +
  geom_text(aes(y = emmean+SE, label = .group, hjust = .6, vjust = -.4), position=position_dodge(width=0.7), show.legend = FALSE) +
  scale_x_discrete(labels= c("AM", "EcM")) +
  facet_grid(~inoc.prov, labeller = as_labeller(c('AM'= 'Residence in AM forest', 'ECM'= 'Residence in EcM forest'))) +
  scale_fill_manual(name="Mycorrhizal\nhyphae\nexcluded", labels=c("Yes", "No"), values = c("#d7191c", "#f7f7f7")) +
  labs(x="Soil provenance", y= expression(paste(delta^{15},'N'))) +
  theme(panel.grid.major.x = element_blank())
# Color strip
g <- ggplot_gtable(ggplot_build(d15N.mod.mult.plot))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("#b2df8a","#33a02c")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g, recording=TRUE)

