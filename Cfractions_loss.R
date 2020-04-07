## Modeling carbon fractions loss

## Load packages and data ####
library(ggpubr)
library(tidyverse); theme_set(theme_bw())
library(nlme)
library(emmeans)
library(grid)

load('data/data_bags.rdata')

## Select litter samples at time 2 only
t2_ave_L<- t1t2 %>% 
  subset(time == 2) %>% 
  select(-bag.id, -time) %>%
  group_by(soil.prov, soil.type.prov, block, inoc.prov, mesh, horizon) %>% 
  summarise_all(funs(mean), na.rm = TRUE) %>%
  subset(horizon == "Litter")

## Soluble content
mod.solub <- lme(gSoluble_loss ~ soil.type.prov*mesh*inoc.prov, random = ~ 1|block,
                 weights = NULL, data = t2_ave_L)
anova(mod.solub)
plot(mod.solub)

## Comparisons
solub.inoc.mesh.prov <- emmeans(mod.solub, pairwise ~ mesh*inoc.prov*soil.type.prov, adjust = "tukey") #only time:inoc.prov significant 
solub.inoc.mesh.prov.mult <- CLD(solub.inoc.mesh.prov, alpha=0.05, Letters=letters, adjust="tukey", reversed = FALSE)
solub.inoc.mesh.prov.mult$.group <- gsub('\\s+', '', solub.inoc.mesh.prov.mult$.group)

solubloss.plot.inoc.mesh.prov <- ggplot(solub.inoc.mesh.prov.mult, aes(x = soil.type.prov, y = emmean, fill = mesh)) + 
  geom_bar(stat = "identity", colour="black", position=position_dodge(width = .7), width = .6) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.7), width = .2, size = .5) +
  geom_text(aes(y = emmean+SE, label = .group, hjust = .6, vjust = -.4), position=position_dodge(width=0.7), show.legend = FALSE) +
  scale_x_discrete(labels= c("AM", "EcM")) +
  facet_grid(~inoc.prov, labeller = as_labeller(c('AM'= 'Residence in AM forest', 'ECM'= 'Residence in EcM forest'))) +
  scale_fill_manual(name="Mycorrhizal\nhyphae\nexcluded", labels=c("Yes", "No"), values = c("#d7191c", "#f7f7f7")) +
  labs(x="Soil provenance", y= "Litter soluble content loss (%)") +
  theme(panel.grid.major.x = element_blank())
# Color strip
g <- ggplot_gtable(ggplot_build(solubloss.plot.inoc.mesh.prov))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("#b2df8a","#33a02c")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g, recording=TRUE)

## Hemicellulose
mod.hemicel <- lme(gHemi_loss ~ soil.type.prov*mesh*inoc.prov, random = ~ 1|block,
                   weights = NULL, data = t2_ave_L)
anova(mod.hemicel)
plot(mod.hemicel)

## Comparisons
hemicel.inoc.mesh.prov <- emmeans(mod.hemicel, pairwise ~ mesh*inoc.prov*soil.type.prov, adjust = "tukey") #only time:inoc.prov significant 
hemicel.inoc.mesh.prov.mult <- CLD(hemicel.inoc.mesh.prov, alpha=0.05, Letters=letters, adjust="tukey", reversed = FALSE)
hemicel.inoc.mesh.prov.mult$.group <- gsub('\\s+', '', hemicel.inoc.mesh.prov.mult$.group)

hemicelloss.plot.inoc.mesh.prov <- ggplot(hemicel.inoc.mesh.prov.mult, aes(x = soil.type.prov, y = emmean, fill = mesh)) + 
  geom_bar(stat = "identity", colour="black", position=position_dodge(width = .7), width = .6) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.7), width = .2, size = .5) +
  geom_text(aes(y = emmean+SE, label = .group, vjust = -.4), position=position_dodge(width=0.7), show.legend = FALSE) +
  scale_x_discrete(labels= c("AM", "EcM")) +
  facet_grid(~inoc.prov, labeller = as_labeller(c('AM'= 'Residence in AM forest', 'ECM'= 'Residence in EcM forest'))) +
  scale_fill_manual(name="Mycorrhizal\nhyphae\nexcluded", labels=c("Yes", "No"), values = c("#d7191c", "#f7f7f7")) +
  labs(x="Soil provenance", y= "Litter hemicellulose loss (%)") +
  theme(panel.grid.major.x = element_blank())
# Color strip
g <- ggplot_gtable(ggplot_build(hemicelloss.plot.inoc.mesh.prov))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("#b2df8a","#33a02c")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g, recording=TRUE)

## Cellulose
mod.cel <- lme(gCellulose_loss ~ soil.type.prov*mesh*inoc.prov, random = ~ 1|block,
               weights = NULL, data = t2_ave_L)
anova(mod.cel)
plot(mod.cel)

## Comparisons
cel.inoc <- emmeans(mod.cel, pairwise ~ inoc.prov, adjust = "tukey")
cel.mesh <- emmeans(mod.cel, pairwise ~ mesh, adjust = "tukey")

cel.inoc.mesh.prov <- emmeans(mod.cel, pairwise ~ mesh*inoc.prov*soil.type.prov, adjust = "tukey") #only time:inoc.prov significant 
cel.inoc.mesh.prov.mult <- CLD(cel.inoc.mesh.prov, alpha=0.05, Letters=letters, adjust="tukey", reversed = FALSE)
cel.inoc.mesh.prov.mult$.group <- gsub('\\s+', '', cel.inoc.mesh.prov.mult$.group)
celloss.plot.inoc.mesh.prov <- ggplot(cel.inoc.mesh.prov.mult, aes(x = soil.type.prov, y = emmean, fill = mesh)) + 
  geom_bar(stat = "identity", colour="black", position=position_dodge(width = .7), width = .6) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.7), width = .2, size = .5) +
  geom_text(aes(y = emmean+SE, label = .group, vjust = -.4), position=position_dodge(width=0.7), show.legend = FALSE) +
  scale_x_discrete(labels= c("AM", "EcM")) +
  facet_grid(~inoc.prov, labeller = as_labeller(c('AM'= 'Residence in AM forest', 'ECM'= 'Residence in EcM forest'))) +
  scale_fill_manual(name="Mycorrhizal\nhyphae\nexcluded", labels=c("Yes", "No"), values = c("#d7191c", "#f7f7f7")) +
  labs(x="Soil provenance", y= "Litter cellulose loss (%)") +
  theme(panel.grid.major.x = element_blank())
# Color strip
g <- ggplot_gtable(ggplot_build(celloss.plot.inoc.mesh.prov))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("#b2df8a","#33a02c")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g, recording=TRUE)

## Lignin
mod.lig <- lme(gLignin_loss ~ soil.type.prov*mesh*inoc.prov, random = ~ 1|block,
               weights = NULL, data = t2_ave_L)
anova(mod.lig)
plot(mod.lig)

## Comparisons
lig.inoc <- emmeans(mod.lig, pairwise ~ inoc.prov, adjust = "tukey")
lig.mesh <- emmeans(mod.lig, pairwise ~ mesh, adjust = "tukey")

lig.inoc.mesh.prov <- emmeans(mod.lig, pairwise ~ mesh*inoc.prov*soil.type.prov, adjust = "tukey") #only time:inoc.prov significant 
lig.inoc.mesh.prov.mult <- CLD(lig.inoc.mesh.prov, alpha=0.05, Letters=letters, adjust="tukey", reversed = FALSE)
lig.inoc.mesh.prov.mult$.group <- gsub('\\s+', '', lig.inoc.mesh.prov.mult$.group)
ligloss.plot.inoc.mesh.prov <- ggplot(lig.inoc.mesh.prov.mult, aes(x = soil.type.prov, y = emmean, fill = mesh)) + 
  geom_bar(stat = "identity", colour="black", position=position_dodge(width = .7), width = .6) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.7), width = .2, size = .5) +
  geom_text(aes(y = emmean+SE, label = .group, vjust = -.4), position=position_dodge(width=0.7), show.legend = FALSE) +
  scale_x_discrete(labels= c("AM", "EcM")) +
  facet_grid(~inoc.prov, labeller = as_labeller(c('AM'= 'Residence in AM forest', 'ECM'= 'Residence in EcM forest'))) +
  scale_fill_manual(name="Mycorrhizal\nhyphae\nexcluded", labels=c("Yes", "No"), values = c("#d7191c", "#f7f7f7")) +
  labs(x="Soil provenance", y= "Litter lignin loss (%)") +
  theme(panel.grid.major.x = element_blank())
# Color strip
g <- ggplot_gtable(ggplot_build(ligloss.plot.inoc.mesh.prov))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("#b2df8a","#33a02c")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g, recording=TRUE)
