## Modeling initial soil quality

## Load packages and data ####
library(ggpubr)
library(tidyverse); theme_set(theme_bw()) 
library(nlme)
library(emmeans)

load('data/site_t0.rdata')
site_t0 <- mutate(site_t0, ligninN_t0 = Lignin/total.N.t0)

# Initial chemistry
observed_initial_table <- site_t0 %>%
  group_by(soil.type.prov, horizon) %>%
  summarize(TotalC = mean(total.C.t0), TotalC.sd = sd(total.C.t0),
            TotalN = mean(total.N.t0), TotalN.sd = sd(total.N.t0),
            Soluble = mean(NDF.soluble), Soluble.sd = sd(NDF.soluble),
            Hemicellulose = mean(ADF.hemi), Hemicellulose.sd = sd(ADF.hemi),
            Cellulose = mean(ADL.cellulose), Cellulose.sd = sd(ADL.cellulose),
            Lignin.mean = mean(Lignin), Lignin.sd = sd(Lignin))

## Soil quality ####
## Lignin:N
ligN <- lme(ligninN_t0 ~ soil.type.prov*horizon, random = ~ 1|block, weights = NULL, data = site_t0)
anova(ligN)
plot(ligN)

## Comparisons
ligN.prov.horiz <- emmeans(ligN, pairwise ~ soil.type.prov*horizon, adjust = "tukey")
mult.ligN.prov.horiz <- CLD(ligN.prov.horiz, alpha=0.05, Letters=letters, adjust="tukey", reversed = FALSE)

## C:N
CN <- lme(CN.t0 ~ soil.type.prov*horizon, random = ~ 1|block, weights = NULL, data = site_t0)
anova(CN)
plot(CN)

CN.prov.horiz <- emmeans(CN, pairwise ~ soil.type.prov*horizon, adjust = "tukey")
mult.CN.prov.horiz <- CLD(CN.prov.horiz, alpha=0.05, Letters=letters, adjust="tukey", reversed = FALSE)

## dN15
dN <- lme(dN15.t0 ~ soil.type.prov + horizon, random = ~ 1|block, weights = NULL, data = site_t0)
anova(dN)
plot(dN)

dN.prov.horiz <- emmeans(dN, pairwise ~ soil.type.prov*horizon, adjust = "tukey")
mult.dN.prov.horiz <- CLD(dN.prov.horiz, alpha=0.05, Letters=letters, adjust="tukey", reversed = FALSE)

## Plots ####
theme_depth <- theme_bw() + theme(panel.border = element_blank(),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  axis.line = element_line(colour = "darkgrey"),
                                  text = element_text(size=10))

CN.prov.horiz.plot <-ggplot(mult.CN.prov.horiz, aes(x = horizon, y = emmean, color = soil.type.prov, shape = soil.type.prov)) +
  geom_line(aes(color=soil.type.prov, group = soil.type.prov, linetype = soil.type.prov)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.05) +
  scale_color_manual(values = c("#b2df8a", "#33a02c")) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x="Horizons", y = 'C:N ratio', linetype = "Forest type", color = "Forest type", shape = "Forest type", size = 1) +
  theme_depth

ligN.prov.horiz.plot <-ggplot(mult.ligN.prov.horiz, aes(x = horizon, y = emmean, color = soil.type.prov, shape = soil.type.prov)) +
  geom_line(aes(color=soil.type.prov, group = soil.type.prov, linetype = soil.type.prov)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.05) +
  scale_color_manual(values = c("#b2df8a", "#33a02c")) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x="", y = 'Lignin:N ratio', linetype = "Forest type", color = "Forest type", shape = "Forest type") +
  theme_depth

dN.prov.horiz.plot <-ggplot(mult.dN.prov.horiz, aes(x = horizon, y = emmean, color = soil.type.prov, shape = soil.type.prov)) +
  geom_line(aes(color=soil.type.prov, group = soil.type.prov, linetype = soil.type.prov)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.05) +
  scale_color_manual(values = c("#b2df8a", "#33a02c")) +
  coord_flip() +
  scale_y_continuous(position = "right") +
  labs(x="", y = expression(paste(delta^{15},'N', ' (‰)')), linetype = "Forest type", color = "Forest type", shape = "Forest type", size = 1) +
  theme_depth

# Arrange in one plot
site.plot <- ggarrange(CN.prov.horiz.plot, ligN.prov.horiz.plot, dN.prov.horiz.plot,
                       labels = c("a", "b", "c"), font.label = list(size = 10),
                       nrow = 1, ncol = 3,
                       common.legend = TRUE, legend = "right")

