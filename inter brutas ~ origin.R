


# si itv se mantiene igual (no proporcional) en todas las especies de una comunidad
# entonces nuestros resultados serian un artefacto de que las espcies invasoras
# tienen menos variabilidad interespecifica



source("scripts/deBello x origin.R")

inter.pc1 <- deBello.ori %>% subset(group=="interspecific" & pc=="PC1") %>% spread("group", "variance")
mod1 <- lmer(interspecific ~ 1 + origin + (1|site), data=inter.pc1)
summary(mod1)
anova(mod1)

inter.pc2 <- deBello.ori %>% subset(group=="interspecific" & pc=="PC2") %>% spread("group", "variance")
mod2 <- lmer(interspecific ~ 1 + origin + (1|site), data=inter.pc2)
summary(mod2)
anova(mod2)


source("scripts/deBello x lifeform.R")

inter.pc1 <- deBello.life %>% subset(group=="interspecific" & pc=="PC1") %>% spread("group", "variance")
mod1 <- lmer(interspecific ~ 1 + life + (1|site), data=inter.pc1)
summary(mod1)
anova(mod1)

inter.pc2 <- deBello.life %>% subset(group=="interspecific" & pc=="PC2") %>% spread("group", "variance")
mod2 <- lmer(interspecific ~ 1 + life + (1|site), data=inter.pc2)
summary(mod2)
anova(mod2)








