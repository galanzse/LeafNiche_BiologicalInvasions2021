


library(lmerTest)
library(MuMIn)
library(multcomp)
library(scales)
library(ggpubr)



# sumo a cada PC su minimo para llevarlo a espacio positivo
source("scripts/pca.R")
mipca$PC1 <- mipca$PC1 + abs(min(mipca$PC1))
mipca$PC2 <- mipca$PC2 + abs(min(mipca$PC2))
mipca$PC3 <- mipca$PC3 + abs(min(mipca$PC3))


# modelo general ####
colnames(mipca)

# origin

# PC1
mod <- lmer(PC1 ~ 1 + origin + (1 | site/species), data = mipca, REML = FALSE)
summary(mod)
anova(mod)
r.squaredGLMM(mod)

# PC2
mod2 <- lmer(PC2 ~ 1 + origin + (1 | site/species), data = mipca, REML = FALSE)
summary(mod2)
anova(mod2)
r.squaredGLMM(mod2)

# figura
PC.origin <- matrix(ncol=4, nrow=4) %>% as.data.frame()
colnames(PC.origin) <- c("PC","origin","estimate","SE")
PC.origin$PC <- c("PC1","PC1","PC2","PC2")
PC.origin$origin <- c(rep(c("invasive","native"),2))

PC.origin[1,c("estimate","SE")] <- c(4.8213, 0.3810)
PC.origin[2,c("estimate","SE")] <- c(4.8213-0.9969 , 0.1879)
PC.origin[3,c("estimate","SE")] <- c(2.0748, 0.2000)
PC.origin[4,c("estimate","SE")] <- c(2.0748+0.2805, 0.1461)

pd <- position_dodge(width=0.7)
a <- ggplot(data = PC.origin, aes(x = PC, y = estimate, group=origin, fill=origin)) +
  geom_bar(stat="identity", position=pd, width = 0.7, colour="black") +
  geom_errorbar(aes(ymin=estimate, ymax=estimate+SE), width=.3, alpha=1, position=pd) +
  scale_y_continuous(limits=c(0,6), labels = scales::number_format(accuracy = 0.01)) +
  theme_bw() +
  xlab(NULL) + ylab("estimate") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text=element_text(size=12),
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.title = element_text(size=14)) +
  scale_fill_manual(values = c("#F8766D","#00BFC4"))



# life form

# PC1
mod <- lmer(PC1 ~ 1 + lifeform + (1 | site/species), data = mipca, REML = FALSE)
summary(mod)
anova(mod)
r.squaredGLMM(mod)
summary(glht(model = mod, linfct = mcp(lifeform = "Tukey"))) # post-hoc

# PC2
mod2 <- lmer(PC2 ~ 1 + lifeform + (1 | site/species), data = mipca, REML = FALSE)
summary(mod2)
anova(mod2)
r.squaredGLMM(mod2)
summary(glht(model = mod2, linfct = mcp(lifeform = "Tukey"))) # post-hoc

# figura
PC.life <- matrix(ncol=4, nrow=6) %>% as.data.frame()
colnames(PC.life) <- c("PC","life","estimate","SE")
PC.life$PC <- c("PC1","PC1","PC1","PC2","PC2","PC2")
PC.life$life <- c(rep(c("annual","herbaceous perennial", "woody"),2))

PC.life[1,c("estimate","SE")] <- c(4.8976, 0.3400)
PC.life[2,c("estimate","SE")] <- c(4.8976-0.8742, 0.2018)
PC.life[3,c("estimate","SE")] <- c(4.8976-1.4054, 0.2056)
PC.life[4,c("estimate","SE")] <- c(1.8108, 0.1832)
PC.life[5,c("estimate","SE")] <- c(1.8108+0.3084, 0.1419)
PC.life[6,c("estimate","SE")] <- c(1.8108+1.0978, 0.1448)

pd <- position_dodge(width=0.7)
b <- ggplot(data = PC.life, aes(x = PC, y = estimate, group=life, fill=life)) +
  geom_bar(stat="identity", position=pd, width = 0.7, colour="black") +
  geom_errorbar(aes(ymin=estimate, ymax=estimate+SE), width=.3, alpha=1, position=pd) +
  scale_y_continuous(limits=c(0,6), labels = scales::number_format(accuracy = 0.01)) +
  theme_bw() +
  xlab(NULL) + ylab("estimate") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text=element_text(size=12),
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.title = element_text(size=14)) +
  scale_fill_manual(values = c("thistle4","palegreen2","khaki3"))


# library(ggpubr)
ggarrange(a, b, labels=c("A","B"))






# modelos x sitios ####

seg_PCxsite <- read_excel("resultados.xlsx", sheet = "segregacionxsitios") %>% as.data.frame()

str(mipca)
mipca$site <- as.factor(mipca$site)

for (com in unique(mipca$site)) {
  for (PCx in c("PC1","PC2","PC3")) {
    
    # selecciono mi sitio
    mycom <- mipca %>% subset(site==com)
    mycom <- mycom[,which(colnames(mycom) %in% c("species","origin",PCx))]
    colnames(mycom)[which(colnames(mycom)==PCx)] <- "PC"
    # modelo
    mod <- lmer(PC ~ 1 + origin + (1 | species), data = mycom, REML = FALSE)
    # estimates
    est <- fixef(mod)
    seg_PCxsite$estimate[seg_PCxsite$site == com & seg_PCxsite$origin=="inv" & seg_PCxsite$PC==PCx] <- est["(Intercept)"] # invasora
    seg_PCxsite$estimate[seg_PCxsite$site == com & seg_PCxsite$origin=="nat" & seg_PCxsite$PC==PCx] <- sum(est) # nativa
    # se
    se <- vcov(mod, useScale=T) %>% diag() %>% sqrt()
    seg_PCxsite$SE[seg_PCxsite$site == com & seg_PCxsite$origin=="inv" & seg_PCxsite$PC==PCx] <- se[1] # invasora
    seg_PCxsite$SE[seg_PCxsite$site == com & seg_PCxsite$origin=="nat" & seg_PCxsite$PC==PCx] <- se[2] # nativa
  }
}


# grafico PC sitios ####

# he copiado las estimas de los modelos a excel
str(seg_PCxsite)
seg_PCxsite$site <- as.factor(seg_PCxsite$site)
seg_PCxsite$origin <- as.factor(seg_PCxsite$origin)
levels(seg_PCxsite$origin) <- c('invasive','native')
seg_PCxsite <- mutate(seg_PCxsite,
                      origin=factor(origin, levels=c('native','invasive')))
seg_PCxsite <- mutate(seg_PCxsite,
                      site=factor(site, levels=c('banksia woodland','coastal banksia woodland',
                                                 'coastal grassland','serpentine grassland','coastal sage scrub',            
                                                 'acid sands fynbos','renosterveld','sclerophyll woodland' )))
seg_PCxsite <- seg_PCxsite[order(seg_PCxsite$site),]
rownames(seg_PCxsite) <- NULL


# voy a haver tres graficos
seg_PCxsite_PC1 <- seg_PCxsite %>% subset(PC=="PC1")
seg_PCxsite_PC2 <- seg_PCxsite %>% subset(PC=="PC2")
seg_PCxsite_PC3 <- seg_PCxsite %>% subset(PC=="PC3")


# funcion para separar los nombres de los sitios por filas
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

pd <- position_dodge(width=0.7)

a <- ggplot(data=seg_PCxsite_PC1,
            aes(x=site, y=estimate, shape=origin, color=origin, group=interaction(origin, PC))) +
  geom_point(aes(group=interaction(origin, PC)), position=pd, size=4, stroke=0.2) +
  geom_errorbar(position=pd, width=0.5, aes(group=interaction(origin, PC),
                ymin = estimate-SE, ymax = estimate+SE), alpha=0.5) +
  xlab(" ") + ylab("PC1") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(color="black", size=12, margin = margin(0, unit = "cm")),
        axis.text.y = element_text(color="black", size=12),
        axis.title=element_text(size=14)) +
  scale_color_manual(values = c("#00BFC4","#F8766D")) +
  scale_shape_manual(values=c(17,16))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_x_discrete(breaks=unique(seg_PCxsite$site),
                   labels=addline_format(as.character(unique(seg_PCxsite_PC1$site))))

b <- ggplot(data=seg_PCxsite_PC2,
            aes(x=site, y=estimate, shape=origin, color=origin, group=interaction(origin, PC))) +
  geom_point(aes(group=interaction(origin, PC)), position=pd, size=4, stroke=0.2) +
  geom_errorbar(position=pd, width=0.5, aes(group=interaction(origin, PC),
                ymin = estimate-SE, ymax = estimate+SE), alpha=0.5) +
  xlab(" ") + ylab("PC2") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(color="black", size=12, margin = margin(0, unit = "cm")),
        axis.text.y = element_text(color="black", size=12),
        axis.title=element_text(size=14)) +
  scale_color_manual(values = c("#00BFC4","#F8766D")) +
  scale_shape_manual(values=c(17,16))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_x_discrete(breaks=unique(seg_PCxsite$site),
                   labels=addline_format(as.character(unique(seg_PCxsite_PC1$site))))

c <- ggplot(data=seg_PCxsite_PC3,
            aes(x=site, y=estimate, shape=origin, color=origin, group=interaction(origin, PC))) +
  geom_point(aes(group=interaction(origin, PC)), position=pd, size=4, stroke=0.2) +
  geom_errorbar(position=pd, width=0.5, aes(group=interaction(origin, PC),
                                            ymin = estimate-SE, ymax = estimate+SE), alpha=0.5) +
  xlab(" ") + ylab("PC3") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(color="black", size=12, margin = margin(0, unit = "cm")),
        axis.text.y = element_text(color="black", size=12),
        axis.title=element_text(size=14)) +
  scale_color_manual(values = c("#00BFC4","#F8766D")) +
  scale_shape_manual(values=c(17,16))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_x_discrete(breaks=unique(seg_PCxsite$site),
                   labels=addline_format(as.character(unique(seg_PCxsite_PC1$site))))


ggarrange(a, b, c, labels=c("A","B","C"), nrow=3, ncol=1)


# leyenda
ggplot(aes(x=site, y=estimate, fill=origin, group=interaction(origin, PC)), data=seg_PCxsite[1:4,]) +
  geom_bar(stat = "identity", colour="black", position = pd, width=0.7) +
  geom_errorbar(position=pd, width=0.5, aes(ymin=estimate, ymax=estimate+SE), alpha=1) +
  xlab(" ") + ylab("estimate") +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(),
        legend.text=element_text(size=12),
        axis.text.x = element_text(color="black", size=12, margin = margin(1, unit = "cm")),
        axis.text.y = element_text(color="black", size=12),
        axis.title=element_text(size=14))



mod <- lmer(PC3 ~ 1 + origin + (1 | species), data = subset(mipca, site=="sclerophyll woodland"), REML = FALSE)
summary(mod)



