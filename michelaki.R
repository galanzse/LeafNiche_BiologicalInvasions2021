

source("scripts/pca.R")

library(lmerTest)

mod1 <- lmer(PC1 ~ 1 + (1|origin/species/site) + (1|site), data = mipca)
mod2 <- lmer(PC2 ~ 1 + (1|origin/species/site) + (1|site), data = mipca)
mod3 <- lmer(PC3 ~ 1 + (1|origin/species/site) + (1|site), data = mipca)


vard <- matrix(nrow=5, ncol=3) %>% as.data.frame()
rownames(vard) <- c('intraspecific', 'interspecific', 'community', 'origin', 'residual')
colnames(vard) <- c('PC1','PC2',"PC3")


var1 <- as.data.frame(summary(mod1)$varcor)
var2 <- as.data.frame(summary(mod2)$varcor)
var3 <- as.data.frame(summary(mod3)$varcor)

vard[,"PC1"] <- var1$vcov/sum(var1$vcov)
vard[,"PC2"] <- var2$vcov/sum(var2$vcov)
vard[,"PC3"] <- var3$vcov/sum(var3$vcov)


vard$factor <- rownames(vard)

vard <- gather(vard, "PC", "variance", 1:3)
vard <- mutate(vard, factor=factor(factor, levels=c("intraspecific","interspecific","origin","community")))
vard <- subset(vard, factor != 'residual')
vard <- mutate(vard,  PC=factor( PC, levels=c("PC3","PC2","PC1")))


# plot
ggplot(data=vard, aes(y=variance, x=PC, fill=factor)) +
  geom_bar(stat="identity", colour="black") +
  
  geom_text(aes(label=round(variance, digits = 2)),
            size = 4,  check_overlap = T,
            position = position_stack(vjust = 1), hjust=1, color="black") +
  
  xlab(NULL) + ylab("proportion of variance") +
  theme_bw() +
  theme(legend.position="bottom", legend.title = element_blank(), legend.text=element_text(size=12),
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.title = element_text(size=14)) +
  coord_flip() +
  scale_fill_manual(values = c("grey20","grey40","grey60", "grey80", "white")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))




