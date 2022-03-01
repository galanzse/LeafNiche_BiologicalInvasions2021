

source("scripts/hypervolumenes life total.R")
source("scripts/hypervolume pc x life.R")


life_total
eff_life_pc
eff_life_pc$group <- paste(eff_life_pc$group, eff_life_pc$PC)
eff_life_pc$PC <- NULL

eff_life <- rbind(life_total,eff_life_pc)

eff_life <- eff_life %>% mutate(group = factor(group, levels =
                                                       c("herbaceous perennial - woody PC3",
                                                         "annual - woody PC3",
                                                         "annual - herbaceous perennial PC3", 
                                                         "herbaceous perennial - woody PC2", 
                                                         "annual - woody PC2",
                                                         "annual - herbaceous perennial PC2",
                                                         "herbaceous perennial - woody PC1",  
                                                         "annual - woody PC1",
                                                         "annual - herbaceous perennial PC1",
                                                         "herbaceous perennial-woody",
                                                         "annual-woody",
                                                         "annual-herbaceous perennial")))


ggplot(data=eff_life) +
  geom_point(aes(y = group, x = EfSize), size=4, shape=19, fill="white") +
  geom_errorbarh(aes(y = group, xmin = lowCI, xmax = UpCI), height=0.3, size=.7, color="black") +
  theme_bw() +
  geom_vline(xintercept = 0) +
  xlab("Hedges' d") + ylab(" ") +
  theme(axis.title.x = element_text(color="black", size=14),
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12))





