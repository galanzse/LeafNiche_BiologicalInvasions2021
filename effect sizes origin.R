

source("scripts/hypervolumenes sites.R")
source("scripts/hypervolumene pc x origin.R")


eff_origin
eff_origin$PC <- NULL
hyp_sites


eff_total <- rbind(eff_origin,hyp_sites)
eff_total$group <- rownames(eff_total)
rownames(eff_total) <- NULL
eff_total$group <- as.factor(eff_total$group)

eff_total <- eff_total %>% mutate(group = factor(group, levels =
                                                c('sclerophyll woodland', 'renosterveld',
                                                'acid sands fynbos', 'coastal sage scrub',
                                                'serpentine grassland', 'coastal grassland',
                                                'coastal banksia woodland', 'banksia woodland',
                                                'PC3','PC2','PC1','total')))

ggplot(data=eff_total) +
  geom_point(aes(y = group, x = EfSize), size=4, shape=19, fill="white") +
  geom_errorbarh(aes(y = group, xmin = lowCI, xmax = UpCI), height=0.3, size=.7, color="black") +
  theme_bw() +
  geom_vline(xintercept = 0) +
  xlab("Hedges' d") + ylab(" ") +
  theme(axis.title.x = element_text(color="black", size=14),
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12))







