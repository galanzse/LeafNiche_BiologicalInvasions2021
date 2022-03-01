

library(tidyverse)
library(lmerTest)
library(ggplot2)
library(ggpubr)
library(scales)


# cargo los datos
source("scripts/deBello x origin.R")
deBello.ori
source("scripts/deBello x lifeform.R")
deBello.life



# por sitios ###
invxsite <- deBello.ori %>% subset(group %in% c("total","intraspecific")) %>% spread("group", "variance")
invxsite$int2total <- invxsite$intraspecific / invxsite$total


# funcion para separar los nombres de los sitios por filas
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

invxsite <- mutate(invxsite,
                      site=factor(site, levels=c('banksia woodland','coastal banksia woodland',
                                                 'coastal grassland','serpentine grassland','coastal sage scrub',            
                                                 'acid sands fynbos','renosterveld','sclerophyll woodland' )))
invxsite <- invxsite[order(invxsite$site),]


pd <- position_dodge(width=0.7)

ggplot(aes(x=site, y=int2total, fill=origin, group=interaction(origin, pc)), data=invxsite) +
  geom_bar(stat = "identity", colour="black", position = pd, width=0.7) +
  xlab(" ") + ylab("intraspecific / total") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(color="black", size=12, margin = margin(0, unit = "cm")),
        axis.text.y = element_text(color="black", size=12),
        axis.title=element_text(size=14)) +
  # scale_fill_manual(values = c("#F8766D","#00BFC4")) +
  scale_fill_manual(values = c("grey41","white")) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_x_discrete(breaks=unique(invxsite$site),
                   labels=addline_format(c("banksia woodland", "coastal banksia woodland",
                                           "coastal grassland","serpentine grassland",
                                           "coastal sage scrub", "acid sands fynbos",
                                           "renosterveld", "sclerophyll woodland")))



