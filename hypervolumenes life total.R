

library(hypervolume)
library(effsize)


# cargamos los datos
source("scripts/pca.R")
colnames(data)
mipca <- mipca %>% select(species,site,origin,lifeform,PC1,PC2,PC3,sitexspp)


# quitamos especies con menos de 3 observaciones
nobs <- as.data.frame(table(mipca$sitexspp))
colnames(nobs) <- c("species", "Freq")
nobs <- nobs %>%  filter (Freq > 3)
mipca <- mipca[which(mipca$sitexspp %in% nobs$species),]


# genero una tabla de resultados provisional
com_hyp <- matrix(nrow=N_analyses, ncol=3)
colnames(com_hyp) <- c("annual", "herbaceous perennial", "woody")
# table(unique(mipca[,c(1,4)])$lifeform)
Nspp <- 33
# genero una dataframe para nativas y otro para invasoras
ss_annual <- mipca %>% subset(lifeform=="annual")
ss_herbaceous <- mipca %>% subset(lifeform=="herbaceous perennial")
ss_woody <- mipca %>% subset(lifeform=="woody")

N_analyses <- 49

for (i in 1:N_analyses) {
  # annual
  # selecciono los individuos
  select_indiv <- list()
  ddd <- ss_annual$species %>% unique() %>% sample(Nspp, replace=F)
  for (s in ddd){
    pre_ss_db <- ss_annual %>%  filter(species==s)
    misite <- pre_ss_db$site %>% sample(1) # selecciono las especies de un unico sitio
    pre_ss_db <- pre_ss_db %>% subset(site == misite)
    ss_db <- pre_ss_db[sample( c(1:nrow(pre_ss_db)), 3, replace=F),]
    select_indiv[[s]] <- ss_db %>% select(PC1,PC2,PC3)
  }
  base_hyp <- do.call(rbind, select_indiv)
  # hago el hypervolumen
  hyp<-hypervolume(base_hyp, method='box')
  # almaceno los resultados
  com_hyp[i,"annual"] <- hyp@Volume
  
  # herbaceous
  select_indiv <- list()
  ddd <- ss_herbaceous$species %>% unique() %>% sample(Nspp, replace=F)
  for (s in ddd){
    pre_ss_db <- ss_herbaceous %>%  filter(species==s)
    misite <- pre_ss_db$site %>% sample(1) # selecciono las especies de un unico sitio
    pre_ss_db <- pre_ss_db %>% subset(site == misite)
    ss_db <- pre_ss_db[sample( c(1:nrow(pre_ss_db)), 3, replace=F),]
    select_indiv[[s]] <- ss_db %>% select(PC1,PC2,PC3)
  }
  base_hyp <- do.call(rbind, select_indiv)
  # hago el hypervolumen
  hyp<-hypervolume(base_hyp, method='box')
  # almaceno los resultados
  com_hyp[i,"herbaceous perennial"] <- hyp@Volume
  
  # woody
  select_indiv <- list()
  ddd <- ss_woody$species %>% unique() %>% sample(Nspp, replace=F)
  for (s in ddd){
    pre_ss_db <- ss_woody %>%  filter(species==s)
    misite <- pre_ss_db$site %>% sample(1) # selecciono las especies de un unico sitio
    pre_ss_db <- pre_ss_db %>% subset(site == misite)
    ss_db <- pre_ss_db[sample( c(1:nrow(pre_ss_db)), 3, replace=F),]
    select_indiv[[s]] <- ss_db %>% select(PC1,PC2,PC3)
  }
  base_hyp <- do.call(rbind, select_indiv)
  # hago el hypervolumen
  hyp<-hypervolume(base_hyp, method='box')
  # almaceno los resultados
  com_hyp[i,"woody"] <- hyp@Volume
}

# ponemos los warnings de nuevo
options(warn=0) 


mean(com_hyp[,"annual"])
SE.mean(com_hyp[,"annual"])
mean(com_hyp[,"herbaceous perennial"])
SE.mean(com_hyp[,"herbaceous perennial"])
mean(com_hyp[,"woody"])
SE.mean(com_hyp[,"woody"])


life_total <- matrix(nrow=3, ncol=3)
rownames(life_total) <- c("annual-herbaceous perennial","annual-woody","herbaceous perennial-woody")
colnames(life_total) <- c("EfSize","lowCI","UpCI")

life_total["annual-herbaceous perennial","EfSize"] <- cohen.d(com_hyp[,"annual"], com_hyp[,"herbaceous perennial"], hedges.correction=TRUE)$estimate
life_total["annual-herbaceous perennial","lowCI"] <- cohen.d(com_hyp[,"annual"], com_hyp[,"herbaceous perennial"], hedges.correction=TRUE)$conf.int[1]
life_total["annual-herbaceous perennial","UpCI"]  <- cohen.d(com_hyp[,"annual"], com_hyp[,"herbaceous perennial"], hedges.correction=TRUE)$conf.int[2]

life_total["annual-woody","EfSize"] <- cohen.d(com_hyp[,"annual"], com_hyp[,"woody"], hedges.correction=TRUE)$estimate
life_total["annual-woody","lowCI"] <- cohen.d(com_hyp[,"annual"], com_hyp[,"woody"], hedges.correction=TRUE)$conf.int[1]
life_total["annual-woody","UpCI"]  <- cohen.d(com_hyp[,"annual"], com_hyp[,"woody"], hedges.correction=TRUE)$conf.int[2]

life_total["herbaceous perennial-woody","EfSize"] <- cohen.d(com_hyp[,"herbaceous perennial"], com_hyp[,"woody"], hedges.correction=TRUE)$estimate
life_total["herbaceous perennial-woody","lowCI"] <- cohen.d(com_hyp[,"herbaceous perennial"], com_hyp[,"woody"], hedges.correction=TRUE)$conf.int[1]
life_total["herbaceous perennial-woody","UpCI"]  <- cohen.d(com_hyp[,"herbaceous perennial"], com_hyp[,"woody"], hedges.correction=TRUE)$conf.int[2]

life_total <- life_total %>% as.data.frame()
life_total$community <- rownames(life_total)

ggplot(data=life_total) +
  geom_point(aes(y = community, x = EfSize), size=4, shape=19, fill="white") +
  geom_errorbarh(aes(y = community, xmin = lowCI, xmax = UpCI), height=0.3, size=.7, color="black") +
  theme_bw() +
  geom_vline(xintercept = 0) +
  xlab("Hedges' d") + ylab(" ") +
  theme(axis.title.x = element_text(color="black", size=14),
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12))



