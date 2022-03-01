

library(hypervolume)
library(effsize)

# genero una funcion para calcular el SE
SE.mean <- function(d) { sd(d)/sqrt(length(d)) }


# cargamos los datos
source("scripts/pca.R")
colnames(data)
mipca <- mipca %>% select(species,site,origin,lifeform,PC1,PC2,PC3,sitexspp)

# quitamos especies con menos de 3 observaciones
nobs <- as.data.frame(table(mipca$sitexspp))
colnames(nobs) <- c("species", "Freq")
nobs <- nobs %>%  filter (Freq > 3)
mipca <- mipca[which(mipca$sitexspp %in% nobs$species),]


# genero mis tablas de resultados
# media y SE de la riqueza funcional, y overlap
div_means <- expand.grid(c(unique(mipca$site),"total"), unique(mipca$origin))
colnames(div_means) <- c("site","origin")
div_means$overlap_mean <- NA
div_means$overlap_SE <- NA
div_means$hyp_mean <- NA
div_means$hyp_SE <- NA


# tabla de resultados de hipervolumenes
hyp_sites <- matrix(nrow=length(unique(mipca$site))+1, ncol=3)
rownames(hyp_sites) <- c(unique(mipca$site),"total")
colnames(hyp_sites) <- c("EfSize","lowCI","UpCI")


# establecemos el numero de replicas por especie, numerod de especies por sitio,
# y numero de analisis
N_analyses <- 49
N_replicates <- 4

# es necesario ajustar el numero de especies en cada grupo de nativas e invasoras por sitio
N_species <- matrix(nrow=length(unique(mipca$site)), ncol=2) %>% as.data.frame()
colnames(N_species) <- c("site", "N_spp")
N_species$site <- unique(mipca$site)
table(unique(mipca[,1:3])$origin, unique(mipca[,1:3])$site)
N_species$N_spp <- c(9,5,8,5,4,3,9,4)


# eliminamos warnings (Quique me dijo que no hay problema con que
# Log number of observations is less than or equal to the number of dimensions, siempre que eset cerca)
options(warn=-1)

# por comunidades ####
for (com in unique(mipca$site)) {
  
  # establezco del numero de especies por cada grupo de nativas e invasoras
  N_spp <- N_species$N_spp[N_species$site==com]
  
  # genero una dataframe para nativas y otro para invasoras
  ss_native <- mipca %>% subset(site==com & origin=="native")
  ss_invasive <- mipca %>% subset(site==com & origin=="invasive")
  
  # genero una tabla de resultados provisional
  com_hyp <- matrix(nrow=N_analyses, ncol=3)
  colnames(com_hyp) <- c("natives", "invasives", "overlap")
  
  
    for (i in 1:N_analyses) {
    
    # NATIVAS
      # selecciono los individuos
      select_indiv <- list()
      ddd <- ss_native$species %>% unique() %>% sample(N_spp, replace=F)
        for (s in ddd){
          pre_ss_db <- ss_native %>%  filter(species==s)
          ss_db <- pre_ss_db[sample( c(1:nrow(pre_ss_db)), N_replicates, replace=F),]
          select_indiv[[s]] <- ss_db %>% select(PC1,PC2,PC3)
        }
      base_hyp <- do.call(rbind, select_indiv)
      # hago el hypervolumen
      hyp_nat<-hypervolume(base_hyp, method='box')
      # almaceno los resultados
      com_hyp[i,"natives"] <- hyp_nat@Volume
    
    # INVASORAS
      select_indiv <- list()
      ddd <- ss_invasive$species %>% unique() %>% sample(N_spp, replace=F)
        for (s in ddd){
          pre_ss_db <- ss_invasive %>%  filter(species==s)
          ss_db <- pre_ss_db[sample( c(1:nrow(pre_ss_db)), N_replicates, replace=F),]
          select_indiv[[s]] <- ss_db %>% select(PC1,PC2,PC3)
        }
      base_hyp <- do.call(rbind, select_indiv)
      # hago el hypervolumen
      hyp_inv<-hypervolume(base_hyp, method='box')
      # almaceno los resultados
      com_hyp[i,"invasives"] <- hyp_inv@Volume
    
    # OVERLAP
    hyp<-hypervolume_set(hyp_nat, hyp_inv, check.memory=FALSE)
    hyp <- hypervolume_overlap_statistics(hyp)
    com_hyp[i,"overlap"] <- hyp["sorensen"]
    
    }
  
  
  # calculo el la diversidad y overlap promedio
  div_means$hyp_mean[div_means$site==com & div_means$origin=="native"] <- mean(com_hyp[,"natives"])
  div_means$hyp_SE[div_means$site==com & div_means$origin=="native"] <- SE.mean(com_hyp[,"natives"])
  div_means$hyp_mean[div_means$site==com & div_means$origin=="invasive"] <- mean(com_hyp[,"invasives"])
  div_means$hyp_SE[div_means$site==com & div_means$origin=="invasive"] <- SE.mean(com_hyp[,"invasives"])
  div_means$overlap_mean[div_means$site==com] <- mean(com_hyp[,"overlap"])
  div_means$overlap_SE[div_means$site==com] <- SE.mean(com_hyp[,"overlap"])
  
  # calulo los effect sizes
  hyp_sites[com,"EfSize"] <- cohen.d(com_hyp[,"natives"], com_hyp[,"invasives"], hedges.correction=TRUE)$estimate
  hyp_sites[com,"lowCI"] <- cohen.d(com_hyp[,"natives"], com_hyp[,"invasives"], hedges.correction=TRUE)$conf.int[1]
  hyp_sites[com,"UpCI"] <- cohen.d(com_hyp[,"natives"], com_hyp[,"invasives"], hedges.correction=TRUE)$conf.int[2]
  
}

# total ####
# genero una tabla de resultados provisional
com_hyp <- matrix(nrow=N_analyses, ncol=3)
colnames(com_hyp) <- c("natives", "invasives", "overlap")
# table(unique(mipca[,c(1,3)])$origin)
Nspp <- 33
# genero una dataframe para nativas y otro para invasoras
ss_native <- mipca %>% subset(origin=="native")
ss_invasive <- mipca %>% subset(origin=="invasive")

for (i in 1:N_analyses) {
  # NATIVAS
  # selecciono los individuos
  select_indiv <- list()
  ddd <- ss_native$species %>% unique() %>% sample(Nspp, replace=F)
  for (s in ddd){
    pre_ss_db <- ss_native %>%  filter(species==s)
    misite <- pre_ss_db$site %>% sample(1) # selecciono las especies de un unico sitio
    pre_ss_db <- pre_ss_db %>% subset(site == misite)
    ss_db <- pre_ss_db[sample( c(1:nrow(pre_ss_db)), 3, replace=F),]
    select_indiv[[s]] <- ss_db %>% select(PC1,PC2,PC3)
  }
  base_hyp <- do.call(rbind, select_indiv)
  # hago el hypervolumen
  hyp_nat<-hypervolume(base_hyp, method='box')
  # almaceno los resultados
  com_hyp[i,"natives"] <- hyp_nat@Volume
  
  # INVASORAS
  select_indiv <- list()
  ddd <- ss_invasive$species %>% unique() %>% sample(Nspp, replace=F)
  for (s in ddd){
    pre_ss_db <- ss_invasive %>%  filter(species==s)
    misite <- pre_ss_db$site %>% sample(1) # selecciono las especies de un unico sitio
    pre_ss_db <- pre_ss_db %>% subset(site == misite)
    ss_db <- pre_ss_db[sample( c(1:nrow(pre_ss_db)), 3, replace=F),]
    select_indiv[[s]] <- ss_db %>% select(PC1,PC2,PC3)
  }
  base_hyp <- do.call(rbind, select_indiv)
  # hago el hypervolumen
  hyp_inv<-hypervolume(base_hyp, method='box')
  # almaceno los resultados
  com_hyp[i,"invasives"] <- hyp_inv@Volume
  
  # OVERLAP
  hyp<-hypervolume_set(hyp_nat, hyp_inv, check.memory=FALSE)
  hyp <- hypervolume_overlap_statistics(hyp)
  com_hyp[i,"overlap"] <- hyp["sorensen"]
}

# ponemos los warnings de nuevo
options(warn=0) 
  

# calculo el la diversidad y overlap promedio
div_means$hyp_mean[div_means$site=="total" & div_means$origin=="native"] <- mean(com_hyp[,"natives"])
div_means$hyp_SE[div_means$site=="total" & div_means$origin=="native"] <- SE.mean(com_hyp[,"natives"])
div_means$hyp_mean[div_means$site=="total" & div_means$origin=="invasive"] <- mean(com_hyp[,"invasives"])
div_means$hyp_SE[div_means$site=="total" & div_means$origin=="invasive"] <- SE.mean(com_hyp[,"invasives"])
div_means$overlap_mean[div_means$site=="total"] <- mean(com_hyp[,"overlap"])
div_means$overlap_SE[div_means$site=="total"] <- SE.mean(com_hyp[,"overlap"])

# calulo los effect sizes
hyp_sites[9,"EfSize"] <- cohen.d(com_hyp[,"natives"], com_hyp[,"invasives"], hedges.correction=TRUE)$estimate
hyp_sites[9,"lowCI"] <- cohen.d(com_hyp[,"natives"], com_hyp[,"invasives"], hedges.correction=TRUE)$conf.int[1]
hyp_sites[9,"UpCI"] <- cohen.d(com_hyp[,"natives"], com_hyp[,"invasives"], hedges.correction=TRUE)$conf.int[2]
