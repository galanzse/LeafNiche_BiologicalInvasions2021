
library(tidyverse)
library(readxl)


# Quantifying the relevance of inter- intra- life forms variability
# De Bello et al. 2011


# cargamos los datos
source("scripts/pca.R")
colnames(mipca)
mipca <- mipca[,c("species","site","origin","lifeform","PC1","PC2","PC3")]



# creo una tabla de resultados
deBello.life <- read_excel("resultados.xlsx", sheet = "de Bello life", range = "A1:E217")



# creo una matriz de medias
mipca.mean <- mipca %>% group_by(site, species, lifeform) %>%
  summarise(PC1=mean(PC1),
            PC2=mean(PC2),
            PC3=mean(PC2))
mipca.mean <- as.data.frame(mipca.mean)



# primero: calculo la varianza total asociada a cada PC x lifeform en cada sitio ###

for (com in unique(deBello.life$site)) {
  # selecciono una comunidad
  micom <- subset(mipca, site==com)
  
  for (lf in c("annual", "herbaceous perennial", "woody")) {
    # selecciono nativas o invasoras
    milife <- subset(micom, lifeform==lf)
    # selecciono las observaciones de la matiz de medias para calcular Xcom
    mimean <- mipca.mean %>% subset(site==com & lifeform==lf)
    
    for (pc in c("PC1", "PC2", "PC3")) {
      # calculo la media del trait que me interesa
      Xcom <- mimean[,pc] %>% mean()
      # miro el numero de especies
      Nsp <- milife$species %>% unique() %>% length()
      # vector de resultados
      vspp <- vector()
      
      for (spp in unique(milife$species)) {
        # selecciono una especie
        myspp <- subset(milife, species == spp)
        # numero de observaciones de la especie
        Nind <- dim(myspp)[1]
        # distancia de las observaciones a la media de la comunidad
        vobs <- vector()
        
        for (obs in rownames(myspp)) {
          # calculo la dispersion de las observaciones respecto a la media del pool de nativas o invasoras
          vobs[obs] <- (myspp[obs,pc] - Xcom)^2/Nind
        }
        vspp <- c(vspp, sum(vobs)/Nsp)
        
      }
      
      # almaceno la varianza total en la tabla de resultados
      deBello.life$variance[deBello.life$site==com & deBello.life$life==lf & deBello.life$group=="total" & deBello.life$pc==pc] <- sum(vspp)
      
    }
  }
}



# segundo: calculo la variabilidad interespecifica

for (com in unique(deBello.life$site)) {
  # selecciono una comunidad
  micom <- subset(mipca, site==com)
  
  for (lf in c("annual", "herbaceous perennial", "woody")) {
    # selecciono nativas o invasoras
    milife <- subset(micom, lifeform==lf)
    # selecciono las observaciones de la matiz de medias para calcular Xcom
    mimean <- mipca.mean %>% subset(site==com & lifeform==lf)
    
    for (pc in c("PC1", "PC2", "PC3")) {
      # calculo la media del trait que me interesa
      Xcom <- mimean[,pc] %>% mean()
      # miro el numero de especies
      Nsp <- milife$species %>% unique() %>% length()
      # vector de resultados
      vspp <- vector()
      
      for (spp in unique(milife$species)) {
        # selecciono una especie
        myspp <- subset(milife, species == spp)
        # calculo la media de ese trait para mi especie
        # y almaceno el resultado
        vspp <- c(vspp, (mean(myspp[,pc]) - Xcom)^2/Nsp)
      }
      
      # almaceno la varianza total en la tabla de resultados
      deBello.life$variance[deBello.life$site==com & deBello.life$life==lf & deBello.life$pc==pc & deBello.life$group=="interspecific"] <- sum(vspp)
      
    }
  }
}



# tercero: calculo la variabilidad intraspecifica

for (com in unique(deBello.life$site)) {
  # selecciono una comunidad
  micom <- subset(mipca, site==com)
  
  for (lf in c("annual", "herbaceous perennial", "woody")) {
    # selecciono nativas o invasoras
    milife <- subset(micom, lifeform==lf)
    
    for (pc in c("PC1", "PC2", "PC3")) {
      # miro el numero de especies
      Nsp <- milife$species %>% unique() %>% length()
      # vector de resultados
      vspp <- vector()
      
      for (spp in unique(milife$species)) {
        # selecciono una especie
        myspp <- subset(milife, species == spp)
        # numero de observaciones de la especie
        Nind <- dim(myspp)[1]
        # calculo la media para pc de la especie
        Xi <- mean(myspp[,pc])
        # distancia de las observaciones a la media de la comunidad
        vobs <- vector()
        
        for (obs in rownames(myspp)) {
          # calculo la dispersion de las observaciones respecto a la media del pool de nativas o invasoras
          vobs[obs] <- (myspp[obs,pc] - Xi)^2/Nind
        }
        vspp <- c(vspp, sum(vobs)/Nsp)
        
      }
      
      # almaceno la varianza total en la tabla de resultados
      deBello.life$variance[deBello.life$site==com & deBello.life$life==lf & deBello.life$group=="intraspecific" & deBello.life$pc==pc] <- sum(vspp)
      
    }
  }
}



# compruebo que lo he hecho bien
for (com in unique(deBello.life$site)) {
  for (lf in unique(deBello.life$life)) {
    for (tr in unique(deBello.life$pc)) {
      mi_ss <- deBello.life %>% subset(site==com & life==lf & pc==tr & group %in% c("interspecific","intraspecific")) 
      mi_tot <- deBello.life %>% subset(site==com & life==lf & pc==tr & group %in% c("total")) 
      print(paste(round(sum(mi_ss$variance), 4)==round(mi_tot$variance, 4)))
    }
  }
}
