
library(tidyverse)
library(readxl)


# Quantifying the relevance of inter- intra- life forms variability
# De Bello et al. 2011


# cargamos los datos
source("scripts/pca.R")
colnames(mipca)
mipca <- mipca[,c("species","site","origin","lifeform","PC1","PC2","PC3")]



# creo una tabla de resultados
deBello.ori <- read_excel("resultados.xlsx", sheet = "de Bello origin", range = "A1:E145")



# creo una matriz de medias
mipca.mean <- mipca %>% group_by(site, species, origin) %>%
  summarise(PC1=mean(PC1),
            PC2=mean(PC2),
            PC3=mean(PC3))
mipca.mean <- as.data.frame(mipca.mean)



# primero: calculo la varianza total asociada a cada PC x origin en cada sitio ###

for (com in unique(deBello.ori$site)) {
  # selecciono una comunidad
  micom <- subset(mipca, site==com)
  
  for (ori in c("native", "invasive")) {
    # selecciono nativas o invasoras
    miori <- subset(micom, origin==ori)
    # selecciono las observaciones de la matiz de medias para calcular Xcom
    mimean <- mipca.mean %>% subset(site==com & origin==ori)

    for (pc in c("PC1", "PC2", "PC3")) {
      # calculo la media del trait que me interesa
      Xcom <- mimean[,pc] %>% mean()
      # miro el numero de especies
      Nsp <- miori$species %>% unique() %>% length()
      # vector de resultados
      vspp <- vector()
      
      for (spp in unique(miori$species)) {
        # selecciono una especie
        myspp <- subset(miori, species == spp)
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
      deBello.ori$variance[deBello.ori$site==com & deBello.ori$origin==ori & deBello.ori$group=="total" & deBello.ori$pc==pc] <- sum(vspp)
      
    }
  }
}



# segundo: calculo la variabilidad interespecifica

for (com in unique(deBello.ori$site)) {
  # selecciono una comunidad
  micom <- subset(mipca, site==com)
  
  for (ori in c("native", "invasive")) {
    # selecciono nativas o invasoras
    miori <- subset(micom, origin==ori)
    # selecciono las observaciones de la matiz de medias para calcular Xcom
    mimean <- mipca.mean %>% subset(site==com & origin==ori)
    
    for (pc in c("PC1", "PC2", "PC3")) {
      # calculo la media del trait que me interesa
      Xcom <- mimean[,pc] %>% mean()
      # miro el numero de especies
      Nsp <- miori$species %>% unique() %>% length()
      # vector de resultados
      vspp <- vector()
      
      for (spp in unique(miori$species)) {
        # selecciono una especie
        myspp <- subset(miori, species == spp)
        # calculo la media de ese trait para mi especie
        # y almaceno el resultado
        vspp <- c(vspp, (mean(myspp[,pc]) - Xcom)^2/Nsp)
      }
      
      # almaceno la varianza total en la tabla de resultados
      deBello.ori$variance[deBello.ori$site==com & deBello.ori$origin==ori & deBello.ori$pc==pc & deBello.ori$group=="interspecific"] <- sum(vspp)
      
    }
  }
}



# tercero: calculo la variabilidad intraspecifica

for (com in unique(deBello.ori$site)) {
  # selecciono una comunidad
  micom <- subset(mipca, site==com)
  
  for (ori in c("native", "invasive")) {
    # selecciono nativas o invasoras
    miori <- subset(micom, origin==ori)
    
    for (pc in c("PC1", "PC2", "PC3")) {
      # miro el numero de especies
      Nsp <- miori$species %>% unique() %>% length()
      # vector de resultados
      vspp <- vector()
      
      for (spp in unique(miori$species)) {
        # selecciono una especie
        myspp <- subset(miori, species == spp)
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
      deBello.ori$variance[deBello.ori$site==com & deBello.ori$origin==ori & deBello.ori$group=="intraspecific" & deBello.ori$pc==pc] <- sum(vspp)
      
    }
  }
}



# compruebo que lo he hecho bien
for (com in unique(deBello.ori$site)) {
  for (ori in unique(deBello.ori$origin)) {
    for (tr in unique(deBello.ori$pc)) {
      mi_ss <- deBello.ori %>% subset(site==com & origin==ori & pc==tr & group %in% c("interspecific","intraspecific")) 
      mi_tot <- deBello.ori %>% subset(site==com & origin==ori & pc==tr & group %in% c("total")) 
      print(paste(round(sum(mi_ss$variance), 4)==round(mi_tot$variance, 4)))
    }
  }
}







