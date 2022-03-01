
library(tidyverse)
library(readxl)


# Quantifying the relevance of inter- intra- life forms variability
# De Bello et al. 2011

# ORIGIN ####

# cargamos los datos
source("scripts/pca.R")
colnames(mipca)
mipca <- mipca[,c("species","site","origin","lifeform","PC1","PC2","PC3")]


# creo una tabla de resultados
dB.origin.total <- read_excel("resultados.xlsx", sheet = "Bello origin total")


# creo una matriz de medias
mipca.mean <- mipca %>% group_by(site, species, origin) %>%
  summarise(PC1=mean(PC1),
            PC2=mean(PC2),
            PC3=mean(PC3))
mipca.mean <- as.data.frame(mipca.mean)


# primero: calculo la varianza total asociada a cada PC x origin en cada sitio ###
for (ori in c("native", "invasive")) {
  # selecciono nativas o invasoras
  miori <- subset(mipca, origin==ori)
  # selecciono las observaciones de la matiz de medias para calcular Xcom
  mimean <- mipca.mean %>% subset(origin==ori)
  
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
    dB.origin.total$variance[dB.origin.total$origin==ori & dB.origin.total$pc==pc& dB.origin.total$group=="total"] <- sum(vspp)
    
  }
}


# segundo: calculo la variabilidad interespecifica
for (ori in c("native", "invasive")) {
  # selecciono nativas o invasoras
  miori <- subset(mipca, origin==ori)
  # selecciono las observaciones de la matiz de medias para calcular Xcom
  mimean <- mipca.mean %>% subset(origin==ori)
  
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
    dB.origin.total$variance[dB.origin.total$origin==ori & dB.origin.total$pc==pc & dB.origin.total$group=="interspecific"] <- sum(vspp)
    
  }
}


# tercero: calculo la variabilidad intraspecifica
for (ori in c("native", "invasive")) {
  # selecciono nativas o invasoras
  miori <- subset(mipca, origin==ori)
  
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
    dB.origin.total$variance[dB.origin.total$origin==ori & dB.origin.total$group=="intraspecific" & dB.origin.total$pc==pc] <- sum(vspp)
    
  }
}


# compruebo que lo he hecho bien
for (ori in unique(dB.origin.total$origin)) {
  for (tr in unique(dB.origin.total$pc)) {
    mi_ss <- dB.origin.total %>% subset(origin==ori & pc==tr & group %in% c("interspecific","intraspecific")) 
    mi_tot <- dB.origin.total %>% subset(origin==ori & pc==tr & group %in% c("total")) 
    print(paste(round(sum(mi_ss$variance), 4)==round(mi_tot$variance, 4)))
  }
}



# LIFE FORM ####

# creo una tabla de resultados
dB.life.total <- read_excel("resultados.xlsx", sheet = "Bello life total")

# creo una matriz de medias
mipca.mean <- mipca %>% group_by(site, species, lifeform) %>%
  summarise(PC1=mean(PC1),
            PC2=mean(PC2),
            PC3=mean(PC3))
mipca.mean <- as.data.frame(mipca.mean)


# primero: calculo la varianza total asociada a cada PC x lifeform en cada sitio ###
for (ori in c("woody", "annual", "herbaceous perennial")) {
  # selecciono nativas o invasoras
  miori <- subset(mipca, lifeform==ori)
  # selecciono las observaciones de la matiz de medias para calcular Xcom
  mimean <- mipca.mean %>% subset(lifeform==ori)
  
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
    dB.life.total$variance[dB.life.total$life==ori & dB.life.total$pc==pc& dB.life.total$group=="total"] <- sum(vspp)
    
  }
}


# segundo: calculo la variabilidad interespecifica
for (ori in c("woody", "annual", "herbaceous perennial")) {
  # selecciono nativas o invasoras
  miori <- subset(mipca, lifeform==ori)
  # selecciono las observaciones de la matiz de medias para calcular Xcom
  mimean <- mipca.mean %>% subset(lifeform==ori)
  
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
    dB.life.total$variance[dB.life.total$life==ori & dB.life.total$pc==pc & dB.life.total$group=="interspecific"] <- sum(vspp)
    
  }
}


# tercero: calculo la variabilidad intraspecifica
for (ori in c("woody", "annual", "herbaceous perennial")) {
  # selecciono nativas o invasoras
  miori <- subset(mipca, lifeform==ori)
  
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
    dB.life.total$variance[dB.life.total$life==ori & dB.life.total$group=="intraspecific" & dB.life.total$pc==pc] <- sum(vspp)
    
  }
}


# compruebo que lo he hecho bien
for (ori in unique(dB.life.total$life)) {
  for (tr in unique(dB.life.total$pc)) {
    mi_ss <- dB.life.total %>% subset(life==ori & pc==tr & group %in% c("interspecific","intraspecific")) 
    mi_tot <- dB.life.total %>% subset(life==ori & pc==tr & group %in% c("total")) 
    print(paste(round(sum(mi_ss$variance), 4)==round(mi_tot$variance, 4)))
  }
}



# PLOTS ####

dB.origin.total
invxori <- dB.origin.total %>% subset(group %in% c("total","intraspecific")) %>% spread("group", "variance")
invxori$int2total <- invxori$intraspecific / invxori$total

origin.bar <- matrix(ncol=3, nrow=6) %>% as.data.frame()
colnames(origin.bar) <- c("PC","origin","intraspecific")
origin.bar$PC <- c("PC1","PC1","PC2","PC2","PC3","PC3")
origin.bar$origin <- c(rep(c("invasive","native"),3))

origin.bar[1,c("intraspecific")] <- 0.254
origin.bar[2,c("intraspecific")] <- 0.107
origin.bar[3,c("intraspecific")] <- 0.239
origin.bar[4,c("intraspecific")] <- 0.185
origin.bar[5,c("intraspecific")] <- 0.266
origin.bar[6,c("intraspecific")] <- 0.189

pd <- position_dodge(width=0.7)
a <- ggplot(data = origin.bar, aes(x = PC, y = intraspecific, group=origin, fill=origin)) +
  geom_bar(stat="identity", position=pd, width = 0.7, colour="black") +
  scale_y_continuous(limits=c(0,0.4), labels = scales::number_format(accuracy = 0.01)) +
  theme_bw() +
  xlab(NULL) + ylab("intraspecific / total variance") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text=element_text(size=12),
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.title = element_text(size=14)) +
  scale_fill_manual(values = c("#F8766D","#00BFC4"))



# barplot de lifeform
dB.life.total
invxlife <- dB.life.total %>% subset(group %in% c("total","intraspecific")) %>% spread("group", "variance")
invxlife$int2total <- invxlife$intraspecific / invxlife$total

life.bar <- matrix(ncol=4, nrow=9) %>% as.data.frame()
colnames(life.bar) <- c("PC","life","intraspecific")
life.bar$PC <- c("PC1","PC1","PC1","PC2","PC2","PC2","PC3","PC3","PC3")
life.bar$life <- c(rep(c("annual","herbaceous perennial", "woody"),3))

life.bar[1,c("intraspecific")] <- 0.313
life.bar[2,c("intraspecific")] <- 0.112
life.bar[3,c("intraspecific")] <- 0.105
life.bar[4,c("intraspecific")] <- 0.277
life.bar[5,c("intraspecific")] <- 0.209
life.bar[6,c("intraspecific")] <- 0.156
life.bar[7,c("intraspecific")] <- 0.370
life.bar[8,c("intraspecific")] <- 0.232
life.bar[9,c("intraspecific")] <- 0.197

pd <- position_dodge(width=0.7)
b <- ggplot(data = life.bar, aes(x = PC, y = intraspecific, group=life, fill=life)) +
  geom_bar(stat="identity", position=pd, width = 0.7, colour="black") +
  scale_y_continuous(limits=c(0,0.4), labels = scales::number_format(accuracy = 0.01)) +
  theme_bw() +
  xlab(NULL) + ylab("intraspecific / total variance") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text=element_text(size=12),
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.title = element_text(size=14)) +
  scale_fill_manual(values = c("thistle4","palegreen2","khaki3"))


ggarrange(a, b, labels=c("A","B"))




