

source("scripts/pca.R")

library(effsize)

colnames(data)
mipca <- mipca %>% select(species,site,origin,lifeform,PC1,PC2,PC3,sitexspp)


# quitamos especies con menos de 3 observaciones
nobs <- as.data.frame(table(mipca$sitexspp))
colnames(nobs) <- c("species", "Freq")
nobs <- nobs %>%  filter (Freq > 3)
mipca <- mipca[which(mipca$sitexspp %in% nobs$species),]


# establecemos el numero de replicas por especie, numero de especies por sitio,
# y numero de analisis
N_analyses <- 99
N_replicates <- 4


# GLOBAL ####

# table(unique(mipca[,c(1,3)])$origin)
Nspp <- 33


# tabla de resultados
mpd_matrix <- matrix(ncol=2, nrow=N_analyses)
colnames(mpd_matrix) <- c("native", "invasive")
  

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
      select_indiv[[s]] <- ss_db
    }
    base_mpd <- do.call(rbind, select_indiv)
    # mpd
    mpd_matrix[i,"native"] <- base_mpd %>% select(PC1,PC2,PC3) %>% dist() %>% mean()

    
    # INVASORAS
    select_indiv <- list()
    ddd <- ss_invasive$species %>% unique() %>% sample(Nspp, replace=F)
    for (s in ddd){
      pre_ss_db <- ss_invasive %>%  filter(species==s)
      misite <- pre_ss_db$site %>% sample(1) # selecciono las especies de un unico sitio
      pre_ss_db <- pre_ss_db %>% subset(site == misite)
      ss_db <- pre_ss_db[sample( c(1:nrow(pre_ss_db)), 3, replace=F),]
      select_indiv[[s]] <- ss_db
    }
    base_mpd <- do.call(rbind, select_indiv)
    # hago el mpdervolumen
    mpd_matrix[i,"invasive"] <- base_mpd %>% select(PC1,PC2,PC3) %>% dist() %>% mean()
    
    print(paste(i,"out of",N_analyses))
    
}


# resultados
boxplot(mpd_matrix)
# las invasoras son significativamente mas diversas
cohen.d(mpd_matrix[,"native"], mpd_matrix[,"invasive"], hedges.correction=TRUE)



# SITIOS ####


# lista de resultados
mpd_site_list <- list()


# es necesario ajustar el numero de especies en cada grupo de nativas e invasoras por sitio
N_species <- matrix(nrow=length(unique(mipca$site)), ncol=2) %>% as.data.frame()
colnames(N_species) <- c("site", "N_spp")
N_species$site <- unique(mipca$site)
table(unique(mipca[,1:3])$origin, unique(mipca[,1:3])$site)
N_species$N_spp <- c(9,5,8,5,4,3,9,4)

for (com in unique(mipca$site)) {
  
  # establezco del numero de especies por cada grupo de nativas e invasoras
  N_spp <- N_species$N_spp[N_species$site==com]
  
  # genero una dataframe para nativas y otro para invasoras
  ss_native <- mipca %>% subset(site==com & origin=="native")
  ss_invasive <- mipca %>% subset(site==com & origin=="invasive")
  
  # genero una tabla de resultados provisional
  mpd_com <- matrix(ncol=2, nrow=N_analyses)
  colnames(mpd_com) <- c("native", "invasive")

    for (i in 1:N_analyses) {
      # NATIVAS
      # selecciono los individuos
      select_indiv <- list()
      ddd <- ss_native$species %>% unique() %>% sample(N_spp, replace=F)
      for (s in ddd){
        pre_ss_db <- ss_native %>%  filter(species==s)
        misite <- pre_ss_db$site %>% sample(1) # selecciono las especies de un unico sitio
        pre_ss_db <- pre_ss_db %>% subset(site == misite)
        ss_db <- pre_ss_db[sample( c(1:nrow(pre_ss_db)), 3, replace=F),]
        select_indiv[[s]] <- ss_db
      }
      base_mpd <- do.call(rbind, select_indiv)
      # mpd
      mpd_com[i,"native"] <- base_mpd %>% select(PC1,PC2,PC3) %>% dist() %>% mean()
      
      
      # INVASORAS
      select_indiv <- list()
      ddd <- ss_invasive$species %>% unique() %>% sample(N_spp, replace=F)
      for (s in ddd){
        pre_ss_db <- ss_invasive %>%  filter(species==s)
        misite <- pre_ss_db$site %>% sample(1) # selecciono las especies de un unico sitio
        pre_ss_db <- pre_ss_db %>% subset(site == misite)
        ss_db <- pre_ss_db[sample( c(1:nrow(pre_ss_db)), 3, replace=F),]
        select_indiv[[s]] <- ss_db
      }
      base_mpd <- do.call(rbind, select_indiv)
      # hago el mpdervolumen
      mpd_com[i,"invasive"] <- base_mpd %>% select(PC1,PC2,PC3) %>% dist() %>% mean()

    }
  
  mpd_site_list[[com]] <- mpd_com

  print(com)
  
}


# resultados
for (com in unique(mipca$site)) {
  mycom <- mpd_site_list[[com]]
  print(com)
  print(cohen.d(mycom[,"native"], mycom[,"invasive"], hedges.correction=TRUE))
}




