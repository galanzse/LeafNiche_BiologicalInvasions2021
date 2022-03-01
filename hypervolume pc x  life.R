

library(hypervolume)
library(effsize)


# genero una funcion para calcular el SE
SE.mean <- function(d) { sd(d)/sqrt(length(d)) }


# cargamos los datos
source("scripts/pca.R")
colnames(data)
mipca <- mipca %>% select(species,site,origin,lifeform,PC1,PC2,PC3,sitexspp)


# genero una tabla de resultados
pc_div_life <- expand.grid(unique(mipca$lifeform), c("PC1","PC2","PC3"))
colnames(pc_div_life) <- c("lifeform","PC")
pc_div_life$hyp_mean <- NA
pc_div_life$hyp_SE <- NA

# quitamos especies con menos de 3 observaciones
nobs <- as.data.frame(table(mipca$sitexspp))
colnames(nobs) <- c("species", "Freq")
nobs <- nobs %>%  filter (Freq > 2)
mipca <- mipca[which(mipca$sitexspp %in% nobs$species),]

# table(unique(mipca[,c(1,4)])$lifeform)
Nspp <- 33
# genero una dataframe para nativas y otro para invasoras
ss_annual <- mipca %>% subset(lifeform=="annual")
ss_herbaceous <- mipca %>% subset(lifeform=="herbaceous perennial")
ss_woody <- mipca %>% subset(lifeform=="woody")

N_analyses <- 49

axes <- c("PC1","PC2","PC3")

com_list <- list()

for (i in 1:length(axes)) {
  pcx <- axes[i]
  
  # genero una tabla de resultados provisional
  com_hyp <- matrix(nrow=N_analyses, ncol=3)
  colnames(com_hyp) <- c("annual", "herbaceous perennial", "woody")
  
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
        select_indiv[[s]] <- ss_db %>% select(pcx)
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
        select_indiv[[s]] <- ss_db %>% select(pcx)
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
        select_indiv[[s]] <- ss_db %>% select(pcx)
      }
    base_hyp <- do.call(rbind, select_indiv)
    # hago el hypervolumen
    hyp<-hypervolume(base_hyp, method='box')
    # almaceno los resultados
    com_hyp[i,"woody"] <- hyp@Volume
    }
  
  com_list[[pcx]] <- com_hyp
  
  
  pc_div_life$hyp_mean[pc_div_life$lifeform=="annual" & pc_div_life$PC==pcx] <- mean(com_hyp[,"annual"])
  pc_div_life$hyp_SE[pc_div_life$lifeform=="annual" & pc_div_life$PC==pcx] <- SE.mean(com_hyp[,"annual"])
  pc_div_life$hyp_mean[pc_div_life$lifeform=="herbaceous perennial" & pc_div_life$PC==pcx] <- mean(com_hyp[,"herbaceous perennial"])
  pc_div_life$hyp_SE[pc_div_life$lifeform=="herbaceous perennial" & pc_div_life$PC==pcx] <- SE.mean(com_hyp[,"herbaceous perennial"])
  pc_div_life$hyp_mean[pc_div_life$lifeform=="woody" & pc_div_life$PC==pcx] <- mean(com_hyp[,"woody"])
  pc_div_life$hyp_SE[pc_div_life$lifeform=="woody" & pc_div_life$PC==pcx] <- SE.mean(com_hyp[,"woody"])

}
  
# ponemos los warnings de nuevo
options(warn=0) 



# figura 

eff_life_pc <- expand.grid(c("annual - herbaceous perennial","annual - woody","herbaceous perennial - woody"), c("PC1","PC2","PC3"))
colnames(eff_life_pc) <- c("group","PC")
eff_life_pc$EfSize <- NA
eff_life_pc$lowCI <- NA
eff_life_pc$UpCI <- NA

com_list[[pcx]] <- com_hyp
str(com_list)

com_hyp <- com_list$PC1
eff_life_pc$EfSize[eff_life_pc$group=="annual - herbaceous perennial" & eff_life_pc$PC=="PC1"] <- cohen.d(com_hyp[,"annual"], com_hyp[,"herbaceous perennial"], hedges.correction=TRUE)$estimate
eff_life_pc$lowCI[eff_life_pc$group=="annual - herbaceous perennial" & eff_life_pc$PC=="PC1"] <- cohen.d(com_hyp[,"annual"], com_hyp[,"herbaceous perennial"], hedges.correction=TRUE)$conf.int[1]
eff_life_pc$UpCI[eff_life_pc$group=="annual - herbaceous perennial" & eff_life_pc$PC=="PC1"]  <- cohen.d(com_hyp[,"annual"], com_hyp[,"herbaceous perennial"], hedges.correction=TRUE)$conf.int[2]
eff_life_pc$EfSize[eff_life_pc$group=="annual - woody" & eff_life_pc$PC=="PC1"] <- cohen.d(com_hyp[,"annual"], com_hyp[,"woody"], hedges.correction=TRUE)$estimate
eff_life_pc$lowCI[eff_life_pc$group=="annual - woody" & eff_life_pc$PC=="PC1"] <- cohen.d(com_hyp[,"annual"], com_hyp[,"woody"], hedges.correction=TRUE)$conf.int[1]
eff_life_pc$UpCI[eff_life_pc$group=="annual - woody" & eff_life_pc$PC=="PC1"]  <- cohen.d(com_hyp[,"annual"], com_hyp[,"woody"], hedges.correction=TRUE)$conf.int[2]
eff_life_pc$EfSize[eff_life_pc$group=="herbaceous perennial - woody" & eff_life_pc$PC=="PC1"] <- cohen.d(com_hyp[,"herbaceous perennial"], com_hyp[,"woody"], hedges.correction=TRUE)$estimate
eff_life_pc$lowCI[eff_life_pc$group=="herbaceous perennial - woody" & eff_life_pc$PC=="PC1"] <- cohen.d(com_hyp[,"herbaceous perennial"], com_hyp[,"woody"], hedges.correction=TRUE)$conf.int[1]
eff_life_pc$UpCI[eff_life_pc$group=="herbaceous perennial - woody" & eff_life_pc$PC=="PC1"]  <- cohen.d(com_hyp[,"herbaceous perennial"], com_hyp[,"woody"], hedges.correction=TRUE)$conf.int[2]

com_hyp <- com_list$PC2
eff_life_pc$EfSize[eff_life_pc$group=="annual - herbaceous perennial" & eff_life_pc$PC=="PC2"] <- cohen.d(com_hyp[,"annual"], com_hyp[,"herbaceous perennial"], hedges.correction=TRUE)$estimate
eff_life_pc$lowCI[eff_life_pc$group=="annual - herbaceous perennial" & eff_life_pc$PC=="PC2"] <- cohen.d(com_hyp[,"annual"], com_hyp[,"herbaceous perennial"], hedges.correction=TRUE)$conf.int[1]
eff_life_pc$UpCI[eff_life_pc$group=="annual - herbaceous perennial" & eff_life_pc$PC=="PC2"]  <- cohen.d(com_hyp[,"annual"], com_hyp[,"herbaceous perennial"], hedges.correction=TRUE)$conf.int[2]
eff_life_pc$EfSize[eff_life_pc$group=="annual - woody" & eff_life_pc$PC=="PC2"] <- cohen.d(com_hyp[,"annual"], com_hyp[,"woody"], hedges.correction=TRUE)$estimate
eff_life_pc$lowCI[eff_life_pc$group=="annual - woody" & eff_life_pc$PC=="PC2"] <- cohen.d(com_hyp[,"annual"], com_hyp[,"woody"], hedges.correction=TRUE)$conf.int[1]
eff_life_pc$UpCI[eff_life_pc$group=="annual - woody" & eff_life_pc$PC=="PC2"]  <- cohen.d(com_hyp[,"annual"], com_hyp[,"woody"], hedges.correction=TRUE)$conf.int[2]
eff_life_pc$EfSize[eff_life_pc$group=="herbaceous perennial - woody" & eff_life_pc$PC=="PC2"] <- cohen.d(com_hyp[,"herbaceous perennial"], com_hyp[,"woody"], hedges.correction=TRUE)$estimate
eff_life_pc$lowCI[eff_life_pc$group=="herbaceous perennial - woody" & eff_life_pc$PC=="PC2"] <- cohen.d(com_hyp[,"herbaceous perennial"], com_hyp[,"woody"], hedges.correction=TRUE)$conf.int[1]
eff_life_pc$UpCI[eff_life_pc$group=="herbaceous perennial - woody" & eff_life_pc$PC=="PC2"]  <- cohen.d(com_hyp[,"herbaceous perennial"], com_hyp[,"woody"], hedges.correction=TRUE)$conf.int[2]

com_hyp <- com_list$PC3
eff_life_pc$EfSize[eff_life_pc$group=="annual - herbaceous perennial" & eff_life_pc$PC=="PC3"] <- cohen.d(com_hyp[,"annual"], com_hyp[,"herbaceous perennial"], hedges.correction=TRUE)$estimate
eff_life_pc$lowCI[eff_life_pc$group=="annual - herbaceous perennial" & eff_life_pc$PC=="PC3"] <- cohen.d(com_hyp[,"annual"], com_hyp[,"herbaceous perennial"], hedges.correction=TRUE)$conf.int[1]
eff_life_pc$UpCI[eff_life_pc$group=="annual - herbaceous perennial" & eff_life_pc$PC=="PC3"]  <- cohen.d(com_hyp[,"annual"], com_hyp[,"herbaceous perennial"], hedges.correction=TRUE)$conf.int[2]
eff_life_pc$EfSize[eff_life_pc$group=="annual - woody" & eff_life_pc$PC=="PC3"] <- cohen.d(com_hyp[,"annual"], com_hyp[,"woody"], hedges.correction=TRUE)$estimate
eff_life_pc$lowCI[eff_life_pc$group=="annual - woody" & eff_life_pc$PC=="PC3"] <- cohen.d(com_hyp[,"annual"], com_hyp[,"woody"], hedges.correction=TRUE)$conf.int[1]
eff_life_pc$UpCI[eff_life_pc$group=="annual - woody" & eff_life_pc$PC=="PC3"]  <- cohen.d(com_hyp[,"annual"], com_hyp[,"woody"], hedges.correction=TRUE)$conf.int[2]
eff_life_pc$EfSize[eff_life_pc$group=="herbaceous perennial - woody" & eff_life_pc$PC=="PC3"] <- cohen.d(com_hyp[,"herbaceous perennial"], com_hyp[,"woody"], hedges.correction=TRUE)$estimate
eff_life_pc$lowCI[eff_life_pc$group=="herbaceous perennial - woody" & eff_life_pc$PC=="PC3"] <- cohen.d(com_hyp[,"herbaceous perennial"], com_hyp[,"woody"], hedges.correction=TRUE)$conf.int[1]
eff_life_pc$UpCI[eff_life_pc$group=="herbaceous perennial - woody" & eff_life_pc$PC=="PC3"]  <- cohen.d(com_hyp[,"herbaceous perennial"], com_hyp[,"woody"], hedges.correction=TRUE)$conf.int[2]












