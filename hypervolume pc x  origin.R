

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
pc_div_ori <- expand.grid(unique(mipca$origin), c("PC1","PC2","PC3"))
colnames(pc_div_ori) <- c("origin","PC")
pc_div_ori$overlap_mean <- NA
pc_div_ori$overlap_SE <- NA
pc_div_ori$hyp_mean <- NA
pc_div_ori$hyp_SE <- NA


# tabla de resultados de hipervolumenes
eff_origin <- matrix(nrow=length(c("PC1","PC2","PC3")), ncol=3)
rownames(eff_origin) <- c("PC1","PC2","PC3")
colnames(eff_origin) <- c("EfSize","lowCI","UpCI")


# establecemos el numero de replicas por especie, numerod de especies por sitio,
# y numero de analisis
N_analyses <- 49
N_replicates <- 4
# table(unique(mipca[,c(1,3)])$origin)
Nspp <- 33


# eliminamos warnings (Quique me dijo que no hay problema con que
# Log number of observations is less than or equal to the number of dimensions, siempre que eset cerca)
options(warn=-1)

# genero una dataframe para nativas y otro para invasoras
ss_native <- mipca %>% subset(origin=="native")
ss_invasive <- mipca %>% subset(origin=="invasive")

axes <- c("PC1","PC2","PC3")

for (i in 1:length(axes)) {
  pcx <- axes[i]
  
  # genero una tabla de resultados provisional
  com_hyp <- matrix(nrow=N_analyses, ncol=3)
  colnames(com_hyp) <- c("native", "invasive", "overlap")
  
  
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
        select_indiv[[s]] <- ss_db %>% select(pcx)
      }
      base_hyp <- do.call(rbind, select_indiv)
      # hago el hypervolumen
      hyp_nat<-hypervolume(base_hyp, method='box')
      # almaceno los resultados
      com_hyp[i,"native"] <- hyp_nat@Volume
      
      # INVASORAS
      select_indiv <- list()
      ddd <- ss_invasive$species %>% unique() %>% sample(Nspp, replace=F)
      for (s in ddd){
        pre_ss_db <- ss_invasive %>%  filter(species==s)
        misite <- pre_ss_db$site %>% sample(1) # selecciono las especies de un unico sitio
        pre_ss_db <- pre_ss_db %>% subset(site == misite)
        ss_db <- pre_ss_db[sample( c(1:nrow(pre_ss_db)), 3, replace=F),]
        select_indiv[[s]] <- ss_db %>% select(pcx)
      }
      base_hyp <- do.call(rbind, select_indiv)
      # hago el hypervolumen
      hyp_inv<-hypervolume(base_hyp, method='box')
      # almaceno los resultados
      com_hyp[i,"invasive"] <- hyp_inv@Volume
      
      # OVERLAP
      hyp<-hypervolume_set(hyp_nat, hyp_inv, check.memory=FALSE)
      hyp <- hypervolume_overlap_statistics(hyp)
      com_hyp[i,"overlap"] <- hyp["sorensen"]
    }
  
  pc_div_ori$hyp_mean[pc_div_ori$PC==pcx & pc_div_ori$origin=="native"] <- mean(com_hyp[,"native"])
  pc_div_ori$hyp_SE[pc_div_ori$PC==pcx & pc_div_ori$origin=="native"] <- SE.mean(com_hyp[,"native"])
  pc_div_ori$hyp_mean[pc_div_ori$PC==pcx & pc_div_ori$origin=="invasive"] <- mean(com_hyp[,"invasive"])
  pc_div_ori$hyp_SE[pc_div_ori$PC==pcx & pc_div_ori$origin=="invasive"] <- SE.mean(com_hyp[,"invasive"])
  pc_div_ori$overlap_mean[pc_div_ori$PC==pcx] <- mean(com_hyp[,"overlap"])
  pc_div_ori$overlap_SE[pc_div_ori$PC==pcx] <- SE.mean(com_hyp[,"overlap"])
  
  eff_origin[pcx,"EfSize"] <- cohen.d(com_hyp[,"native"], com_hyp[,"invasive"], hedges.correction=TRUE)$estimate
  eff_origin[pcx,"lowCI"] <- cohen.d(com_hyp[,"native"], com_hyp[,"invasive"], hedges.correction=TRUE)$conf.int[1]
  eff_origin[pcx,"UpCI"] <- cohen.d(com_hyp[,"native"], com_hyp[,"invasive"], hedges.correction=TRUE)$conf.int[2]

}

# ponemos los warnings de nuevo
options(warn=0) 

eff_origin <- eff_origin %>% as.data.frame()
eff_origin$PC <- rownames(eff_origin)
ggplot(data=eff_origin) +
  geom_point(aes(y = PC, x = EfSize), size=4, shape=19, fill="white") +
  geom_errorbarh(aes(y = PC, xmin = lowCI, xmax = UpCI), height=0.3, size=.7, color="black") +
  theme_bw() +
  geom_vline(xintercept = 0) +
  xlab("Hedges' d") + ylab(" ") +
  theme(axis.title.x = element_text(color="black", size=14),
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12))



