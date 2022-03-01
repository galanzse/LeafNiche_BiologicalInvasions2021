
library(tidyverse)

# Quantifying the relevance of intraspecific trait variability for functional diversity
# De Bello et al. 2011



# cargamos los datos
source("scripts/pca.R")



# como la varianza puede depender del sitio vamos a seleccionar por especies y sitio
mipca$SpeciesSite <- paste(mipca$species,mipca$site)



# variance ####

variances_spp <- matrix(ncol=3, nrow=length(unique(mipca$SpeciesSite)))
rownames(variances_spp) <- unique(mipca$SpeciesSite)
colnames(variances_spp) <- c('PC1','PC2','PC3')

# calculamos la varianza para los pares unicos de especies y sitio

for (sp in unique(mipca$SpeciesSite)) {
  
  pre_ss_db <- subset(mipca, SpeciesSite==sp)
  
  variances_spp[sp,'PC1'] <- var(pre_ss_db[,'PC1'])
  variances_spp[sp,'PC2'] <- var(pre_ss_db[,'PC2'])
  variances_spp[sp,'PC3'] <- var(pre_ss_db[,'PC3'])
  
  # rm(pre_ss_db, sp)
  
}



# pegamos la tabla de varianzas a la original

variances_spp <- as.data.frame(variances_spp)
variances_spp$SpeciesSite <- rownames(variances_spp)

# quito los valores originales de PC en mipca porque se van a confundir con las varianzas de variances_spp
mipca$PC1 <- NULL
mipca$PC2 <- NULL
mipca$PC3 <- NULL

variances_spp <- left_join(variances_spp, unique(mipca), by='SpeciesSite')

variances_spp$site <- as.factor(variances_spp$site)

levels(variances_spp$site) <- c("ASF","BW","CBW","CG","CSS","R","SW","SG")

par(mfrow=c(3,1), mar=c(2,5,1,2))
boxplot(variances_spp$PC1 ~ variances_spp$site, ylab="ITV PC1", xlab=NULL)
boxplot(variances_spp$PC2 ~ variances_spp$site, ylab="ITV PC2", xlab=NULL)
boxplot(variances_spp$PC3 ~ variances_spp$site, ylab="ITV PC3", xlab=NULL)


lmer(PC1 ~ site + (1|species), data=variances_spp) %>% anova()
lmer(PC2 ~ site + (1|species), data=variances_spp) %>% anova()
lmer(PC3 ~ site + (1|species), data=variances_spp) %>% anova()




