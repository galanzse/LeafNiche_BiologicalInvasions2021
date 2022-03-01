

library(tidyverse)
library(readxl)
library(ggplot2)
library(ggfortify)
library(lmerTest)


# cARGO Y CURO LOS DATOS ####

replicates <- read_excel("C:/Users/Javier/OneDrive/TESIS Y PUBLICACIONES/TESIS/2. Leaf-Niche Funk/Analyses Biological Invasions 250320/data.xlsx", 
                   sheet = "replicates")
str(replicates)
factors <- read_excel("C:/Users/Javier/OneDrive/TESIS Y PUBLICACIONES/TESIS/2. Leaf-Niche Funk/Analyses Biological Invasions 250320/data.xlsx", 
                      sheet = "factors")
str(factors)
# pego las tablas
data <- merge(replicates, factors, by=c("species","site"))
str(data)

# quitamos las observaciones con NAs
data <- na.omit(data)

# quitamos las variables basadas en area
data <- subset(data, select=-c(Aarea,Narea,Parea)) # ,PNUE,PPUE

# quitamos especies con 1 y 2 observaciones
data$sitexspp <- paste(data$site, data$species)
nobs <- as.data.frame(table(data$sitexspp))
colnames(nobs) <- c("species", "Freq")
nobs <- nobs %>%  filter (Freq > 2)
data <- data[which(data$sitexspp %in% nobs$species),]

# pairs
pairs(data[,3:10])

# eliminamos puntos muy altos
data <- data[-which(data$WUE>11), ]
data <- data[-which(data$Pmass>1), ]
data <- data[-which(data$Amass>1150), ]
# data$height <- log(data$height)
data$height[which(data$height>3500)] <- 3500
data <- data[-which(data$LMA>650), ]
data <- data[-which(data$PPUE>15), ]
data <- data[-which(data$PNUE>480), ]
# reviso como queda
pairs(data[,3:10])



# PCA ######
# como sabemos que neustro espacio funcional es robusto, vamos a usar todas las observaciones

library(factoextra)
library(vegan)

# eliminamos las variables que no nos interesan
mipca.obj <- prcomp(data[,c("LMA","height","Amass","WUE","Nmass","Pmass","PNUE","PPUE")], scale. = TRUE, center = TRUE)
fviz_eig(mipca.obj)
eig <- eigenvals(mipca.obj)
summary(mipca.obj)
fviz_pca_var(mipca.obj, col.var = "contrib")
mipca.obj$rotation[,1:3]

# preparo la tabla que usare a partir de ahora
mipca <- cbind(data, mipca.obj$x[,1:3])
head(mipca)

# descriptivos para paper
sppxorixlf <- unique(mipca[,c("site","species","origin","lifeform")])
table(sppxorixlf$site, sppxorixlf$origin)
table(sppxorixlf$site, sppxorixlf$lifeform)
rm(sppxorixlf)

# boxplot
par(mfrow=c(1,2))
boxplot(mipca$PC1, main="PC1")
boxplot(mipca$PC2, main="PC2")
boxplot(mipca$PC2, main="PC3")

# espacio funcional
autoplot(mipca.obj, data = mipca,
         shape = "origin", size = 1.5, colour = "origin",
         label.color = "black",
         loadings = TRUE, loadings.colour = 'black',
         loadings.label = TRUE, loadings.label.size = 5,
         loadings.label.colour = "black",
         loadings.label.vjust = -0.9) + # loadings.label.hjust = 1.2) +
  # geom_point(aes(shape=origin)) +
  scale_shape_manual(values=c(16,17))+
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.title = element_text(size=14)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))
  

# modelos
mod <- lmer(PC1 ~ 1 + origin * site + (1 | origin/species), data = mipca, REML = FALSE)
summary(mod)
anova(mod)
# r.squaredGLMM(mod)

mod2 <- lmer(PC2 ~ 1 + origin * site  + (1 | origin/species), data = mipca, REML = FALSE)
summary(mod2)
anova(mod2)

mod3 <- lmer(PC3 ~ 1 + origin * site + (1 | origin/species), data = mipca, REML = FALSE)
summary(mod3)
anova(mod3)







