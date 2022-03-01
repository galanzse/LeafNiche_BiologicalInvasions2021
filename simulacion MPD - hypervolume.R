
library(readxl)
library(hypervolume)
library(tidyverse)

data <- read_excel("C:/Users/Javier/OneDrive/TESIS Y PUBLICACIONES/TESIS/2. Leaf-Niche Funk/Original Data/database_JG.xlsx") %>% as.data.frame()
colnames(data)
traits <- c("LMA","Height","N_mass","P_mass","A_area","A_mass","WUE","N_area","PNUE","P_area","PPUE")

# correlaciones
pairs(data[,traits])

data_reduced <- data %>% select("Species","LMA","Height","N_mass") %>% na.omit()
data_reduced$LMA <- scale(data_reduced$LMA, center=T, scale=T)
data_reduced$Height <- scale(data_reduced$Height, center=T, scale=T)
data_reduced$N_mass <- scale(data_reduced$N_mass, center=T, scale=T)

# hypervolume - MPD
nruns <- 100
myresults <- matrix(ncol=3, nrow=nruns)
colnames(myresults) <- c("hypervolume","mpd","outlier")

for (l in 1:nruns) {
  mysamp <- data_reduced[sample(rownames(data_reduced), 25, replace=F), c("LMA", "Height", "N_mass")]
  
  hyp <- hypervolume(mysamp, method="box")
  myresults[l,"hypervolume"] <- hyp@Volume
  
  mydist <- mysamp %>% dist()
  myresults[l,"mpd"] <- mydist %>% mean()
  myout <- mydist %>% as.vector() %>% sort(decreasing=T)
  myresults[l,"outlier"] <- mean(myout[1:10])
  
  print(paste(rep("-",l)))
}

myresults <- myresults %>% as.data.frame()

myresults$hypervolume <- myresults$hypervolume + min(myresults$hypervolume)
myresults$hypervolume <- scale(myresults$hypervolume, center=F, scale=T)

myresults$mpd <- myresults$mpd + min(myresults$mpd)
myresults$mpd <- scale(myresults$mpd, center=F, scale=T)

hist(myresults$outlier, 20)
myresults$outlier_category <- NA
for (l in 1:nruns) {
  if (myresults$outlier[l] < 5) { myresults$outlier_category[l] <- "close" }
  if (myresults$outlier[l] > 10) { myresults$outlier_category[l] <- "far" }
}
myresults$outlier_category[is.na(myresults$outlier_category)] <- "medium"

ggplot(data=myresults, aes(y=hypervolume, x=mpd, colour=outlier_category)) +
  geom_point() +
  geom_smooth(method = "lm")

# hypervolume - richness
nruns <- 1000
myresults <- matrix(ncol=2, nrow=nruns)
colnames(myresults) <- c("hypervolume","richness")
myresults[,"richness"] <- c(rep(20,200), rep(30,200), rep(40,200), rep(50,200), rep(60,200)) 

for (l in 1:nruns) {
  mysamp <- data_reduced[sample(rownames(data_reduced), myresults[l,"richness"], replace=F), c("LMA", "Height", "N_mass")]
  
  hyp <- hypervolume(mysamp, method="box")
  myresults[l,"hypervolume"] <- hyp@Volume

  print(paste(rep("-",l)))
}

myresults <- myresults %>% as.data.frame()

boxplot(myresults$hypervolume ~ myresults$richness)

# bandwidth
nruns <- 500
myresults <- matrix(ncol=5, nrow=nruns)
colnames(myresults) <- c("sd","bandwith","hypervolume","mpd","outlier")
for (l in 1:nruns) {
  mytrait <- rnorm(25, mean=10, sd=runif(1, min = 1, max = 10))
  # mytrait <- runif(25, min = 0, max = 25)
  myresults[l,"sd"] <- sd(mytrait)
  
  myresults[l,"bandwith"] <- estimate_bandwidth(as.data.frame(mytrait), method="silverman")
  
  hyp <- hypervolume(as.data.frame(mytrait), method="box")
  myresults[l,"hypervolume"] <- hyp@Volume
  
  mydist <- mytrait %>% dist()
  myresults[l,"mpd"] <- mydist %>% mean()
  myresults[l,"outlier"] <- max(mydist)
  
  print(paste(rep("-",l)))
}

ggplot(data=as.data.frame(myresults), aes(y=hypervolume, x=mpd, colour=outlier)) +
  geom_point() +
  geom_smooth(method="lm")

pairs(myresults)

# hipervolumenes sin rarefaccion
source("pca.R")

hyp_total <- matrix(ncol = 2, nrow = length(unique(mipca$site)))
rownames(hyp_total) <- unique(mipca$site)
colnames(hyp_total) <- c("native","invasive")

mpd_total <- matrix(ncol = 4, nrow = length(unique(mipca$site)))
rownames(mpd_total) <- unique(mipca$site)
colnames(mpd_total) <- c("native","invasive","out_nat","out_inv")


for (s in unique(mipca$site)) {
  
  nat <- mipca %>% subset(site==s & origin=="native") %>% select(PC1,PC2,PC3)
  inv <- mipca %>% subset(site==s & origin=="invasive") %>% select(PC1,PC2,PC3)

  hyp_total[s,"native"] <- hypervolume(nat, method="box")@Volume
  hyp_total[s,"invasive"] <- hypervolume(inv, method="box")@Volume
  
  mpd_total[s,"native"] <- nat %>% dist() %>% mean()
  mpd_total[s,"out_nat"] <- nat %>% dist() %>% max()
  mpd_total[s,"invasive"] <- inv %>% dist() %>% mean()
  mpd_total[s,"out_inv"] <- inv %>% dist() %>% max()

}

plot(hyp_total[,"invasive"]-hyp_total[,"native"], mpd_total[,"invasive"]-mpd_total[,"native"])
abline(h=0,v=0)

hyp_total <- hyp_total %>% as.data.frame()
hyp_total$site <- rownames(hyp_total)
rownames(hyp_total) <- NULL
hyp_total <- gather(hyp_total, "origin", "volume", 1:2)

mpd_total <- mpd_total %>% as.data.frame()
mpd_total$site <- rownames(mpd_total)
rownames(mpd_total) <- NULL
mpd_total <- gather(mpd_total, "origin", "mpd", 1:2)
mpd_total$outlier <- c(mpd_total$out_nat[1:8], mpd_total$out_inv[9:16])
mpd_total$out_nat <- NULL
mpd_total$out_inv <- NULL

total_div <- merge(hyp_total, mpd_total, by=c("site","origin"))

plot(total_div$volume, total_div$outlier)
plot(total_div$mpd, total_div$outlier)

ggplot(data=total_div, aes(y=volume, x=mpd, colour=origin)) +
  geom_point() +
  geom_smooth(method="lm")

      