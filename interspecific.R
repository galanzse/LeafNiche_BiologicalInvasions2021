

inter_ori <- dB.origin.total %>% subset(group %in% c("interspecific")) %>% spread("group", "variance")
colnames(inter_ori)[1] <- ""

inter_life <- dB.life.total %>% subset(group %in% c("interspecific")) %>% spread("group", "variance")
colnames(inter_life)[1] <- ""

rbind(inter_ori, inter_life) %>% write.table("interspecific.txt")




