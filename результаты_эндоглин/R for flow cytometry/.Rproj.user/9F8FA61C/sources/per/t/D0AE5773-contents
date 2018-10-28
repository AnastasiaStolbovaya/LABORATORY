library(flowCore)

guide_am<-read.table("guide_ad_mol.csv", header=T, sep=";", dec=",", stringsAsFactors=FALSE)
str(guide_am)
is.character(guide_am$file)

 df <- read.FCS.as.data_frame(filename= guide_am$file[1], x = "SSC.H", y = "FSC.H",
                             log.transform = T, default.sd = 2.5,
                             gate.filter ="default", show = F)

for(i in 2:nrow(guide_am)){
  x <- read.FCS.as.data_frame(filename= guide_am$file[i], x = "SSC.H", y = "FSC.H",
                              log.transform = T, default.sd = 2.5,
                              gate.filter ="default", show = F)
  df<- rbind(df, x)
}

unique(df$file)

head(df)


df_2 <- subset(df, select=c("SSC.H", "FSC.H", "FITC.A", "PE.A", "file"))

unique(guide_am$file)
guide_am$file<- factor(guide_am$file)

df_2$file <- factor(df_2$file)

guide_am$maker <- factor(guide_am$maker)
levels(guide_am$maker)


f.data <- merge(guide_am, df_2, by="file")

head(f.data)
names(f.data)
 levels(f.data$maker)

#логарифмирование переменной РЕ - интенсивность флуоресценции фикоэритрина
f.data$log_PE <- log10(f.data$PE)
sum(is.na(f.data$log_PE), na.rm = T)




#создаем переменную comb, в каждой ячейке записаны файлы и маркер
f.data$comb <- paste(f.data$file, f.data$maker)
#делаем выравнивание образцов по переменной comb, одинаковое количество событий и в изотипе и в опыте/контроле. По File выравнивание не происходит так так файлы с изотипами повторяются несколько раз, программа их объединяет и не производит выравнивание. Comb содержит не только файл но и маркер, поэтому программа воспринимает файлы изотипов как разные.
f.data <- equal.samples(f.data, sample.var = "comb")

aggregate(maker ~ comb, f.data , length)

head(f.data)
# задаем порядок графиков по маркеру при фасетировании
# f.data$Marker <- factor(f.data$Marker, levels=c("CD29", "CD62E", "CD54", "CD105"))
#график вот
plot <- ggplot(data=f.data  , aes(x=log_PE, color = mab, linetype= isotype))+
  geom_freqpoly(bins = 30)+
  facet_grid(mab ~ maker)+
  theme_bw()+
  scale_x_continuous(limits = c(1, 5)) + 
  # scale_color_manual(values = c("black","red", "blue"), name="")+
  # scale_linetype_manual(values = c(1,2), name="")+
  theme(legend.position = "",
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.text = element_text()) +
  labs(x = expression("Интенсивность флуоресценции," ~log[10]), y = "Число событий")
 ggsave("Plot_ad_mol.jpeg", width = 20, height = 20, units = "cm", dpi = 600)


#  # #######################

guide_am_2<-read.table("guide_ad_mol_2.csv", header=T, sep=";", dec=",", stringsAsFactors=FALSE)

f.data.2 <- merge(guide_am_2, df_2, by="file")
names(f.data.2)

f.data.2$log_PE <- log10(f.data.2$PE)
sum(is.na(f.data.2$log_PE), na.rm = T)


f.data.2$comb <- paste(f.data.2$file,  f.data.2$maker)

f.data.2 <- equal.samples(f.data.2, sample.var = "comb")

aggregate(mab ~ comb, f.data.2 , length)

f.data.2$maker <- factor(f.data.2$maker, levels=c("CD29","CD44","CD49E","CD54", "CD144","CD146","CD166" ))




plot_2 <- ggplot(data=f.data.2  , aes(x=log_PE, color=control, linetype=isotype_hp))+
  geom_freqpoly(bins = 30)+
  facet_grid(mab ~ maker)+
  theme_bw()+
  scale_x_continuous(limits = c(1, 5)) + 
  scale_color_manual(values = c("red","blue"), name=""+
  scale_linetype_manual(values = c(1,2), name=""))+
  theme(legend.position = "",
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text()) +
  labs(x = expression("Интенсивность флуоресценции," ~log[10]), y = "Число событий")

ggsave("Plot2_ad_mol.jpeg", width = 25, height = 20, units = "cm", dpi = 600)





# ############
guide_am_3<-read.table("guide_ad_mol_3.csv", header=T, sep=";", dec=",", stringsAsFactors=FALSE)


f.data.3 <- merge(guide_am_3, df_2, by="file")

head(f.data.3)
names(f.data.3)


f.data.3$log_PE <- log10(f.data.3$PE)
sum(is.na(f.data.3$log_PE), na.rm = T)


f.data.3$comb <- paste(f.data.3$file, f.data.3$mab)

f.data.3 <- equal.samples(f.data.3, sample.var = "comb")

# aggregate(maker ~ comb, f.data.3 , length)

f.data.3$maker <- factor(f.data.3$maker, levels=c("CD29","CD44","CD49E","CD54", "CD144","C146", "CD166", "HP3A7"))
levels(f.data.3$maker) <- c("CD29","CD44","CD49E","CD54", "CD144","CD146", "CD166", "Ig1")

plot_3 <- ggplot(data=f.data.3  , aes(x=log_PE, color=control))+
  geom_freqpoly(bins = 30)+
  facet_grid(mab ~ maker)+
  theme_bw()+
  scale_x_continuous(limits = c(1, 5)) + 
  # scale_color_manual(values = c("red","blue"), name="")+
  # scale_linetype_manual(values = c(1,2), name=""))+
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text()) +
  labs(x = expression("Интенсивность флуоресценции," ~log[10]), y = "Число событий")

ggsave("Plot2_ad_mol.jpeg", width = 25, height = 20, units = "cm", dpi = 600)



