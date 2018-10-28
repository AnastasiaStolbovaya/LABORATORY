## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("flowCore")

library(flowCore)

guide<-read.table("guide.csv", header=T, sep=",", dec=",", stringsAsFactors=FALSE)
str(guide)
is.character(guide$File)

df <- read.FCS.as.data_frame(filename= guide$File[1], x = "SSC.H", y = "FSC.H",
                             log.transform = T, default.sd = 2.5,
                             gate.filter ="default", show = F)

for(i in 2:nrow(guide)){
x <- read.FCS.as.data_frame(filename= guide$File[i], x = "SSC.H", y = "FSC.H",
                             log.transform = T, default.sd = 2.5,
                             gate.filter ="default", show = F)
df<- rbind(df, x)
}

unique(df$file)

head(df)



df_2 <- subset(df, select=c("SSC.H", "FSC.H", "FITC.A", "PE.A", "file"))

unique(guide$File)
guide$File<- factor(guide$File)

df_2$File <- factor(df_2$file)



f.data <- merge(guide, df_2, by="File")

head(f.data)
names(f.data)


#логарифмирование переменной РЕ - интенсивность флуоресценции фикоэритрина
f.data$log_PE <- log10(f.data$PE)
sum(is.na(f.data$log_PE), na.rm = T)




#создаем переменную comb, в каждой ячейке записаны файлы и маркер
f.data$comb <- paste(f.data$File, f.data$Marker)
#делаем выравнивание образцов по переменной comb, одинаковое количество событий и в изотипе и в опыте/контроле. По File выравнивание не происходит так так файлы с изотипами повторяются несколько раз, программа их объединяет и не производит выравнивание. Comb содержит не только файл но и маркер, поэтому программа воспринимает файлы изотипов как разные.
f.data <- equal.samples(f.data, sample.var = "comb")

aggregate(Marker ~ comb, f.data , length)

head(f.data)
# задаем порядок графиков по маркеру при фасетировании
f.data$Marker <- factor(f.data$Marker, levels=c("CD29", "CD62E", "CD54", "CD105"))
#график вот
plot <- ggplot(data=f.data  , aes(x=log_PE, color = treatment, linetype= Isotype))+
  geom_freqpoly(bins = 30)+
  facet_grid(treatment ~ Marker)+
  theme_bw()+
  scale_x_continuous(limits = c(1, 5)) + 
  scale_color_manual(values = c("black","red", "blue"), name="")+
  scale_linetype_manual(values = c(1,2), name="")+
  theme(legend.position = "",
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.text = element_text()) +
  labs(x = expression("Fluorescence intensity," ~log[10]), y = "Count")

ggsave("Plot.jpeg", width = 18, height = 20, units = "cm", dpi = 600)


f.data.var2 <- f.data[f.data$Isotype == "FALSE", ]
unique(f.data.var2$Isotype)

plot_2 <- ggplot(data=f.data.var2  , aes(x=log_PE, color = treatment))+
  geom_freqpoly(bins = 30)+
  facet_grid( ~ Marker)+
  theme_bw()+
  scale_x_continuous(limits = c(1, 5)) + 
  scale_color_manual(values = c("black","red", "blue"), name="")+
  scale_linetype_manual(values = c(1,2), name="")+
  theme(legend.position = "",
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.text = element_text()) +
  labs(x = expression("Fluorescence intensity," ~log[10]), y = "Count")

ggsave("Plot_2.jpeg", width = 18, height = 7, units = "cm", dpi = 600)
