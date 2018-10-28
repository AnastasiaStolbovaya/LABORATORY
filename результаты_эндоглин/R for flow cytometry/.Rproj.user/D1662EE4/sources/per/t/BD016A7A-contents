library(ggplot2)
library("cowplot")

path <- "D:/temp/"
setwd(path)

source("D:/temp/flow.functions2.R")


guide_L <- read.table("guide_L_2.csv", sep = ",", dec=",", header = T, stringsAsFactors = F)
################
# reading

f.data <- data.frame(SSC.H = NA,
                     FSC.H = NA,
                     FITC = NA,
                     PE = NA,
                     file = NA)

for(i in 1:nrow(guide_L)){
        df <- read.FCS.as.data_frame(guide_L$File[i], x = "SSC.H", y = "FSC.H",
                                     log.transform = T, default.sd = 2.5,
                                     gate.filter = guide_L$filter[i], show = F)
        if(i <= 56) {
                df2 <- subset(df, select = c("SSC.H", "FSC.H", "FITC.A", "PE.A", "file"))
                names(df2)[c(3, 4)] <- c("FITC", "PE")
        } 
        else {
                df2 <- subset(df, select = c("SSC.H", "FSC.H", "FL1.H", "FL2.H", "file"))
                names(df2)[c(3, 4)] <- c("FITC", "PE")
        }
        f.data <- rbind(f.data, df2)
}


sum(f.data$PE <0, na.rm = T)
sum(is.na(f.data$PE), na.rm = T)

head(f.data)
str(f.data)

# удаление первой строки с пустыми значениями таблицы f.data
f.data <- f.data[complete.cases(f.data),]
head(f.data)

#при помощи функции merge соединяем прочитанные файлы fcs с таблицой guide по общему столбцу File_F
unique(guide_L$File)
guide_L$File_F <- factor(guide_L$File)

f.data$File_F <- factor(f.data$file)

#merge(guide_L, f.data)

f.data_2 <- merge(guide_L, f.data, by="File_F")
head(f.data_2)
names(f.data_2)

#логарифмирование переменной РЕ - интенсивность флуоресценции фикоэритрина
f.data_2$log_PE <- log10(f.data_2$PE)
sum(is.na(f.data_2$log_PE), na.rm = T)

############# строим графики
f.data_hepg2 <- f.data_2[f.data_2$Cell %in% "Hep_G2", ]
#графики 4 маркера Hep G2 через 48 часов после удаления препарата
f.data_hepg2_48h <- f.data_hepg2[f.data_hepg2$Time %in% c("48"), ]
f.data_hepg2_48h <- f.data_hepg2_48h[f.data_hepg2_48h$Marker %in% c("CD29", "CD44", "CD54", "CD166"), ]
unique(f.data_hepg2_48h$Marker)
#создаем переменную comb, в каждой ячейке записаны файлы и маркер
f.data_hepg2_48h$comb <- paste(f.data_hepg2_48h$File_F, f.data_hepg2_48h$Marker)
#делаем выравнивание образцов по переменной comb, одинаковое количество событий и в изотипе и в опыте/контроле. По File_F выравнивание не происходит так так файлы с изотипами повторяются несколько раз, программа их объединяет и не производит выравнивание. Comb содержит не только файл но и маркер, поэтому программа воспринимает файлы изотипов как разные.
f.data_hepg2_48h <- equal.samples(f.data_hepg2_48h, sample.var = "comb")

aggregate(Marker ~ comb, f.data_hepg2_48h , length)
# задаем порядок графиков по маркеру при фасетировании
f.data_hepg2_48h$Marker <- factor(f.data_hepg2_48h$Marker, levels=c("CD29", "CD44", "CD54", "CD166"))
#график вот
Hep_48h <- ggplot(data=f.data_hepg2_48h  , aes(x=log_PE, color = L_concentration, linetype= Isotype))+
        geom_freqpoly(bins = 30)+
        facet_grid(L_concentration ~ Marker)+
        theme_bw()+
        scale_x_continuous(limits = c(1, 5)) + 
        scale_color_manual(values = c("red", "black"), name="")+
        scale_linetype_manual(values = c(1,2), name="")+
        theme(legend.position = "",
              strip.background = element_blank(),
              strip.text = element_text(size = 12),
              axis.title = element_text(size = 10),
              axis.text = element_text(size = 8),
              legend.text = element_text()) +
        labs(x = expression("Интенсивность флуоресценции," ~log[10]), y = "Число событий")

ggsave("Hep_48h.jpeg", width = 16, height = 10, units = "cm", dpi = 600)

#######
f.data_ea <- f.data_2[f.data_2$Cell %in% "EA.hy926", ]
f.data_ea_48h <- f.data_ea[f.data_ea$Time %in% c("48"), ]
f.data_ea_48h <- f.data_ea_48h[f.data_ea_48h$Marker %in% c("CD29", "CD44", "CD105", "CD166"), ]

unique(f.data_ea_48h$Time)

f.data_ea_48h$comb <- paste(f.data_ea_48h$File_F, f.data_ea_48h$Marker)

f.data_ea_48h <- equal.samples(f.data_ea_48h, sample.var = "comb")
aggregate(Marker ~ comb, f.data_ea_48h , length)
# задаем порядок графиков по маркеру при фасетировании
f.data_ea_48h$Marker <- factor(f.data_ea_48h$Marker, levels=c("CD29", "CD44", "CD105", "CD166"))
#график вот
EA_48h <- ggplot(data=f.data_ea_48h  , aes(x=log_PE, color = L_concentration, linetype= Isotype))+
        geom_freqpoly(bins = 30)+
        facet_grid(L_concentration ~ Marker)+
        theme_bw()+
        scale_x_continuous(limits = c(1, 5)) + 
        scale_color_manual(values = c("red", "black"), name="")+
        scale_linetype_manual(values = c(1,2), name="")+
        theme(legend.position = "",
              strip.background = element_blank(),
              strip.text = element_text(size = 12),
              axis.title = element_text(size = 10),
              axis.text = element_text(size = 8),
              legend.text = element_text()) +
        labs(x = expression("Интенсивность флуоресценции," ~log[10]), y = "Число событий")

ggsave("EA_48h.jpeg", width = 16, height = 10, units = "cm", dpi = 600)

#######
f.data_ecv <- f.data_2[f.data_2$Cell %in% "ECV304", ]
f.data_ecv_48h <- f.data_ecv[f.data_ecv$Time %in% c("24"), ]
f.data_ecv_48h <- f.data_ecv_48h[f.data_ecv_48h$Marker %in% c("CD29", "CD44", "CD105", "CD166", "CD62E"), ]

unique(f.data_ecv_48h$Marker)

f.data_ecv_48h$comb <- paste(f.data_ecv_48h$File_F, f.data_ecv_48h$Marker)

f.data_ecv_48h <- equal.samples(f.data_ecv_48h, sample.var = "comb")
aggregate(Marker ~ comb, f.data_ecv_48h , length)
# задаем порядок графиков по маркеру при фасетировании
f.data_ecv_48h$Marker <- factor(f.data_ecv_48h$Marker, levels=c("CD29", "CD44", "CD105", "CD166","CD62E"))
#график вот
ECV_24h <- ggplot(data=f.data_ecv_48h  , aes(x=log_PE, color = L_concentration, linetype= Isotype))+
        geom_freqpoly(bins = 30)+
        facet_grid(L_concentration ~ Marker)+
        theme_bw()+
        scale_x_continuous(limits = c(0, 5)) + 
        scale_color_manual(values = c("red", "black"), name="")+
        scale_linetype_manual(values = c(1,2), name="")+
        theme(legend.position = "",
              strip.background = element_blank(),
              strip.text = element_text(size = 12),
              axis.title = element_text(size = 10),
              axis.text = element_text(size = 8),
              legend.text = element_text()) +
        labs(x = expression("Интенсивность флуоресценции," ~log[10]), y = "Число событий")

ggsave("ECV_48h.jpeg", width = 16, height = 10, units = "cm", dpi = 600)

#соединяем три объекта
aligned_plots <- align_plots(Hep_48h, EA_48h, ECV_24h, plotlist = NULL, align = 'v')
Cells_48h <- plot_grid(Hep_48h, EA_48h, ECV_24h, ncol=1, labels=c('Hep G2','EA.hy926', 'ECV304'), label_size = 10, scale=1)

save_plot("phen_cells_48h.jpg", Cells_48h, ncol = 1, base_height = 7, base_aspect_ratio = 1.1, base_width = NULL)




#########

f.data_ea <- f.data_2[f.data_2$Cell %in% "EA.hy926", ]
f.data_ea_0h <- f.data_ea[f.data_ea$Time %in% c("0"), ]
f.data_ea_0h <- f.data_ea_0h[f.data_ea_0h$Marker %in% c("CD29", "CD44", "CD105", "CD166"), ]

unique(f.data_ea_0h$Time)

f.data_ea_0h$comb <- paste(f.data_ea_0h$File_F, f.data_ea_0h$Marker)

f.data_ea_0h <- equal.samples(f.data_ea_0h, sample.var = "comb")
aggregate(Marker ~ comb, f.data_ea_0h , length)

f.data_ea_0h$Marker <- factor(f.data_ea_0h$Marker, levels=c("CD29", "CD44", "CD105", "CD166"))
#график вот
EA_0h <- ggplot(data=f.data_ea_0h  , aes(x=log_PE, color = L_concentration, linetype= Isotype))+
        geom_freqpoly(bins = 30)+
        facet_grid(L_concentration ~ Marker)+
        theme_bw()+
        scale_x_continuous(limits = c(1, 5)) + 
        scale_color_manual(values = c("red", "black"), name="")+
        scale_linetype_manual(values = c(1,2), name="")+
        theme(legend.position = "",
              strip.background = element_blank(),
              strip.text = element_text(size = 12),
              axis.title = element_text(size = 10),
              axis.text = element_text(size = 8),
              legend.text = element_text()) +
        labs(x = expression("Интенсивность флуоресценции," ~log[10]), y = "Число событий")

ggsave("EA_0h.jpeg", width = 16, height = 10, units = "cm", dpi = 600)


#######
f.data_hepg2 <- f.data_2[f.data_2$Cell %in% "Hep_G2", ]

f.data_hepg2_0h <- f.data_hepg2[f.data_hepg2$Time %in% c("0"), ]
unique(f.data_hepg2_0h$Marker)

f.data_hepg2_0h$comb <- paste(f.data_hepg2_0h$File_F, f.data_hepg2_0h$Marker)

f.data_hepg2_0h <- equal.samples(f.data_hepg2_0h, sample.var = "comb")

aggregate(Marker ~ comb, f.data_hepg2_0h , length)

#график вот
Hep_0h <- ggplot(data=f.data_hepg2_0h  , aes(x=log_PE, color = L_concentration, linetype= Isotype))+
        geom_freqpoly(bins = 30)+
        facet_grid(L_concentration ~ Marker)+
        theme_bw()+
        scale_x_continuous(limits = c(1, 5)) + 
        scale_color_manual(values = c("red", "black"), name="")+
        scale_linetype_manual(values = c(1,2), name="")+
        theme(legend.position = "",
              strip.background = element_blank(),
              strip.text = element_text(size = 12),
              axis.title = element_text(size = 10),
              axis.text = element_text(size = 8),
              legend.text = element_text()) +
        labs(x = expression("Интенсивность флуоресценции,"
                            ~log[10]), y = "Число событий")

ggsave("Hep_0h.jpeg", width = 7, height = 10, units = "cm", dpi = 600)


#соединяем два объекта

aligned_plots <- align_plots(Hep_0h, EA_0h, plotlist = NULL, align = 'hv')
Cells_0h <- plot_grid(Hep_0h, EA_0h, nrow=1, ncol=2, labels=c('Hep G2', 'EA.hy926'), align="hv", label_size = 12, scale=1, rel_widths=c(1,2))

save_plot("phen_cells_0h.jpg", Cells_0h, ncol = 2, base_height =7, base_aspect_ratio = 1.1, base_width = NULL)














