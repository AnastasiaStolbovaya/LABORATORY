
# Распакуй архив на диск D в папку temp

library(ggplot2)

path <- "D:/temp/"
setwd(path)

source("D:/Dropbox/!Scripts/flow.functions2.R")

# файл с описанием экспериментов
guide <- read.csv("file.guide_culture.conditions3a.csv", sep = ";", dec = ",", stringsAsFactors = F)

# демонстрация работы функции и гейтирования:

read.FCS.as.data_frame(guide$file[1], x = "SSC.H", y = "FSC.H", gate.filter = "default", default.sd = 2.5, show = T)

# на выходе имеем дата-фрейм с ГЕЙТИРОВАННЫМИ данными


flow.data <-
        do.call(rbind,
                lapply(guide$file, read.FCS.as.data_frame,
                       x = "SSC.H", y = "FSC.H", log.transform = T, default.sd = 2.5,
                       gate.filter = "default", # морфологическое гейтирование по умолчанию
                       show = F)
        )

flow.data2 <- merge(guide, flow.data)             # объединение двух дата-фреймов в один по именам файлов
flow.data2$log_FITC.A <- log10(flow.data2$FITC.A) # лог-трансформацию лучше делать до выравнивания образцов, т.к. будут значения NaN
head(flow.data2)
flow.data3 <- flow.data2[!is.nan(flow.data2$log_FITC.A),] # Убираю лишнее
flow.data3$comb <- paste(flow.data3$file.name, flow.data3$string) # файлы с изотипами читаются несклько раз, поэтмоу в этом случае выравниваю образцы 
                                                                  # не по файлу, а по созданной уникальной переменной comb
head(flow.data3)

flow.data4 <- equal.samples(flow.data3, sample.var = "comb")      # Выравниевание образцов
flow.data4$cell <- factor(flow.data4$cell, levels = c("EA.hy926", "MSC-6", "MSC-16"))
levels(flow.data4$cell)[2:3] <- c("МСК-6", "МСК-16")
head(flow.data4)
# Ввожу русские подписи

flow.data4$regimen <- factor(flow.data4$regimen, levels = c("no change", "change", "re-seeding", "isotype"))
levels(flow.data4$regimen) <- c("Без смены среды", "С заменой среды", "С пересевом клеток", "Изотипический контроль")

names(flow.data4)

figure <-
        ggplot(flow.data4, aes(log_FITC.A, group = file)) +
        geom_freqpoly(binwidth = 0.15, aes(linetype = regimen, color = regimen)) +
        scale_linetype_manual(values = c(1, 1, 1, 2), name = "") +
        scale_color_manual(values = c("red", "green4", "blue", "black"), name = "") +
        scale_x_continuous(limits = c(0.5, 4.5)) +
        facet_grid(Mab.tested ~ cell) +
        theme_bw() +
        theme(legend.position = "bottom",
              strip.background = element_blank(),
              panel.border = element_rect(color = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.text = element_text(size = 14),
              axis.title = element_text(size = 14),
              axis.text = element_text(size = 12),
              legend.text = element_text(size = 12)) +
        labs(x = expression("Интенсивность флуоресценции," ~log[10]), y = "Количество событий")


fig <- figure + guides(linetype = guide_legend(nrow = 2)) # коррекция вида легенды

fig

ggsave("culture.conditions3a.jpeg", fig, width = 18, height = 13, units = "cm", dpi = 600)


