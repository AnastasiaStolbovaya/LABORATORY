#Read files-reports after image analisys2

getwd()
path <-"E:/Лаборатория гибридомной технологии/результаты_эндоглин/R for monocyte adhesion/2018-10-25_co-culture_U937_TGF"
setwd(path)
getwd()

read_my_csv <- function(f, sep){
        data <- read.csv(file=f, sep=sep, dec=",", stringsAsFactors = F)
        data$file <- f
        data
}

files_cells <- list.files("E:/Лаборатория гибридомной технологии/результаты_эндоглин/R for monocyte adhesion/2018-10-25_co-culture_U937_TGF", pattern = "*.csv", full.names = F)

 f1 <- lapply(files_cells, read_my_csv, sep=";")
 class(f1)
 lapply(f1, class)
 table <- do.call(plyr::rbind.fill, f1)
 class(table)
 table$date <- "2018-10-25"
 table$exp <- 3
 write.table(table, file="count_U937_TGF_exp3.csv", sep=";", row.names = F)
 