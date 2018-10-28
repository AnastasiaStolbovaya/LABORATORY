library(doBy)
library(ggplot2)
library(multcomp)
library(nlme)
library(cowplot)

getwd()
path <-"D:/Анастасия Столбовая/результаты_эндоглин/R for monocyte adhesion/"



# read_my_csv <- function(f, sep){
#         data <- read.csv(file=f, sep=sep, dec=",", stringsAsFactors = F)
#         data$file <- f
#         data
# }
# 
# files_cells_2 <- list.files("D:/Анастасия Столбовая/результаты_эндоглин/R for monocyte adhesion/2018-09-20_co-culture_U937", pattern = "*.csv", full.names = F)
#  path <- "D:/Анастасия Столбовая/результаты_эндоглин/R for monocyte adhesion/2018-09-20_co-culture_U937"
#  setwd(path)
# 
# f1 <- lapply(files_cells_2, read_my_csv, sep=";")
# class(f1)
# lapply(f1, class)
# count_U937 <- do.call(plyr::rbind.fill, f1)
# class(count_U937)
# write.table(count_U937, file="count_U937.csv", sep=";", row.names = F)


count_U937 <-read.csv(file="count_U937.csv", sep=";", dec=",", stringsAsFactors = F)

count_U937 <- count_U937[, -c(2:10,12)]
count_U937$mab <- "1"
count_U937[1:18, 4] <- "2C8"
count_U937[19:36, 4] <- "HP3A7"
count_U937[37:54, 4] <- "4C9"
count_U937[55:72, 4] <- "4E4"
count_U937[73:90, 4] <- "5H7"
count_U937[91:108, 4] <- "control"


count_U937$well <- stringr::str_split_fixed(count_U937$name, pattern="\\_", n=2)[,1]
count_U937$n <- stringr::str_split_fixed(count_U937$name, pattern="\\_", n=3)[,2]

length(count_U937$n)
length(count_U937$well)

count_2 <- summaryBy(left.circles ~ mab+well, data=count_U937, FUN=mean, keep.names = T)


ggplot(data=count_2, aes(x=mab, y=left.circles))+geom_point()

count_2$mab <- factor(count_2$mab) 

mod_1 <- glm(left.circles ~ mab, data=count_2, family="poisson")
mod <- glm.nb(left.circles ~ mab, data=count_2)
summary(mod)

test <- (glht(model=mod, linfct= mcp(mab = "Tukey")))
summary(test)



mydata_count_2 <- data.frame(mab = levels(count_2$mab))

mydata_count_2$predict <- predict(mod, newdata = mydata_count_2, re.form = NA, type = "response", se.fit=T)$fit

mydata_count_2$se <- predict(mod, newdata = mydata_count_2, re.form = NA, type = "response", se.fit=T)$se

mydata_count_2$mab <- factor(mydata_count_2$mab, levels=c( "control", "HP3A7" , "2C8" ,   "4C9" , "4E4" , "5H7"))
levels(mydata_count_2$mab) <- c("К", "IgG1", "2C8" ,   "4C9" , "4E4" , "5H7")

fig2 <- ggplot(mydata_count_2, aes(mab, predict)) +
        geom_bar(stat = "identity", width = 0.45, color = "black", fill = "grey") +
        geom_errorbar(width = 0.25, aes(ymin = predict-1.96*se, ymax =predict+1.96*se)) +
        theme_bw()+
        ylim(0,70)+
        xlab("")+
        ylab("Число клеток U937")

ggsave("fig2_count_U937.jpeg", fig2, width=10, height=12, dpi=600, units="cm")




fig3<- gridExtra::grid.arrange(fig1, fig2, nrow=1)

ggsave("fig3.jpeg", fig3, width=15, height=12, dpi=600, units="cm")

