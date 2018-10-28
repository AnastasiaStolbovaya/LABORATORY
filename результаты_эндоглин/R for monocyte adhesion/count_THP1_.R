library(doBy)
library(ggplot2)
library(multcomp)
library(nlme)
library(boot)

getwd()
path <-"D:/Анастасия Столбовая/результаты_эндоглин/R for monocyte adhesion/"






count_THP1 <-read.csv(file="count_THP1.csv", sep=";", dec=",", stringsAsFactors = F)
count_THP1 <- count_THP1[, -c(2:10,12)]
count_THP1$mab <- "1"
count_THP1[1:18, 4] <- "2C8"
count_THP1[19:36, 4] <- "HP3A7"
count_THP1[37:54, 4] <- "4C9"
count_THP1[55:72, 4] <- "4E4"
count_THP1[73:90, 4] <- "5H7"
count_THP1[91:108, 4] <- "control"


count_THP1$well <- stringr::str_split_fixed(count_THP1$name, pattern="\\_", n=2)[,1]
count_THP1$n <- stringr::str_split_fixed(count_THP1$name, pattern="\\_", n=3)[,2]

length(count_THP1$n)
length(count_THP1$well)

 count <- summaryBy(left.circles ~ mab+well, data=count_THP1, FUN=mean, keep.names = T)

 
 ggplot(data=count, aes(x=mab, y=left.circles))+geom_point()

count$mab <- factor(count$mab) 
 
M1 <- glm.nb(left.circles ~ mab, data=count)
summary(M1)

test <- (glht(model=M1, linfct= mcp(mab = "Tukey")))
summary(test)



mydata_count <- data.frame(mab = levels(count$mab))
 
mydata_count$predict <- predict(M1, newdata = mydata_count, re.form = NA, type = "response", se.fit=T)$fit

mydata_count$se <- predict(M1, newdata = mydata_count, re.form = NA, type = "response", se.fit=T)$se


# confint(M1, parm=, level=0.95)

mydata_count$mab <- factor(mydata_count$mab, levels=c( "control", "HP3A7" , "2C8" ,   "4C9" , "4E4" , "5H7"))
levels(mydata_count$mab) <- c("К", "IgG1", "2C8" ,   "4C9" , "4E4" , "5H7")


fig1 <- ggplot(mydata_count, aes(mab, predict)) +
        geom_bar(stat = "identity", width = 0.45, color = "black", fill = "grey") +
        geom_errorbar(width = 0.25, aes(ymin = predict-1.96*se, ymax =predict+1.96*se)) +
        geom_segment(aes(x=c("IgG1"), xend=c("5H7"), y=70, yend=70)) +
        geom_text(aes(x = "2C8", y = 71, label = "**"), size = 7) +
        theme_bw()+
        # scale_y_continuous(breaks = c("0", "10", "20", "30", "40", "50", "60"))+
        xlab("")+
        ylab("Число клеток THP-1")

ggsave("fig1_count_THP1.jpeg", fig1, width=10, height=12, dpi=600, units="cm")









