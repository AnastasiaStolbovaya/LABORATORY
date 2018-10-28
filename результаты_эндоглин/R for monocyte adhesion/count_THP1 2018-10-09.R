library(doBy)
library(ggplot2)
library(multcomp)
library(nlme)
library(boot)
library(lme4)
library(emmeans)
getwd()
path <-"D:/Анастасия Столбовая/результаты_эндоглин/R for monocyte adhesion/"
setwd(path)


# # # # # # # # # # # # # 
count_THP1_exp1 <-read.csv(file="count_THP1_exp1.csv", sep=";", dec=",", stringsAsFactors = F)
count_THP1_exp2 <-read.csv(file="count_THP1_exp2.csv", sep=";", dec=",", stringsAsFactors = F)



count_thp1 <- rbind(count_THP1_exp1, count_THP1_exp2)

count_thp1 <- count_thp1[, -c(2:10,12)]
count_thp1$mab<- stringr::str_split_fixed(count_thp1$file, pattern="\\_", n=2)[,2]
count_thp1$mab<- stringr::str_split_fixed(count_thp1$mab, pattern="\\.", n=2)[,1]

count_thp1$well <- stringr::str_split_fixed(count_thp1$name, pattern="\\_", n=2)[,1]
count_thp1$n <- stringr::str_split_fixed(count_thp1$name, pattern="\\_", n=3)[,2]

length(count_thp1$n)
length(count_thp1$well)

count_thp1_sum <- summaryBy(left.circles ~exp+well+mab, data=count_thp1, FUN=mean, keep.names = T)



ggplot(data=count_thp1_sum, aes(x=mab, y=left.circles, color=factor(exp)))+geom_point()

count_thp1_sum$mab <- factor(count_thp1_sum$mab) 
count_thp1_sum$exp <- factor(count_thp1_sum$exp) 
count_thp1_sum$mab<- relevel(count_thp1_sum$mab, ref="3A7") 
count_thp1_sum$left.circles <-round(count_thp1_sum$left.circles)

M1 <- glm(left.circles ~ mab*exp, data=count_thp1_sum, family="poisson")
M2 <- glmer(left.circles ~ mab+(1|exp), data=count_thp1_sum, family="poisson")
M3 <- glm.nb(left.circles ~ mab*exp, data=count_thp1_sum)
M4 <- glm.nb(left.circles ~ mab+exp, data=count_thp1_sum)
M5 <- glmer.nb(left.circles ~ mab+(1|exp), data=count_thp1_sum)
AIC(M1, M2, M3, M4, M5)
summary(M5)

test_M5 <- (glht(model=M5, linfct= mcp(mab = "Tukey")))
summary(test_M5)


M5_means <- emmeans(M5, specs = ~ mab)

M5_means<- as.data.frame(M5_means)


M5_means$emmean_resp <- exp(M5_means$emmean) 

M5_means$asymp.LCL_resp <- exp(M5_means$asymp.LCL) 

M5_means$asymp.UCL_resp <- exp(M5_means$asymp.UCL) 

M5_means$mab <- factor(M5_means$mab, levels=c( "control", "3A7" , "2C8" ,   "4C9" , "4E4" , "5H7"))
levels(M5_means$mab) <- c("К", "IgG1", "2C8" ,   "4C9" , "4E4" , "5H7")


fig4_count_THP1 <- ggplot(M5_means, aes(mab, emmean_resp)) +
        geom_col(color = "black", fill = "grey") +
        geom_errorbar(width = 0.25, aes(ymin = asymp.LCL_resp, ymax =asymp.UCL_resp)) +
        theme_bw()+
        theme(axis.text.x = element_text(angle=30, hjust=1, size=10),
              axis.title.y=element_text(size=12),
              plot.title = element_text(hjust=0.5, face="bold"))+
        xlab("")+
        ylab("Число прикрепленных клеток")+
        ylim(0,2000)+
        ggtitle("THP-1")




# mydata_M5 <- expand.grid(mab = levels(count_thp1_sum$mab))
#  
# mydata_M5$predict <- predict(M5, newdata = mydata_M5, re.form=NA)
# 
# X <- model.matrix(~ mab, data =count_thp1_sum)
#  
# mydata_M5$SE <- sqrt(diag(X %*% vcov(M5) %*% t(X)))








mydata_M2$mab <- factor(mydata_M2$mab, levels=c( "control", "3A7" , "2C8" ,   "4C9" , "4E4" , "5H7"))
levels(mydata_M2$mab) <- c("К", "IgG1", "2C8" ,   "4C9" , "4E4" , "5H7")


fig6 <- ggplot(mydata_M2, aes(mab, predict)) +
        geom_bar(stat = "identity", width = 0.45, color = "black", fill = "grey") +
        geom_errorbar(width = 0.25, aes(ymin = predict-1.96*se, ymax =predict+1.96*se)) +
        # geom_segment(aes(x=c("IgG1"), xend=c("5H7"), y=70, yend=70)) +
        # geom_text(aes(x = "2C8", y = 71, label = "**"), size = 7) +
        theme_bw()+
        # scale_y_continuous(breaks = c("0", "10", "20", "30", "40", "50", "60"))+
        xlab("")+
        ylab("Число клеток THP-1")

ggsave("fig1_count_THP1.jpeg", fig1, width=10, height=12, dpi=600, units="cm")









