library(doBy)
library(ggplot2)
library(multcomp)
library(nlme)
library(boot)
library(cowplot)

getwd()
# path <-"D:/Анастасия Столбовая/результаты_эндоглин/R for monocyte adhesion/"



count_THP1_TGF_exp1<- read.csv("count_THP1_TGF_exp1.csv", sep=";", dec=",", stringsAsFactors = F)
count_THP1_TGF_exp2<- read.csv("count_THP1_TGF_exp2.csv", sep=";", dec=",", stringsAsFactors = F)
count_THP1_TGF <-rbind(count_THP1_TGF_exp1, count_THP1_TGF_exp2)

count_THP1_TGF <- count_THP1_TGF[, -c(2:10,12)]

count_THP1_TGF$mab<- stringr::str_split_fixed(count_THP1_TGF$file, pattern="\\_", n=2)[,2]
count_THP1_TGF$mab<- stringr::str_split_fixed(count_THP1_TGF$mab, pattern="\\.", n=2)[,1]
# sum(is.na(count_THP1_TGF$mab))

count_THP1_TGF$well <- stringr::str_split_fixed(count_THP1_TGF$name, pattern="\\_", n=2)[,1]
count_THP1_TGF$n <- stringr::str_split_fixed(count_THP1_TGF$name, pattern="\\_", n=2)[,2]

length(count_THP1_TGF$n)
length(count_THP1_TGF$well)

count_THP1_TGF_sum <- summaryBy(left.circles ~ mab+well+exp+date, data=count_THP1_TGF, FUN=mean, keep.names = T)
# nrow(count_THP1_TGF_sum)

ggplot(data=count_THP1_TGF_sum, aes(x=mab, y=left.circles, color=factor(exp)))+geom_point()

count_THP1_TGF_sum$mab <- factor(count_THP1_TGF_sum$mab) 
count_THP1_TGF_sum$mab <- relevel(count_THP1_TGF_sum$mab, ref="3A7+TGF-beta") 
count_THP1_TGF_sum$exp <- factor(count_THP1_TGF_sum$exp) 

count_THP1_TGF_sum$left.circles <- round(count_THP1_TGF_sum$left.circles, 0)

mod_THP1_TGF <- glmer.nb(left.circles ~ mab+(1|exp), data=count_THP1_TGF_sum)
summary(mod_THP1_TGF)

test_THP1_TGF <- (glht(model=mod_THP1_TGF, linfct= mcp(mab = "Tukey")))
summary(test_THP1_TGF)


mod_THP1_TGF_means <- emmeans(mod_THP1_TGF, specs = ~ mab)

mod_THP1_TGF_means<- as.data.frame(mod_THP1_TGF_means)

mod_THP1_TGF_means$emmean_resp <- exp(mod_THP1_TGF_means$emmean) 

mod_THP1_TGF_means$asymp.LCL_resp <- exp(mod_THP1_TGF_means$asymp.LCL) 

mod_THP1_TGF_means$asymp.UCL_resp <- exp(mod_THP1_TGF_means$asymp.UCL) 

mod_THP1_TGF_means$mab <- factor(mod_THP1_TGF_means$mab, levels=c( "control", "TGF-beta", "3A7+TGF-beta" ,  "2C8+TGF-beta" ,   "4C9+TGF-beta" , "4E4+TGF-beta" , "5H7+TGF-beta"))

levels(mod_THP1_TGF_means$mab) <- c( "К", "TGF-beta", "Ig1+TGF-beta" , "2C8+TGF-beta" ,   "4C9+TGF-beta" , "4E4+TGF-beta" , "5H7+TGF-beta")


fig3_count_THP1_TGF<- ggplot(mod_THP1_TGF_means, aes(x = mab, y = emmean_resp)) + 
        geom_col(fill = "gray", color = "black") + 
        geom_errorbar(aes(ymin = asymp.LCL_resp, ymax = asymp.UCL_resp), width = 0.2)+
        theme_bw()  +
        theme(axis.text.x = element_text(angle = 30, hjust=1, size=10),
              axis.title.y=element_text(size=12),
              plot.title = element_text(hjust=0.5, face="bold"))+
        ylab("Число прикрепленных клеток")+
        xlab("")+
        ggtitle("THP-1")
        


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



mod_1 <- glm(left.circles ~ mab*exp, data=count_THP1_TGF_sum, family="poisson")
mod_2 <- glm.nb(left.circles ~ mab*exp, data=count_THP1_TGF_sum)
mod_3 <- glm.nb(left.circles ~ mab+exp, data=count_THP1_TGF_sum)
mod_4 <- glm.nb(left.circles ~ mab, data=count_THP1_TGF_sum)

summary(mod_4)
AIC(mod_2, mod_3, mod_4)

test <- (glht(model=mod, linfct= mcp(mab = "Tukey")))
summary(test)

test_Dunnet <- (glht(model=mod_4, linfct= mcp(mab = "Dunnet")))
summary(test_Dunnet)


mydata_ <- data.frame(mab = levels(count_THP1_TGF_sum$mab))

mydata_$predict <- predict(mod_4, newdata = mydata_, re.form = NA, type = "response", se.fit=T)$fit

mydata_$se <- predict(mod_4, newdata = mydata_, re.form = NA, type = "response", se.fit=T)$se


mydata_$mab <- factor(mydata_$mab, levels=c( "control", "TGF-beta", "3A7+TGF-beta" ,  "2C8+TGF-beta" ,   "4C9+TGF-beta" , "4E4+TGF-beta" , "5H7+TGF-beta"))
levels(mydata_$mab) <- c( "К", "TGF-beta", "Ig1+TGF-beta" , "2C8+TGF-beta" ,   "4C9+TGF-beta" , "4E4+TGF-beta" , "5H7+TGF-beta")


fig3a <- ggplot(mydata_, aes(mab, predict)) +
        geom_bar(stat = "identity", width = 0.45, color = "black", fill = "grey") +
        geom_errorbar(width = 0.25, aes(ymin = predict-1.96*se, ymax =predict+1.96*se)) +
        geom_segment(aes(x=c("К"), xend=c("Ig1+TGF-beta"), y=2000, yend=2000)) +
        geom_text(aes(x = "TGF-beta", y = 2050, label = "***"), size = 7) +
        theme_bw()+
        xlab("")+
        ylab("Число клеток THP-1")+
        theme(axis.text.x = element_text(angle=30, hjust=1))
       
       

ggsave("fig3a_count_THP1_TGF.jpeg", fig3a, width=10, height=12, dpi=600, units="cm")






fig3 <- ggplot(mydata_, aes(mab, predict)) +
        geom_bar(stat = "identity", width = 0.45, color = "black", fill = "grey") +
        geom_errorbar(width = 0.25, aes(ymin = predict-1.96*se, ymax =predict+1.96*se)) +
        geom_segment(aes(x=c("TGFb"), xend=c("4E4+TGFb"), y=164, yend=164)) +
        geom_segment(aes(x=c("Ig1+TGFb"), xend=c("2C8+TGFb"), y=155, yend=155)) +
        geom_segment(aes(x=c("Ig1+TGFb"), xend=c("4C9+TGFb"), y=158, yend=158)) +
        geom_text(aes(x = "Ig1+TGFb", y = 165, label = "*"), size = 7) +
        geom_text(aes(x = "2C8+TGFb", y = 159, label = "*"), size = 7) +
        geom_segment(aes(x=c("К"), xend=c("TGFb"), y=169, yend=169)) +
        geom_segment(aes(x=c("К"), xend=c("Ig1+TGFb"), y=171, yend=171)) +
        geom_segment(aes(x=c("К"), xend=c("2C8+TGFb"), y=173, yend=173))+
        geom_segment(aes(x=c("К"), xend=c("4C9+TGFb"), y=175, yend=175)) +
        geom_segment(aes(x=c("К"), xend=c("4E4+TGFb"), y=177, yend=177))+
        geom_segment(aes(x=c("К"), xend=c("5H7+TGFb"), y=179, yend=179))+
        geom_text(aes(x = "TGFb", y = 181, label = "***"), size = 7) +
        theme_bw()+
        # scale_y_continuous(breaks = c("0", "20", "40", "60", "80","100", "120", "140", "160","180", "200"))+
        xlab("")+
        ylab("Число клеток THP-1")+
        theme(axis.text.x = element_text(angle=30, hjust=1))



ggsave("fig3_count_THP1_TGF.jpeg", fig3, width=10, height=12, dpi=600, units="cm")
