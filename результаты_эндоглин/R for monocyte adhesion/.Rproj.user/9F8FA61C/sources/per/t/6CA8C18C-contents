library(doBy)
library(ggplot2)
library(multcomp)
library(nlme)
library(cowplot)

getwd()
path <-"D:/Анастасия Столбовая/результаты_эндоглин/R for monocyte adhesion/"
setwd(path)




count_U937_exp1 <- read.csv("count_U937_exp1.csv", sep=";", dec=",", stringsAsFactors = F)
count_U937_exp2 <- read.csv("count_U937_exp2.csv", sep=";", dec=",", stringsAsFactors = F)

count_U937 <- rbind(count_U937_exp1, count_U937_exp2)


count_U937 <- count_U937[, -c(2:10,12)]
count_U937$mab<- stringr::str_split_fixed(count_U937$file, pattern="\\_", n=2)[,2]
count_U937$mab<- stringr::str_split_fixed(count_U937$mab, pattern="\\.", n=2)[,1]
count_U937$well <- stringr::str_split_fixed(count_U937$name, pattern="\\_", n=2)[,1]
count_U937$n <- stringr::str_split_fixed(count_U937$name, pattern="\\_", n=2)[,2]


length(count_U937$n)
length(count_U937$well)

count_U937_sum <- summaryBy(left.circles ~ mab+well+exp+date, data=count_U937, FUN=mean, keep.names = T)

count_U937_sum$left.circles <- round(count_U937_sum$left.circles, 0)

ggplot(data=count_U937_sum, aes(x=mab, y=left.circles, color=factor(exp)))+
        geom_point()

count_U937_sum$mab <- factor(count_U937_sum$mab) 
count_U937_sum$mab <- relevel(count_U937_sum$mab, ref="3A7") 
count_U937_sum$exp <- factor(count_U937_sum$exp) 


# mod_1 <- glm(left.circles ~ mab, data=, family="poisson")
mod_U937 <- glmer.nb(left.circles ~ mab+(1|exp), data=count_U937_sum)

summary(mod_U937)

test_mod_U937 <- (glht(model=mod_U937, linfct= mcp(mab = "Tukey")))
summary(test_mod_U937)

mod_U937_means <- emmeans(mod_U937, specs = ~ mab)

mod_U937_means<- as.data.frame(mod_U937_means)

mod_U937_means$emmean_resp <- exp(mod_U937_means$emmean) 

mod_U937_means$asymp.LCL_resp <- exp(mod_U937_means$asymp.LCL) 

mod_U937_means$asymp.UCL_resp <- exp(mod_U937_means$asymp.UCL) 

mod_U937_means$mab <- factor(mod_U937_means$mab, levels=c( "control",  "3A7" ,  "2C8" ,   "4C9" , "4E4" , "5H7"))

levels(mod_U937_means$mab) <- c( "К",  "Ig1" , "2C8" ,   "4C9" , "4E4" , "5H7")


fig2_count_U937_sum <- ggplot(mod_U937_means, aes(x = mab, y = emmean_resp)) + 
        geom_col(fill = "gray", color = "black") + 
        geom_errorbar(aes(ymin = asymp.LCL_resp, ymax = asymp.UCL_resp), width = 0.2)+
        theme_bw()  +
        theme(axis.text.x = element_text(angle = 30, hjust=1))+
        ylab("Число клеток U-937")+
        xlab("")+
        geom_segment(aes(x=c("Ig1"), xend=c("4C9"), y=1200, yend=1200))+
        geom_segment(aes(x=c("Ig1"), xend=c("4E4"), y=1210, yend=1210))+
        geom_segment(aes(x=c("Ig1"), xend=c("5H7"), y=1220, yend=1220))+
        geom_text(aes(x = "2C8", y = 1230, label = "***"), size = 7)



# mydata_count_2 <- data.frame(mab = levels(count_2$mab))
# 
# mydata_count_2$predict <- predict(mod, newdata = mydata_count_2, re.form = NA, type = "response", se.fit=T)$fit
# 
# mydata_count_2$se <- predict(mod, newdata = mydata_count_2, re.form = NA, type = "response", se.fit=T)$se
# 
# mydata_count_2$mab <- factor(mydata_count_2$mab, levels=c( "control", "HP3A7" , "2C8" ,   "4C9" , "4E4" , "5H7"))
# levels(mydata_count_2$mab) <- c("К", "IgG1", "2C8" ,   "4C9" , "4E4" , "5H7")
# 
# fig2 <- ggplot(mydata_count_2, aes(mab, predict)) +
#         geom_bar(stat = "identity", width = 0.45, color = "black", fill = "grey") +
#         geom_errorbar(width = 0.25, aes(ymin = predict-1.96*se, ymax =predict+1.96*se)) +
#         theme_bw()+
#         ylim(0,70)+
#         xlab("")+
#         ylab("Число клеток U937")
# 
# ggsave("fig2_count_U937.jpeg", fig2, width=10, height=12, dpi=600, units="cm")
# 
# 
# 
# 
# fig3<- gridExtra::grid.arrange(fig1, fig2, nrow=1)
# 
# ggsave("fig3.jpeg", fig3, width=15, height=12, dpi=600, units="cm")

