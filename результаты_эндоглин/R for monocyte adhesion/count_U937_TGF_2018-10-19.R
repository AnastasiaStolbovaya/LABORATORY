library(doBy)
library(ggplot2)
library(multcomp)
library(nlme)
library(boot)
library(cowplot)
library(lme4)
library(MASS)
library(emmeans)
getwd()
# path <-"D:/Анастасия Столбовая/результаты_эндоглин/R for monocyte adhesion/"
# setwd(path)
count_U937_TGF_exp1 <- read.csv("count_U937_TGF_exp1.csv", sep=";", dec=",", stringsAsFactors = F)
count_U937_TGF_exp2 <- read.csv("count_U937_TGF_exp2.csv", sep=";", dec=",", stringsAsFactors = F)

count_U937_TGF <- rbind(count_U937_TGF_exp1, count_U937_TGF_exp2)


count_U937_TGF <- count_U937_TGF[, -c(2:10,12)]
count_U937_TGF$mab<- stringr::str_split_fixed(count_U937_TGF$file, pattern="\\_", n=2)[,2]
count_U937_TGF$mab<- stringr::str_split_fixed(count_U937_TGF$mab, pattern="\\.", n=2)[,1]
count_U937_TGF$well <- stringr::str_split_fixed(count_U937_TGF$name, pattern="\\_", n=2)[,1]
count_U937_TGF$n <- stringr::str_split_fixed(count_U937_TGF$name, pattern="\\_", n=2)[,2]

length(count_U937_TGF$n)
length(count_U937_TGF$well)

count_U937_TGF_sum <- summaryBy(left.circles ~ mab+well+exp+date, data=count_U937_TGF, FUN=mean, keep.names = T)

count_U937_TGF_sum$left.circles <- round(count_U937_TGF_sum$left.circles, 0)

ggplot(data=count_U937_TGF_sum, aes(x=mab, y=left.circles, color=factor(exp)))+geom_point()


ggplot(data=count_U937_TGF_sum, aes(x=mab, y=left.circles, color=factor(exp)))+geom_boxplot()




count_U937_TGF_sum$mab <- factor(count_U937_TGF_sum$mab) 
count_U937_TGF_sum$mab <- relevel(count_U937_TGF_sum$mab, ref="3A7+TGF-beta") 
count_U937_TGF_sum$exp <- factor(count_U937_TGF_sum$exp) 



mod_3 <- glm(left.circles ~ mab, data=count_U937_TGF_sum, family="poisson")

###################################
mod <- glm.nb(left.circles ~  exp / mab, data=count_U937_TGF_sum)

summary(mod)

drop1(mod, test = "Chi")

plot(mod)


mydata <- expand.grid(mab = unique(count_U937_TGF_sum$mab), exp = unique((count_U937_TGF_sum$exp)))


mydata$Predict <- predict(mod, newdata = mydata, type = "response", se.fit = T)$fit

mydata$SE <- predict(mod, newdata = mydata, type = "response", se.fit = T)$se.fit


mydata$mab <- factor(mydata$mab, levels=c( "control", "TGF-beta", "3A7+TGF-beta" ,  "2C8+TGF-beta" ,   "4C9+TGF-beta" , "4E4+TGF-beta" , "5H7+TGF-beta"))
levels(mydata$mab) <- c( "К", "TGF-beta", "Ig1+TGF-beta" , "2C8+TGF-beta" ,   "4C9+TGF-beta" , "4E4+TGF-beta" , "5H7+TGF-beta")






ggplot(mydata, aes(x = mab, y = Predict)) + 
        geom_col(fill = "gray", color = "black") + 
        facet_wrap(~exp) + theme_bw() +
        geom_errorbar(aes(ymin = Predict - 1.96*SE, ymax = Predict + 1.96*SE), width = 0.2) + 
        theme(axis.text.x = element_text(angle = 30, hjust=1))+
        xlab("")+
        ylab("Число клеток U-937")


# Удобный пакет, который корректно считает CI, но все совпадает с тем, что делаем через predict. Множитель везде 1.96

m_means <- emmeans(mod, specs = ~ mab)
m_means <- as.data.frame(m_means)
# Можно сразу всовывать в ggplot, но не забывать делать обратное преобразование.Для negnb - это exp()



# Пробую упрятать exp в случайный фактор

mod_lme <- glmer.nb(left.circles ~ mab + (1|exp), data=count_U937_TGF_sum)

summary(mod_lme)

m_means <- emmeans(mod_lme, specs = ~ mab)

m_means <- as.data.frame(m_means)

m_means$emmean_resp <- exp(m_means$emmean) 

m_means$asymp.LCL_resp <- exp(m_means$asymp.LCL) 

m_means$asymp.UCL_resp <- exp(m_means$asymp.UCL) 

m_means$mab <- factor(m_means$mab, levels=c( "control", "TGF-beta", "3A7+TGF-beta" ,  "2C8+TGF-beta" ,   "4C9+TGF-beta" , "4E4+TGF-beta" , "5H7+TGF-beta"))

levels(m_means$mab) <- c( "К", "TGF-beta", "Ig1+TGF-beta" , "2C8+TGF-beta" ,   "4C9+TGF-beta" , "4E4+TGF-beta" , "5H7+TGF-beta")



fig1_count_U937_TGF_sum <- ggplot(m_means, aes(x = mab, y = emmean_resp)) + 
        geom_col(fill = "gray", color = "black") + 
        geom_errorbar(aes(ymin = asymp.LCL_resp, ymax = asymp.UCL_resp), width = 0.2)+
        theme_bw()  +
        theme(axis.text.x = element_text(angle = 30, hjust=1, size=10),
              axis.title.y = element_text(size=12),
              plot.title = element_text(hjust=0.5, face="bold"))+
        ylab("Число прикрепленных клеток")+
        xlab("")+
        geom_segment(aes(x=c("Ig1+TGF-beta"), xend=c("4C9+TGF-beta"), y=1850, yend=1850))+
        geom_segment(aes(x=c("Ig1+TGF-beta"), xend=c("4E4+TGF-beta"), y=1900, yend=1900))+
        geom_segment(aes(x=c("Ig1+TGF-beta"), xend=c("5H7+TGF-beta"), y=1960, yend=1960))+
        geom_text(aes(x = "2C8+TGF-beta", y = 1990, label = "***"), size = 7)+
        ylim(0,2000)+
        ggtitle("U-937")
       

#########################







mod_1 <- glm.nb(left.circles ~ mab*exp, data=count_U937_TGF_sum)
mod_2 <- glmer.nb(left.circles ~ mab+(1|exp), data=count_U937_TGF_sum)
summary(mod_2)
AIC(mod,mod_1 )

drop1(mod_1, test="Chi")


test <- (glht(model=mod_2, linfct= mcp(mab = "Tukey")))
summary(test)

# test_Dunnet <- (glht(model=mod, linfct= mcp(mab = "Dunnet")))
# summary(test_Dunnet)


mydata_ <- data.frame(mab = levels(count_U937_TGF_sum$mab))

mydata_$predict <- predict(mod_2, newdata = mydata_, re.form = NA, type="response")

# mydata_$se <- predict(mod, newdata = mydata_, re.form = NA, type = "response", se.fit=T)$se

X <- model.matrix(~ mab, data =count_U937_TGF_sum)

mydata_$se <- sqrt(diag(vcov(mod_2)))
# mydata_$SE <- sqrt(diag(X %*% vcov(mod_2) %*% t(X)))


mydata_$mab <- factor(mydata_$mab, levels=c( "control", "TGF-beta", "3A7+TGF-beta" ,  "2C8+TGF-beta" ,   "4C9+TGF-beta" , "4E4+TGF-beta" , "5H7+TGF-beta"))
levels(mydata_$mab) <- c( "К", "TGF-beta", "Ig1+TGF-beta" , "2C8+TGF-beta" ,   "4C9+TGF-beta" , "4E4+TGF-beta" , "5H7+TGF-beta")

is.numeric(mydata_$predict)

fig4 <- ggplot(mydata_, aes(mab, predict)) +
        geom_bar(stat = "identity", width = 0.45, color = "black", fill = "grey") +
        # geom_errorbar(width = 0.25, aes(ymin = predict-1.96*se, ymax =predict+1.96*se)) +
        # scale_y_continuous(breaks = c("0", "20", "40", "60", "80","100", "120", "140", "160","180", "200"))+
        xlab("")+
        ylab("Число клеток U-937")+
        theme_bw()+
        theme(axis.text.x = element_text(angle=30, hjust=1))+
        ylim(0,150)



ggsave("fig4_count_U937_TGF.jpeg", fig4, width=10, height=12, dpi=600, units="cm")


fig5<- gridExtra::grid.arrange(fig3, fig4, nrow=1, heights=2:1)

ggsave("fig5_THP1_U937_TGF.jpeg", fig5, width=15, height=12, dpi=600, units="cm")
