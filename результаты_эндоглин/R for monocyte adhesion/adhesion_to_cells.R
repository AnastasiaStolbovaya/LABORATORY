
library(stringr)
library(ggplot2)
library(lmerTest)

path <- "D:/Google Drive/материалы для отчета по гранту/adhesion_THP1_U937/Илья/"
setwd(path)

## ----The Theme-----------------------------------------------------------------------------------

source("D:/Google Drive/материалы для статьи функциональные тесты/Figure_theme.R")

## ----Constants-----------------------------------------------------------------------------------

MAbs <- c("control", "3A7", "2C8", "4C9", "4E4", "5H7")

## ----Fucntions-----------------------------------------------------------------------------------

read_file <- function(filename) {
        df <- read.csv(file = filename, sep = ";", dec = ",")
        df$exp <- basename(filename)
        df
}

read_files <- function(path = getwd(), files) {
        do.call(rbind, lapply(file.path(path, files), read_file))
}

strip_filename <- function(filename) {
        stripped <- str_remove(str_remove(filename, "report_"), ".csv")
        stripped
}

read_data <- function(files) {
        df <- read_files(files = files)
        df$mab <- strip_filename(df$file)
        df$mab <- factor(df$mab, levels = MAbs)
        df
}

contrast_readable <- function(contrasts, round) {
        df <- data.frame(Contrasts = contrasts$contrasts$Contrast,
                         p.value = round(contrasts$contrasts$p.value, round))
        df
}

text_to_df <- function(text){
        read.table(textConnection(text), header = T)
}

## ----Data----------------------------------------------------------------------------------------

THP.files <- list.files("Short incubation/THP-1/", pattern = "*.csv", full.names = T)
THP_short <- read_data(THP.files)

U937.files <- list.files("Short incubation/U937/", pattern = "*.csv", full.names = T)
U937_short <- read_data(U937.files)

## ----Raw data figures (Short incubation)---------------------------------------------------------

ggplot(THP_short, aes(mab, left.circles)) +
        geom_jitter(aes(color = exp), width = 0.15, size = 2.5, alpha = 0.75) +
        theme_bw()


ggplot(U937_short, aes(mab, left.circles)) +
        geom_jitter(aes(color = exp), width = 0.15, size = 2.5, alpha = 0.75) +
        theme_bw()

## ----Statistics----------------------------------------------------------------------------------

#       THP-1 short incubation

THP_short.model <- glmer(left.circles ~ mab + (1|exp), data = THP_short, family = "poisson")

summary(THP_short.model)
plot(THP_short.model)

THP_short_contr <- psycho::get_contrasts(THP_short.model, formula = "mab", adjust = "bonferroni")
contrast_readable(THP_short_contr, 5)

THP_short_result <- 
data.frame(
        mab = THP_short_contr$means$mab,
        mean = exp(THP_short_contr$means$Mean),
        CI.lower = exp(THP_short_contr$means$CI_lower),
        CI.upper = exp(THP_short_contr$means$CI_higher))

THP_short_result$mab <- factor(THP_short_result$mab, levels = MAbs)


#       U937 short incubation

U937_short.model <- glm(left.circles ~ mab, U937_short, family = "poisson")

summary(U937_short.model)

contrast::contrast(U937_short.model, list("mab"))

summary(multcomp::glht(U937_short.model, mcp(rank="Tukey")))

## ----Figures-------------------------------------------------------------------------------------

#       THP-1 short incubation

THP_short_lines <-
text_to_df(
"line  path     y
    1     1  1800
    1     2  1800
    2     2  2000
    2     6  2000
    3     2  1950
    3     5  1950
    4     2  1900
    4     4  1900")

THP_short_stars <-
text_to_df(
"  x    y stars
 1.5 1850   ***
 3.0 2050   ***"
)

fig_THP_short <-
ggplot(THP_short_result, aes(mab, mean)) +
        geom_bar(stat = "identity", width = 0.45, color = "black", fill = "grey") +
        geom_errorbar(width = 0.25, aes(ymin = CI.lower, ymax = CI.upper)) +
        geom_line(data = THP_short_lines, aes(path, y, group = line)) +
        geom_text(data = THP_short_stars, aes(x, y, label = stars), size = 7) +
        my_theme +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        ggtitle("THP-1")

fig_THP_short_rus <- fig_THP_short +
        scale_x_discrete(labels = c("К", "IgG1", "2C8", "4C9", "4E4", "5H7")) +
        labs(x = "", y = "Кол-во прикрепившихся клеток")

fig_THP_short_eng <- fig_THP_short +
        scale_x_discrete(labels = c("Control", "IgG1", "2C8", "4C9", "4E4", "5H7")) +
        labs(x = "", y = "Number of cells attached")



