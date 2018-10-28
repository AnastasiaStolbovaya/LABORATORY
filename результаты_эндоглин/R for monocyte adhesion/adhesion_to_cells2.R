
library(stringr)
library(ggplot2)
library(cowplot)
library(lmerTest)
library(multcomp)

path <- "D:/Google Drive/материалы для отчета по гранту/adhesion_THP1_U937/Илья/"
setwd(path)

## ----The Theme-----------------------------------------------------------------------------------

source("D:/Google Drive/материалы для статьи функциональные тесты/Figure_theme.R")

## ----Constants-----------------------------------------------------------------------------------

MAbs <- c("control", "3A7", "2C8", "4C9", "4E4", "5H7")
MAbs_TGF <- c("control", "TGF-beta", "3A7+TGF-beta", "2C8+TGF-beta", "4C9+TGF-beta", "4E4+TGF-beta", "5H7+TGF-beta")

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

read_data <- function(files, type.mab = "mab") {
        df <- read_files(files = files)
        df$mab <- strip_filename(df$file)
        if (type.mab == "mab") {
                df$mab <- factor(df$mab, levels = MAbs)
        } else {
                df$mab <- factor(df$mab, levels = MAbs_TGF)
        }
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

add_first <- function(x) {
        x.last <- x[1] + x[-1]
        new.x <- c(x[1], x.last)
        new.x
}

glm_stat <- function(data, type.mab = "mab") {
        model <- glm(left.circles ~ mab, data, family = "poisson")
        summary.model <- summary(model)
        CI <- confint(model)
        mab <- if(type.mab == "mab") factor(MAbs, levels = MAbs) else factor(MAbs_TGF, levels = MAbs_TGF)
        result <- data.frame(mab = mab,
                             mean = exp(add_first(summary.model$coef[,1])),
                             CI.lower = exp(add_first(CI[,1])),
                             CI.upper = exp(add_first(CI[,2])))
        contr <- summary(glht(model, linfct = mcp(mab = "Tukey")))
        list(result = result, contrasts = contr, model = model)
}


## ----Data----------------------------------------------------------------------------------------

THP.files <- list.files("Short incubation/THP-1/", pattern = "*.csv", full.names = T)
THP_short <- read_data(THP.files)

U937.files <- list.files("Short incubation/U937/", pattern = "*.csv", full.names = T)
U937_short <- read_data(U937.files)

THP.files2 <- list.files("Long incubation/THP-1/", pattern = "*.csv", full.names = T)
THP_long <- read_data(THP.files2, type.mab = "mab.tgf")

U937.files2 <- list.files("Long incubation/U937/", pattern = "*.csv", full.names = T)
U937_long <- read_data(U937.files2, type.mab = "mab.tgf")


## ----Raw data figures (Short incubation)---------------------------------------------------------

ggplot(THP_short, aes(mab, left.circles+removed.circles)) +
        geom_jitter(aes(color = exp), width = 0.15, size = 2.5, alpha = 0.75) +
        theme_bw()


ggplot(U937_short, aes(mab, left.circles+removed.circles)) +
        geom_jitter(aes(color = exp), width = 0.15, size = 2.5, alpha = 0.75) +
        theme_bw()


ggplot(THP_long, aes(mab, left.circles+removed.circles)) +
        geom_jitter(aes(color = exp), width = 0.15, size = 2.5, alpha = 0.75) +
        theme_bw()


ggplot(U937_long, aes(mab, left.circles+removed.circles)) +
        geom_jitter(aes(color = exp), width = 0.15, size = 2.5, alpha = 0.75) +
        theme_bw()

U937_long[U937_long$mab == "control", c("mab", "left.circles")]

U937_long <- U937_long[-c(99, 102),]

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

U937_short_list <- glm_stat(U937_short)
# plot(U937_short_list[[3]])


THP_long_list <- glm_stat(THP_long, type.mab = "mab.tgf")
# plot(THP_long_list[[3]])

#       U937 long incubation

U937_long_list <- glm_stat(U937_long, type.mab = "mab.tgf")
# plot(U937_long_list[[3]])


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

#       U937 short incubation

U937_short_lines <-
text_to_df(
"line  path     y
    1     1   450
    1     2   450
    2     2   600
    2     6   600
    3     2   580
    3     5   580
    4     2   560
    4     4   560
    5     2   540
    5     3   540")

U937_short_stars <-
text_to_df(
"  x    y stars
 1.5  470     *
 2.5  620   ***"
)

fig_U937_short <-
ggplot(U937_short_list[[1]], aes(mab, mean)) +
        geom_bar(stat = "identity", width = 0.45, color = "black", fill = "grey") +
        geom_errorbar(width = 0.25, aes(ymin = CI.lower, ymax = CI.upper)) +
        geom_line(data = U937_short_lines, aes(path, y, group = line)) +
        geom_text(data = U937_short_stars, aes(x, y, label = stars), size = 7) +
        my_theme +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        ggtitle("U937")

fig_U937_short_rus <- fig_U937_short +
        scale_x_discrete(labels = c("К", "IgG1", "2C8", "4C9", "4E4", "5H7")) +
        labs(x = "", y = "Кол-во прикрепившихся клеток")

fig_U937_short_eng <- fig_U937_short +
        scale_x_discrete(labels = c("Control", "IgG1", "2C8", "4C9", "4E4", "5H7")) +
        labs(x = "", y = "Number of cells attached")


#       THP-1 long incubation

THP_long_lines <-
text_to_df(
"line  path     y
    1     1  1400
    1     2  1400
    2     2  1500
    2     3  1500
    3     3  1750
    3     7  1750
    4     3  1700
    4     6  1700
    5     3  1650
    5     4  1650")

THP_long_stars <-
text_to_df(
"  x    y stars
 1.5 1450   ***
 2.5 1550   ***
 3.5 1800   ***"
)

fig_THP_long <- 
ggplot(THP_long_list[[1]], aes(mab, mean)) +
        geom_bar(stat = "identity", width = 0.45, color = "black", fill = "grey") +
        geom_errorbar(width = 0.25, aes(ymin = CI.lower, ymax = CI.upper)) +
        geom_line(data = THP_long_lines, aes(path, y, group = line)) +
        geom_text(data = THP_long_stars, aes(x, y, label = stars), size = 7) +
        my_theme +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1)) +
        ggtitle("THP-1")

fig_THP_long_rus <- fig_THP_long +
        scale_x_discrete(labels = c("К", "+TGF", "IgG1+TGF", "2C8+TGF", "4C9+TGF", "4E4+TGF", "5H7+TGF")) +
        labs(x = "", y = "Кол-во прикрепившихся клеток")

fig_THP_long_eng <- fig_THP_long +
        scale_x_discrete(labels = c("Control", "+TGF", "IgG1+TGF", "2C8+TGF", "4C9+TGF", "4E4+TGF", "5H7+TGF")) +
        labs(x = "", y = "Number of cells attached")

#       U937 long incubation
U937_long_lines <-
text_to_df(
"line  path     y
    1     1   800
    1     2   800
    2     2   750
    2     3   750
    3     3  1050
    3     7  1050
    4     3  1000
    4     6  1000
    5     3   950
    5     5   950
    6     3   900
    6     4   900")

U937_long_stars <-
text_to_df(
"  x    y stars
 1.5  850   ***
 2.5  800    **
 3.5 1100   ***"
)

fig_U937_long <-
ggplot(U937_long_list[[1]], aes(mab, mean)) +
        geom_bar(stat = "identity", width = 0.45, color = "black", fill = "grey") +
        geom_errorbar(width = 0.25, aes(ymin = CI.lower, ymax = CI.upper)) +
        geom_line(data = U937_long_lines, aes(path, y, group = line)) +
        geom_text(data = U937_long_stars, aes(x, y, label = stars), size = 7) +
        my_theme +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1)) +
        ggtitle("U937")

fig_U937_long_rus <- fig_U937_long +
        scale_x_discrete(labels = c("К", "+TGF", "IgG1+TGF", "2C8+TGF", "4C9+TGF", "4E4+TGF", "5H7+TGF")) +
        labs(x = "", y = "Кол-во прикрепившихся клеток")

fig_U937_long_eng <- fig_U937_long +
        scale_x_discrete(labels = c("Control", "+TGF", "IgG1+TGF", "2C8+TGF", "4C9+TGF", "4E4+TGF", "5H7+TGF")) +
        labs(x = "", y = "Number of cells attached")


## ----Final plots---------------------------------------------------------------------------------

fig_rus <-
ggdraw() +
        draw_plot(fig_THP_short_rus, x = 0, y = 0.55, width = 0.5, height = 0.45) +
        draw_plot(fig_U937_short_rus, x = 0.5, y = 0.55, width = 0.5, height = 0.45) +
        draw_plot(fig_THP_long_rus, x = 0, y = 0, width = 0.5, height = 0.55) +
        draw_plot(fig_U937_long_rus, x = 0.5, y = 0, width = 0.5, height = 0.55) +
        draw_plot_label(label = c("А", "Б", "В", "Г"), x = c(0, 0.5, 0, 0.5),
                        y = c(1, 1, 0.55, 0.55), size = 12)

ggsave("adhesion_to_cells-rus.jpeg", fig_rus, width = 18, height = 16, units = "cm")

fig_eng <-
        ggdraw() +
        draw_plot(fig_THP_short_eng, x = 0, y = 0.55, width = 0.5, height = 0.45) +
        draw_plot(fig_U937_short_eng, x = 0.5, y = 0.55, width = 0.5, height = 0.45) +
        draw_plot(fig_THP_long_eng, x = 0, y = 0, width = 0.5, height = 0.55) +
        draw_plot(fig_U937_long_eng, x = 0.5, y = 0, width = 0.5, height = 0.55) +
        draw_plot_label(label = c("А", "B", "C", "D"), x = c(0, 0.5, 0, 0.5),
                        y = c(1, 1, 0.55, 0.55), size = 12)

ggsave("adhesion_to_cells-eng.jpeg", fig_eng, width = 18, height = 16, units = "cm")

