
library(flowCore)
#library(flowViz)
library(ggplot2)
library(grid)

read.FCS.as.data_frame <- function(filename, x = "SSC.A", y = "FSC.A",
                                   log.transform = FALSE , gate.filter = "no", default.sd = 2,
                                   show = FALSE) {
        raw.data <- read.FCS(filename, alter.names = TRUE)
        if(log.transform == TRUE){
                logTrans <- logTransform(transformationId="log10-transformation", logbase=10, r=1, d=1)
                trans <- transformList(c(x, y), logTrans)
                raw.data <- transform(raw.data, trans)
        }
        if(is.character(gate.filter)){
                if(gate.filter == "default"){
                        morphGate <- norm2Filter(x, y, filterId = "default", scale = default.sd)
                        sub.data <- Subset(raw.data, morphGate)
                } else if(gate.filter == "no"){
                        sub.data <- raw.data
                }
        } else  sub.data <- Subset(raw.data, gate.filter)
        raw.df <- data.frame(exprs(raw.data), subsetting = "Original")
        subset.df  <- data.frame(exprs(sub.data), subsetting = "Subset")
        df.common <- rbind(raw.df, subset.df)
        df.common$subsetting <- factor(df.common$subsetting, levels = c("Original", "Subset"))
        if(show == TRUE){
                eval(parse(text = paste0("df.common$X <- df.common$", x)))
                eval(parse(text = paste0("df.common$Y <- df.common$", y)))
                pict <-
                ggplot(df.common, aes(X, Y)) +
                        geom_point(alpha = 0.3) +
                        facet_wrap(~subsetting) +
                        theme_bw() +
                        theme(strip.background = element_blank(),
                              panel.border = element_rect(color = "black"),
                              panel.grid.major = element_line(size = 0.5),
                              panel.grid.minor = element_blank(),
                              strip.text = element_text(size = 14),
                              axis.title = element_text(size = 14),
                              axis.text = element_text(size = 12)) +
                        labs(x = x, y = y)
                print(pict)
        }
        subset.df$file <- filename
        return(subset(subset.df, select = -subsetting))
}


combine2df <- function(df1, df2){
        df1.to.df2 <- match(names(df1), names(df2))
        col.to.add1 <- names(df1)[which(is.na(df1.to.df2))]
        if(length(col.to.add1) > 0) df2[col.to.add1] <- NA
        df2.to.df1 <- match(names(df2), names(df1))
        col.to.add2 <- names(df2)[which(is.na(df2.to.df1))]
        if(length(col.to.add2) > 0) df1[col.to.add2] <- NA
        return(rbind(df1, df2))
}

equal.samples <- function(data, sample.var, size = "min"){
        fact <- as.factor(subset(data, select = sample.var)[,1])
        if(size == "min") size <- min(table(fact))
        rows <- numeric(0)
        for(f in levels(fact)){
                rows <- c(rows, sample(which(fact == f), size = size, replace = F))
        }
        return(data[rows,])
}

put.mark <-  function(x, y, units = "npc", length, angle = 45, text = "", vjust = "centre", hjust = "centre" ){
        radians <- angle*pi/180
        x1 <- length * cos(radians) + x
        y1 <- length * sin(radians) + y
        line <- grid.lines(x = unit(c(x, x1), units),
                           y = unit(c(y, y1), units))
        mark <- grid.text(label = text,
                          x = unit(x1, units), y = unit(y1, units),
                          just = c(hjust, vjust), rot = 0,
                          gp = gpar(fontsize = 12))
        
        return(line)
        return(mark)
}
