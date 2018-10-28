
library(EBImage)

filenamer <- function(path = getwd(), base.name, ext, sep = "_", prefix = "", suffix = "",
                      replace = FALSE) {
        short.name <- paste0(ifelse(missing(prefix), "", paste0(prefix, sep)),
                             base.name,
                             ifelse(missing(suffix), "", paste0(sep, suffix)))
        full.name <- paste(short.name, ext, sep = ".")
        if(!missing(path)) {
                file.name <- file.path(path, full.name)
                if(!replace){
                        n <- 0
                        while(file.exists(file.name)){
                                n <- n + 1
                                file.name <- file.path(path,
                                                       paste0(short.name, " (", n, ")", ".", ext))
                        }
                }
                return(file.name)
        } else return(full.name)
}

operator_fun <- function(a, b, operator){
        if(operator == "+") return(a + b)
        if(operator == "-") return(a - b)
        if(operator == "*") return(a * b)
        if(operator == "/") return(a / b)
        if(operator == "^") return(a ^ b)
        else return(NA)
}

br_contrast <- function(img, coef, FUN){
        if(length(coef) != length(FUN)) stop("coef and FUN have different lengths")
        if(length(FUN) > 1) {
                img_2 <- br_contrast(img, coef = coef[length(coef)], FUN = FUN[length(FUN)] )
                return(br_contrast(img_2, coef = coef[1:(length(coef) - 1)], FUN = FUN[1:(length(FUN) - 1)] ))
        } else {
                img_corrected <- operator_fun(img, coef, operator = FUN)
                img_corrected <- ifelse(img_corrected < 0, 0, img_corrected)
                img_corrected <- ifelse(img_corrected > 1, 1, img_corrected)
                return(img_corrected)
        }
}


color_extraction <- function(img, channel = c("R", "G", "B")){
        w <- dim(img)[1]
        h <- dim(img)[2]
        col.val <- c(R = 1, G = 2, B = 3)
        grey.img <- img[1:w, 1:h, col.val[channel]]
        colorMode(grey.img) <- Grayscale
        return(grey.img)
}


binarization <- function(img, otsu.offset = 0.05){
        thr <- otsu(img)
        img.binary <- img > thr + otsu.offset
        return(img.binary)
}


detect_Objects <- function(img.binary, dist.metric = "euclidean", water.tolerance = 1, water.ext = 1){
        img.water <- watershed(distmap(img.binary, metric = dist.metric), tolerance = water.tolerance, ext = water.ext)
        return(img.water)
}

object_filter <- function(objects, size.threshold){
        obj.table <- table(objects)
        objects.left = rmObjects(objects, index = which(obj.table < size.threshold), reenumerate = TRUE)
        objects.removed = rmObjects(objects, index = which(obj.table >= size.threshold), reenumerate = TRUE)
        filtered <- list(obj.table = obj.table,
                         obj.left = objects.left,
                         obj.removed = objects.removed,
                         N.left = sum(obj.table >= size.threshold) - 1,
                         N.removed = sum(obj.table < size.threshold))
        return(filtered)
}

overlay <- function(colorless.ori, bin.img, channel, col.cor = 0.25){
        switch(channel,
               "R" = rgbImage(red = bin.img, green = colorless.ori - col.cor, blue = colorless.ori - col.cor),
               "G" = rgbImage(red = colorless.ori - col.cor, green = bin.img, blue = colorless.ori - col.cor),
               "B" = rgbImage(red = colorless.ori - col.cor, green = colorless.ori - col.cor, blue = bin.img))
}



#-------------------------------------------------------

circle_counter_core <- function(image, return.images = 1:8,
                                channel,
                                br.con.coef, br.con.FUN,
                                blur.int, otsu.offset,
                                dist.metric = "euclidean",
                                water.tolerance = 1, water.ext = 1,
                                size.threshold = 2000){
        # создание индекса img.report
        if (missing(return.images)) index <- numeric(0) else index <- return.images
        img.report <- list()
        
        img.report[[1]] <- if (1 %in% index) image else NA
        cat("step 1/9 \n")
        
        colorless <- color_extraction(image, channel = channel)
        pixels <- length(colorless)                             # Общий размер фотографии 
        if (!8 %in% index) rm(image)                            # Освобождаю память, если 8 не понадобится
        img.report[[2]] <- if (2 %in% index) colorless else NA
        cat("step 2/9 \n")
        
        br.contrasted <- br_contrast(colorless, coef = br.con.coef, FUN = br.con.FUN)
        if (!6 %in% index) rm(colorless)
        img.report[[3]] <- if (3 %in% index) br.contrasted else NA
        cat("step 3/9 \n")
        
        blured <- gblur(br.contrasted, sigma = blur.int)
        rm(br.contrasted)
        img.report[[4]] <- if (4 %in% index) blured else NA
        cat("step 4/9 \n")
        
        bin <- binarization(blured, otsu.offset = otsu.offset)
        rm(blured)
        img.report[[5]] <- if (5 %in% index) bin else NA
        cat("step 5/9 \n")
        
        obj <- detect_Objects(bin, dist.metric = dist.metric, water.tolerance = water.tolerance,
                              water.ext = water.ext)
        
        colored.circles <- colorLabels(obj)
        img.report[[6]] <- if (6 %in% index) colored.circles else NA
        cat("step 6/9 \n")
        
        filtered.obj <- object_filter(obj, size.threshold = size.threshold)
        img.report[[7]] <- if (7 %in% index) filtered.obj$obj.left else NA
        cat("step 7/9 \n")
        
        img.report[[8]] <- if (8 %in% index) filtered.obj$obj.removed else NA
        cat("step 8/9 \n")
        
        obj <- detect_Objects(bin, dist.metric = dist.metric,
                              water.tolerance = water.tolerance,
                              water.ext = water.ext)
        over <- overlay(colorless.ori = colorless,
                        bin.img = filtered.obj$obj.left,
                        channel = channel)
        img.report[[9]] <- if (7 %in% index) over else NA
        cat("step 9/9 \n")
        
        num.report <- data.frame(channel = channel,
                                 circles = filtered.obj$N.left,
                                 br.cont = paste(br.con.FUN, br.con.coef, collapse = " "),
                                 blur.int = blur.int,
                                 otsu.offset = otsu.offset,
                                 dist.metric = dist.metric,
                                 water.tolerance = water.tolerance,
                                 water.ext = water.ext,
                                 size.threshold = size.threshold,
                                 left.circles = filtered.obj$N.left,
                                 removed.circles = filtered.obj$N.removed)
        return(list(img.report = img.report, num.report = num.report))
}

circle_counter <- function(image.folder = getwd(), files,
                           output.folder = getwd(),
                           return.images, channel,
                           br.con.coef, br.con.FUN, blur.int,
                           otsu.offset = 0,
                           dist.metric = "euclidean",
                           water.tolerance = 1, water.ext = 1,
                           size.threshold = 2000) {
        if (!dir.exists(output.folder)) dir.create(output.folder)
        circles <- data.frame()
        for (f in files) {
                cat(f, "\n")
                img <- readImage(file.path(image.folder, f))
                result <- circle_counter_core(img, return.images, channel,
                                               br.con.coef, br.con.FUN, blur.int,
                                               otsu.offset, dist.metric, water.tolerance,
                                               water.ext, size.threshold)
                rm(img)
                bn <- tools::file_path_sans_ext(basename(f))
                if (!all(is.na(result$img.report))) { # запись изображений на диск
                        cat("Saving images...")
                        saving.folder <- file.path(output.folder, bn)
                        if (!dir.exists(saving.folder)) dir.create(saving.folder)
                        image_saver(result$img.report, "circles", saving.folder, base.name = bn)
                        cat("Done!\n\n")
                }
                new.circles.line <- cbind(data.frame(name = bn), result$num.report)
                circles <- plyr::rbind.fill(circles, new.circles.line)
                write.table(x = circles, file = filenamer(output.folder, base.name = "report", ext = "csv", replace = TRUE),
                            sep = ";", dec = ",", row.names = FALSE)
                cat("\n")
                print(result$num.report[,c("left.circles", "removed.circles")])
                cat("\n")
                rm(result)
        }
        return(circles)
}


scratch_area_core <- function(image, return.images = 1:8,
                              channel,
                              br.con.coef, br.con.FUN,
                              blur.int, otsu.offset,
                              brush.size, brush.shape,
                              largest) {
        # создание индекса img.report
        if (missing(return.images)) index <- numeric(0) else index <- return.images
        
        img.report <- list()
        
        img.report[[1]] <- if (1 %in% index) image else NA
        cat("step 1/8 \n")
        
        colorless <- color_extraction(image, channel = channel)
        pixels <- length(colorless)                             # Общий размер фотографии 
        if (!8 %in% index) rm(image)                            # Освобождаю память, если 8 не понадобится
        img.report[[2]] <- if (2 %in% index) colorless else NA
        cat("step 2/8 \n")
        
        br.contrasted <- br_contrast(colorless, coef = br.con.coef, FUN = br.con.FUN)
        rm(colorless)
        img.report[[3]] <- if (3 %in% index) br.contrasted else NA
        cat("step 3/8 \n")
        
        blured <- gblur(br.contrasted, sigma = blur.int)
        rm(br.contrasted)
        img.report[[4]] <- if (4 %in% index) blured else NA
        cat("step 4/8 \n")
        
        bin <- binarization(blured, otsu.offset = otsu.offset)
        rm(blured)
        img.report[[5]] <- if (5 %in% index) bin else NA
        cat("step 5/8 \n")
        
        kern <- makeBrush(brush.size, shape = brush.shape)
        open.img <- opening(bin, kern)
        rm(bin)
        img.report[[6]] <- if (6 %in% index) open.img else NA
        cat("step 6/8 \n")
        
        bwl <- bwlabel(open.img)
        rm(open.img)
        bwtab <- table(bwl)
        no.bg <- bwtab[-1]                              # Первый член вектора - фон, его удаляю
        ord.bwtab <- order(no.bg, decreasing = TRUE)
        scratch.obj <- names(no.bg) %in% ord.bwtab[1:largest]
        remove.obj <- which(!scratch.obj)
        scratch <- rmObjects(bwl, remove.obj)
        img.report[[7]] <- if (7 %in% index) scratch else NA # scratch
        cat("step 7/8 \n")
        
        if (8 %in% index) {
                over <- rgbImage(red   = ifelse(scratch > 0, 1, image[,,1]), # При largest = 2, по некоторым пикселям может быть значение 2
                                 green = ifelse(scratch > 0, 0, image[,,2]),
                                 blue  = ifelse(scratch > 0, 0, image[,,3]))
                img.report[[8]] <- over
                cat("step 8/8 \n")
        } else {
                img.report[[8]] <- NA
                cat("step 8/8 skipped \n")
        }
        
        
        num.report <- data.frame(channel = channel,
                                 br.cont = paste(br.con.FUN, br.con.coef, collapse = " "),
                                 blur.int = blur.int,
                                 otsu.offset = otsu.offset,
                                 brush.size = brush.size,
                                 brush.shape = brush.shape,
                                 largest = largest,
                                 scr.area = sum(no.bg[scratch.obj]),
                                 perc = round(sum(no.bg[scratch.obj]) * 100 / pixels, 1))
        return(list(img.report = img.report, num.report = num.report))
}


image_saver <- function(img.report, FUN, output.dir = getwd(), base.name) {
        dictionary <- switch(FUN,
                             "scratch_area" = c("original", "channel", "brightness_contrast", "blured",
                                                "binary", "open", "scratch", "overlaid"),
                             "circles" = c("original", "channel", "brightness_contrast", "blured",
                                           "binary", "colored_circles", "objects_left", "objects_removed", "overlaid"))
        for(i in 1:length(img.report)) {
                if (is.na(img.report[i])) next
                writeImage(img.report[[i]], filenamer(path = output.dir, base.name = base.name,
                                                      prefix = paste0(0, i), suffix = dictionary[i],
                                                      ext = "jpeg", replace = FALSE))
        }
}


scratch_area <- function(image.folder = getwd(), files,
                         output.folder = getwd(),
                         return.images, channel,
                         br.con.coef, br.con.FUN, blur.int,
                         otsu.offset, brush.size, brush.shape, largest = 1) {
        if (!dir.exists(output.folder)) dir.create(output.folder)
        areas <- data.frame()
        for (f in files) {
                cat(f, "\n")
                img <- readImage(file.path(image.folder, f))
                result <- scratch_area_core(img, return.images, channel,
                                            br.con.coef, br.con.FUN, blur.int,
                                            otsu.offset, brush.size, brush.shape, largest)
                rm(img)
                bn <- tools::file_path_sans_ext(basename(f))
                if (!all(is.na(result$img.report))) { # запись изображений на диск
                        cat("Saving images...")
                        saving.folder <- file.path(output.folder, bn)
                        if (!dir.exists(saving.folder)) dir.create(saving.folder)
                        image_saver(result$img.report, "scratch_area", saving.folder, base.name = bn)
                        cat("Done!\n\n")
                }
                new.areas.line <- cbind(data.frame(name = bn), result$num.report)
                areas <- plyr::rbind.fill(areas, new.areas.line)
                write.table(x = areas, file = filenamer(output.folder, base.name = "report", ext = "csv", replace = TRUE),
                          sep = ";", dec = ",", row.names = FALSE)
                rm(result)
        }
        return(areas)
}

