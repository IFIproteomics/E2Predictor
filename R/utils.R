

#' countCharOccurrences counts the presence of a given character into a string
#'
#' @param s string where chars should be counted
#' @param char character you want to count
#'
#' @export
countCharOccurrences <- function(s, char) {
    s2 <- gsub(char,"",s)
    return (nchar(s) - nchar(s2))
}

#' read_the_file is a wrapper of read_excel and read_delim (from readr and readxl packages) that
#' checks which is the proper delimiter for the file before reading it.
#'
#' @param thefile file to be read
#'
#' @export
#'
read_the_file <- function(thefile, ...){
    seps <- c(",", ";", "\t", ":")

    if(file_ext(thefile) == "xls" | file_ext(thefile) == "xlsx" )
    {
        df <- read_excel(thefile, ...)
        return(df)
    }

    # try different delims
    trydf <- suppressWarnings(read_delim(thefile, delim = "", n_max = 5))

    tryDelim <- function(delim){
        numOccurs <- sapply(trydf[, 1], countCharOccurrences, delim)
        if(! length(unique(numOccurs)) == 1 ) return(FALSE)
        if(max(numOccurs) == 0) return(FALSE)

        return(TRUE)
    }

    ll <- sapply(seps, tryDelim)

    if(table(ll)["TRUE"] > 1)
        stop(paste("Ambiguous file (", thefile ,"). Can't find a right delimiter to parse it! (several delimiters work)"))

    if(table(ll)["TRUE"] == 0)
        stop(paste("File", thefile, "is not delimited by any of the common delimiters!"))

    df = read_delim(thefile, delim= names(which(ll)), ... )

    return(df)
}


#' read_file_add_filename reads a file and attaches a column with the file name or the folder name
#'
#' @param thefile the file to be read
#' @param keepFolderName uses the folder name instead of the file name
#' @param UseGenericHeaders substitute headers by generic headers (useful when concatenating several files and headers are different)
#'
#' @export
#'
read_file_add_filename <- function(thefile, keepFolderName = F, UseGenericHeaders=F){

    print(paste("reading file:", thefile))

    if(length(thefile)==0){
        return(NULL)
    }

    csvdata <- read_file(thefile)

    foldername <- strsplit(thefile, "/")[[1]]
    foldername <- foldername[length(foldername) - 1]

    filename <- gsub(".*/", "", thefile)

    if(UseGenericHeaders){
        names(csvdata) <- paste0("X", 1:length(names(csvdata)) )
    }

    csvdata$FileName <-  filename

    if(keepFolderName){
        csvdata$FileName <- foldername
    }

    return(csvdata)
}

#' read_batch_in_folder reads several files from a folder, and concatenates them
#'
#' @param folder folder where files will be read
#' @param filepattern regular expression to filter files in the folder by file name
#' @param keepFolderName keeps the folder name in a column instead of the file name
#' @param UseGenericHeaders substitute headers by generic headers (useful when concatenating several files and headers are different)
#'
#' @export
#'
read_batch_in_folder <- function(folder,filepattern=".", keepFolderName = F, UseGenericHeaders=T){

    thefiles <- list.files(folder, filepattern, full.names = T)

    if(length(thefiles)==0){
        return(NULL)
    }

    csvdata <- do.call("rbind", lapply(thefiles, read_file_add_filename, keepFolderName, UseGenericHeaders))

    return(csvdata)
}


cal.AverageSample <- function(sample, df){

    df.sample <- df[, grep(sample, names(df))]
    df.avg <- data.frame( avg = rowMeans(df.sample, na.rm = T) )

    return(df.avg)
}

cv <- function(x) ( sd(x, na.rm = T) / mean(x, na.rm = T))

cal.CVs <- function(sample, df){

    df.values <- df[, grep(sample, names(df))]
    cvs <- data.frame( CV = apply(df.values, 1, cv) )

    return(cvs)
}

cal.FDR <- function(values, w, slicing = NULL){
    #values <- log2AB

    log2.mean.weighted = sum(values*w, na.rm = T) / sum(w, na.rm = T)
    #slicing = 500
    phi = 0.6745

    robust_var <- function(x){

        rv <- median((x - mean(x, na.rm = T))^2, na.rm = T) / phi^2
    }

    robust_sd <- function(x){ sqrt(robust_var(x)) }

    df <- data.frame( index = 1:length(values), value = values, weight = w ) %>%
        arrange(desc(weight))

    df$mean <- log2.mean.weighted
    df$stddev <-  sd(df$value, na.rm = T)

    if(!is.null(slicing)){

        stdev2 <-  rollapply(df$value, slicing, robust_sd )
        stdev2 <- c( rep(stdev2[1], round(slicing/2, 0)), stdev2)
        stdev2 <- c( stdev2, rep(stdev2[length(stdev2)], round(slicing/2,0) ) )
        stdev2 <- stdev2[1:nrow(df)]

        mean_roll <- rollapply(df$value, slicing, mean, na.rm=T)
        mean_roll <- c( rep(mean_roll[1], round(slicing/2, 0)), mean_roll)
        mean_roll <- c( mean_roll, rep(mean_roll[length(mean_roll)], round(slicing/2, 0)))
        mean_roll <- mean_roll[1:nrow(df)]

        df$stddev <- stdev2
        df$mean <- mean_roll
    }

    #plot(log2(df$w), df$stddev2, type="l")

    df$pvalue <- pnorm( abs(df$value), mean = df$mean, sd = df$stddev, lower.tail = F)

    df <- df %>% arrange(pvalue)

    total_elements = nrow(df[!is.na(df$pvalue),])  # do not count NA values

    df <- df %>% mutate(rank = min_rank(pvalue)) %>%
        mutate(FDR = pvalue * total_elements / rank  ) %>%
        arrange(index)

    return(df$FDR)
}


panel.logAB <- function (x, y, graphtitle ="", reportOutliers=T, saveIndividualPlots=F, w_threshold= 0.0, fdr_threshold=0.01, ...) {
    # x = dt.Conditions$Bortezomib.avg
    # y = dt.Conditions$Control.avg

    w = x + y
    counter <<- counter + 1
    #print(counter)
    log2AB = log(x/y, base = 2)
    log2AB[!is.finite(log2AB)]<- NA

    fdr = cal.FDR(log2AB, w, slicing = 250)

    colors.log = ifelse(fdr <= fdr_threshold & w > w_threshold, "red", "grey")
    colors.log[log2AB >= 0 & fdr <= fdr_threshold & w > w_threshold] = "green"
    pch.log = rep(20, length(colors.log))
    pch.log[colors.log == "red" | colors.log == "green"] = 19

    outliersID <- as.data.frame(cbind(proteinIDs, colors.log, log2AB))

    names(outliersID) <- c("proteinID", "regulation", "log2.ratio")
    outliersID$regulation <- gsub("red", "down-regulated", outliersID$regulation)
    outliersID$regulation <- gsub("green", "up-regulated", outliersID$regulation)

    outliersID <- outliersID %>%
        filter(regulation == "down-regulated"| regulation == "up-regulated")

    if(reportOutliers){
        write_tsv(outliersID,
                  path = file.path(proteins.outputBaseDir, paste(Cell.Line, "outliers", graphtitle, counter , "tsv", sep = ".")))
    }

    if(saveIndividualPlots){
        pdf(file.path(proteins.figuresBaseDir, paste(Cell.Line, "plot", graphtitle, counter, "pdf", sep = ".")))
        currdev = dev.cur()
        par.old = par()
        par(mar=c(5,5,3,1))
        plot(log2AB, log(w, base = 2),
             col = colors.log,
             xlab=expression('log'[2]*'(A:B)'),
             ylab=expression('log'[2]*'(A+B)'),
             pch= pch.log, cex = 1.5, cex.lab = 1.5)
        par = par.old
        dev.off(which = currdev)
    }


    points(log2AB, log(w, base = 2), col = colors.log, ...)
    #abline(lm(y~x), col="red")
    #lines(stats::lowess(y~x), col="blue")
}


modsInSequence <- function(sequence, mods){

    seq_chars <- str_match_all(sequence, ".{1}")[[1]]

    mods_positions <- str_extract_all(mods, "\\([^()]+\\)")[[1]]
    mods_positions <- as.numeric(substring(mods_positions, 2, nchar(mods_positions)-1))
    mods <- gsub("\\(.*?\\)", "", mods)
    mods <- strsplit(mods, ",")[[1]]

    mods <- gsub("Carbamidomethyl C", "(+57.02)", mods)
    mods <- gsub("Oxidation M", "(+15.99)", mods)
    mods <- gsub("Phosphoryl STY", "(+79.97)", mods)
    mods <- gsub("Phosphorylation S", "(+79.97)", mods)
    mods <- gsub("Phosphorylation T", "(+79.97)", mods)
    mods <- gsub("Phosphorylation Y", "(+79.97)", mods)
    mods <- gsub(" ", "", mods)

    modifications <- data.frame( mod_pos = as.numeric(mods_positions), mod = mods)
    modlist <- 1:nrow(modifications)

    sapply(modlist, function(x){ seq_chars[modifications[x,1]] <<- paste0(seq_chars[modifications[x,1]], modifications[x,2]) } )

    seq_with_mod <- ""
    sapply(seq_chars, function(x){ seq_with_mod <<- paste0(seq_with_mod, x) } )

    return(seq_with_mod)
}
