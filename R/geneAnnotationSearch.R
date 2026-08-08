#' Search for gene annotations on positively classified contigs
#'
#' Search contigs classified with TrIdent for gene-annotations that match a provided
#' key-word(s). Outputs read coverage plots for contigs with matching annotations.
#'
#' @param TrIdentResults The output from `TrIdent()`.
#' @param VLPpileup A .txt file containing mapped sequencing read coverages averaged over
#' 100 bp windows/bins.
#' @param inPatMat TRUE or FALSE. If TRUE, only search for gene-annotations in
#' the pattern-match region. Default is FALSE (i.e search the
#' entire contig for the gene annotation key-words)
#' @param bpRange If `inPatMat` = TRUE, the user may specify the region (in base pairs) that should
#' be searched to the left and right of the pattern-match region. Default is 0.
#' @param gff A .gff file imported using Bioc.gff::readGFF that is 
#' associated with the whole-community contigs
#' @param searchCol The exact column name (in quotes) in the gff file that you'd like to search for specific 
#' `keyWords`. Commonly annotation, product, gene, or acc columns, for example, however the specific
#' name may change based on annotation tool used. For example: "pfam_desc". 
#' @param keyWords The keyWord(s) to search for. Case independent. Searches will return the string
#' that contains the matching keyWord. KeyWord(s) must be in quotes, comma-separated, and surrounded by
#' c() i.e( c("antibiotic", "resistance", "drug") )
#' @param verbose TRUE or FALSE. Print progress messages to console. Default is TRUE.
#' @param saveFilesTo Optional, Provide a path to the directory you wish to save
#' output to. A folder will be made within the provided directory to store
#' results.
#' @returns list of ggplot objects
#' @importFrom dplyr %>%
#' @importFrom S4Vectors DataFrame
#' @export
#' @examples
#' data("VLPFractionSamplePileup")
#' data("TrIdentSampleOutput")
#' data("gffSample")
#'
#' patternMatches <- geneSearch(
#'   TrIdentResults = TrIdentSampleOutput,
#'   VLPpileup = VLPFractionSamplePileup,
#'   gff = gffSample,
#'   searchCol="pfam_desc",
#'   keyWords=c("toxin", "drug", "resistance", "phage")
#' )
#' 
geneSearch <- function(TrIdentResults, VLPpileup, gff,
                       searchCol, keyWords, inPatMat = FALSE,
                        bpRange = 0, saveFilesTo, verbose = TRUE) {
    start <- NULL
    if(bpRange != 0 & inPatMat == FALSE) {stop("Cannot set bpRange if inPatMat = FALSE")}
    if(searchCol %in% colnames(gff) == FALSE) {stop("The `searchCol` does not exist in the provided gff file")}
    summaryTable <- TrIdentResults[[1]]
    patternMatches <- TrIdentResults[[3]]
    windowSize <- TrIdentResults[[5]]
    if(verbose){message("Cleaning pileup file...")}
    VLPpileup <- pileupFormatter(VLPpileup)
    if(verbose){message("Searching for matching annotations...")}
    gff <- as.data.frame(gff)
    plots <- lapply(seq_along(patternMatches), function(i){
        refName <- patternMatches[[i]][[8]]
        classification <- patternMatches[[i]][[7]]
        if (inPatMat & classification == "NoPattern") {return(NULL)}
        VLPpileupSubset <- VLPpileup[which(VLPpileup[,1] == refName),]
        VLPpileupRegion <- seq(VLPpileupSubset[1,3], 
                               VLPpileupSubset[nrow(VLPpileupSubset),3], 1)
        gff <- gff[which(gff$seqid == refName), ]
        geneAnnotSubset <- gff[which(gff$end %in% VLPpileupRegion), ]
        colIdx <- which(colnames(geneAnnotSubset) == searchCol)
        if(typeof(geneAnnotSubset[,colIdx])=="list"){
        geneAnnotSubset <- unnest(geneAnnotSubset, searchCol)
        geneAnnotSubset <- as.data.frame(distinct(geneAnnotSubset, start, .keep_all=TRUE))
        }
        if (TRUE %in% (str_detect(geneAnnotSubset[,colIdx], 
                                  regex(paste(keyWords, collapse="|"), 
                                        ignore_case = TRUE)))) {
            startbpRange <- endbpRange <- NULL
            if (inPatMat){
                startPos <- patternMatches[[i]][[5]] * windowSize
                endPos <- patternMatches[[i]][[6]] * windowSize
                bpRange <- ifelse (classification == "HighCovNoPattern", 0, bpRange)
                endbpRange <- ifelse ((endPos + bpRange > (VLPpileupSubset[nrow(VLPpileupSubset), 3])),
                                      (VLPpileupSubset[nrow(VLPpileupSubset), 3]), (endPos + bpRange))
                startbpRange <- ifelse((startPos - bpRange < VLPpileupSubset[1, 3] ), 
                                       VLPpileupSubset[1, 3], (startPos - bpRange))
                matchRegion <- seq(startbpRange, endbpRange, 1)
                geneAnnotSubset <- as.data.frame(gff[which(gff$end %in% matchRegion), ])
                if (!(TRUE %in% (str_detect(geneAnnotSubset[,colIdx],
                                            regex(paste(keyWords, collapse="|"),
                                                  ignore_case = TRUE))))) {
                    return(NULL)
                }
            }
            plot <- geneAnnotationPlot(geneAnnotSubset, keyWords,
                                       VLPpileupSubset, colIdx, startbpRange, endbpRange,
                                       patternMatches[[i]], windowSize)
        } else {return(NULL)}
        plot
    })
    plots <- (plots[!vapply(plots, is.null, logical(1))])
    refNames <- vapply(seq_along(plots), function(i){
        plot <- plots[[i]]
        plotdata <- ggplot_build(plot)
        contig <- plotdata$plot$labels$title %>% str_extract(".+?(?=[:space:])")
        contig
    }, character(1))
    names(plots) <- refNames
    if(verbose){message(length(plots), 
                        " contigs have gene annotations that match one or more of the provided keyWords")}
    if (missing(saveFilesTo) == FALSE) {
        ifelse(!dir.exists(paths = paste0(saveFilesTo, "\\TrIdentGeneAnnotMatchPlots")),
               dir.create(paste0(saveFilesTo, "\\TrIdentGeneAnnotMatchPlots")),
               stop("'TrIdentGeneAnnotMatchPlots' already exists in the provided directory")
        )
        lapply(
            names(plots),
            function(X) {
                ggsave(
                    filename = paste0(
                        saveFilesTo,
                        "\\TrIdentGeneAnnotMatchPlots\\", X, ".png"
                    ),
                    plot = plots[[X]],
                    width = 8,
                    height = 4
                )
            }
        )
        return(plots)
    } else {
        return(plots)
    }
}






#' Gene annotation plot
#'
#' Plot read coverage and location of gene annotations that match the keywords and
#' search criteria for contig currently being assessed
#'
#' @param geneAnnotSubset Subset of gene annotations to be plotted
#' @param keywords The key-word(s) used for the search.
#' @param VLPpileupSubset A subset of the pileup associated with the contig/chunk being assessed
#' @param colIdx The gff column being searched
#' @param startbpRange The basepair at which the search is started if a 'specific' search is used
#' @param endbpRange The basepair at which the search is ended if a 'specific' search is used
#' @param pattern The pattern-match information associated with the contig/chunk being assessed
#' @param windowSize The number of basepairs to average read coverage values over.
#' @keywords internal
#' @importFrom stringr str_which
geneAnnotationPlot <- function(geneAnnotSubset, keywords, VLPpileupSubset,
                               colIdx, startbpRange, endbpRange,
                               pattern, windowSize){
    position <- coverage <- start <- NULL
    classification <- pattern[[7]]
    refName <- VLPpileupSubset[1, 1]
    startPos <- pattern[[5]] * windowSize
    endPos <- pattern[[6]] * windowSize
    matchIdxs <- str_which(geneAnnotSubset[,colIdx], regex(paste(keywords, collapse = "|"), ignore_case = TRUE))
    geneAnnotMatches <- geneAnnotSubset[matchIdxs,]
    geneStartPos <- geneAnnotMatches$start
    geneAnnotLabels <- paste0("#", c(1:nrow(geneAnnotMatches)), ": ", geneAnnotMatches[,colIdx],
                              sep = " ", collapse = " \n ")
    VLPpileupSubset <- changeWindowSize(VLPpileupSubset, windowSize)
    plot <- ggplot(data = VLPpileupSubset, aes(x = position, y = coverage)) +
        geom_area(fill = "#009E73") +
        geom_vline(xintercept = geneStartPos, linewidth = 1, alpha=0.5) +
        geom_vline(xintercept = c(startPos, endPos), color = "#D55E00", linewidth = 1, alpha=0.5) +
        geom_vline(xintercept = c(startbpRange, endbpRange), color = "#D55E00", linewidth = 1, linetype = "dotted", alpha=0.5) +
        geom_label(data = geneAnnotMatches, aes(x = start, y = (max(VLPpileupSubset$coverage) / 2),
                                                label = paste0("#", c(1 : nrow(geneAnnotMatches)))),
                   size = 2.75) +
        scale_x_continuous(expand = c(0, 0)) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              text = element_text(size = 15),
              plot.margin = margin(
                  t = 0,
                  r = 10,
                  b = 0,
                  l = 2
              )) +
        labs(title = paste(refName, classification),
             x = "Basepair position",
             caption = geneAnnotLabels,
             y = "Read coverage")
    return(plot)
}
