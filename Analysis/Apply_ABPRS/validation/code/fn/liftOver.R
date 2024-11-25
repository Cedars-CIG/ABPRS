library(data.table)
library(tidyverse)

liftOver <- function(
    input_bed,
    chain = NULL,
    tmpdir = NULL
) {
    if (is.null(chain)) {
        message("`chain`: defaulting to hg19ToHg38.over.chain")
        chain <- "hg19ToHg38.over.chain"
    }
    if (is.null(tmpdir)) {
        message("`tmpdir`: defaulting `getwd()`")
        tmpdir <- getwd()
    }
    input_bed <- unique(input_bed)
    tmpdir <- paste0(tmpdir, "tmp/")
    dir.create(tmpdir)
    write.table(
        input_bed, paste0(tmpdir, "tmp.bed"),
        col.names = FALSE, row.names = FALSE, quote = FALSE
    )
    system(paste0(
        "/home/mmsal/group/datasets/ucsc-chainfiles/liftOver ",
            tmpdir, "tmp.bed ",
            "/home/mmsal/group/datasets/ucsc-chainfiles/", chain, " ",
            tmpdir, "lift ",
            tmpdir, "exclude "
    ))

    lift    <- suppressWarnings(data.table::fread(paste0(tmpdir, "lift")))
    exclude <- suppressWarnings(data.table::fread(paste0(tmpdir, "exclude")))

    bed <- data.table::as.data.table(input_bed)
    names(bed)  <- c("CHR_old", "POS_old", "POS1_old")
    names(lift) <- c("CHR", "POS", "POS1")
    lift[, POS := as.character(POS)]

    on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
    list(
        lift    = lift,
        exclude = exclude,
        map     = cbind(
            lift[, 1:2],
            bed[, 1:2]
        )
    )
}
