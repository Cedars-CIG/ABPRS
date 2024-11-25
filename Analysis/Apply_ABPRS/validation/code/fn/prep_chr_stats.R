require(data.table)
require(dplyr)
require(cli)

prep_chr_stats <- function(chromosome, snplist = "data/SNP_list.lifted", geno_dir = "geno/pmbb/", prefix = "", suffix = "_small") {
    cli::cli_alert_info("chromosome {chromosome}")
    # read and process snplist
    snplist <- data.table::fread(snplist) |>
        dplyr::filter(V1 == paste0("chr", chromosome)) |>
        dplyr::mutate(snp = paste0(V1, ":", V2))

    # read and process the hardy file
    hardy <- data.table::fread(paste0(geno_dir, prefix, "chr", chromosome, suffix, ".hardy")) |>
        dplyr::mutate(
            snp     = gsub("(:[^:]+){2}$", "", ID),
            hardy_p = as.numeric(P)
        ) |> dplyr::select(snp, hardy_p)
    multiallelic_snps <- hardy[, .N, snp][N > 1, snp]
    hardy <- hardy[!snp %in% multiallelic_snps,]

    # read and process .vmiss file
    vmiss <- data.table::fread(paste0(geno_dir, prefix, "chr", chromosome, suffix, ".vmiss")) |>
        dplyr::mutate(
            snp = gsub("(:[^:]+){2}$", "", ID)
        ) |>
        dplyr::select(snp, vmiss = F_MISS)
    multiallelic_snps <- vmiss[, .N, snp][N > 1, snp]
    vmiss <- vmiss[!snp %in% multiallelic_snps, ]


    # read and process .raw file
    d   <- data.table::fread(paste0(geno_dir, prefix, "chr", chromosome, suffix, ".raw")) |>
        select(-c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")) 
    d2 <- d |>
        dplyr::mutate_if(is.numeric, \(x) round(x, 0))
    diff <- data.table::as.data.table(abs(as.matrix(d) - as.matrix(d2)))
    diff_sum <- lapply(names(diff), \(x) {
        data.table::data.table(
            id         = x,
            mean_dist   = mean(diff[[x]], na.rm = TRUE),
            pct_75_dist = quantile(diff[[x]], 0.75),
            pct_90_dist = quantile(diff[[x]], 0.9),
            max_dist    = max(diff[[x]], na.rm = TRUE))
    }) |> data.table::rbindlist() |>
        dplyr::mutate(
            snp = gsub("(:[^:]+){2}$", "", gsub("_.*", "", id))
        )
    multiallelic_snps <- diff_sum[, .N, snp][N > 1, snp]
    diff_sum <- diff_sum[!(snp %in% multiallelic_snps), ]

    # read and process .afreq file
    maf <- data.table::fread(paste0(geno_dir, prefix, "chr", chromosome, suffix, ".afreq")) |>
        dplyr::mutate(
            var_id = paste0(ID, "_", REF)
        ) |>
        dplyr::mutate(
            ref_len = nchar(REF),
            alt_len = nchar(ALT),
            len = as.numeric(ref_len > 1 | alt_len > 1),
            snp = gsub("(:[^:]+){2}$", "", ID),
            aaf = ALT_FREQS,
            maf = ifelse(aaf <= 0.5, aaf, 1 - aaf),
            maf_cat = fcase(
                maf > 0.05, "Common",
                maf > 0.01, "Low-frequency",
                default = "Rare"
            )
        ) |>
        dplyr::filter(len == 0) 
    multiallelic_snps <- maf[, .N, snp][N > 1, snp]
    maf <- maf[!snp %in% multiallelic_snps, ]

    tmp <- data.table::merge.data.table(
        diff_sum,
        maf,
        by = "snp"
    )
    tmp <- data.table::merge.data.table(
        tmp,
        hardy,
        by = "snp"
    )
    tmp <- data.table::merge.data.table(
        tmp,
        vmiss,
        by = "snp"
    )
    tmp[snp %in% snplist[, snp], ] |>
        dplyr::select(snp, id, mean_dist, pct_75_dist, pct_90_dist, max_dist,
               REF, ALT, aaf, maf, maf_cat, info = MINIMAC3_R2,
               hardy_p, vmiss)
}
