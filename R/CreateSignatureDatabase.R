#' Create TCGA pan-cancer-atlas signature database
#'
#' Creates an sqlite database describing the mutational signatures present in every TCGA patient
#'
#' @param outdir output directory
#' @param ids a character vector of TCGA IDs to analyse. If NULL will run on all TCGA samples. See [tcga_fetch_metadata()] for a full list of possible IDs
#' @param strict throw an error when supplied IDs are not actually present in mutation dataset. If falws, will just warn and continue creating the signature database on those IDs which are present.
#' @return run for its side-effects. Invisibly returns NULL
#' @export
#'
#' @examples
#' \dontrun{
#' create_tcga_pancan_database(ids = c('TCGA-CA-6717-01', 'TCGA-A2-A0T5-01', 'TCGA-CF-A9FF-01'))
#' }
create_tcga_pancan_database <- function(outdir = getwd(), ids = NULL, strict = FALSE){

  if (!requireNamespace("R.utils", quietly = TRUE))
    cli::cli_abort("Package \"R.utils\" must be installed to use this function.")

  cli::cli_progress_step('Streaming In TCGA PanCancerAtlas data')
  df_pancan_meta <- tcga_fetch_metadata(ids = ids)
  df_pancan_meta <- df_pancan_meta |>
    dplyr::filter(sampleId %in% ids)

  maf <- tcga_fetch_pancan_maf(ids = ids)


  cli::cli_progress_step('Running Signature Analysis and Creating SQLITE database')
  create_database(maf = maf, outdir = outdir, prefix = "TCGA_mc3_pancanatlas", ref = "hg19", metadata = df_pancan_meta)

  return(invisible(NULL))
}

#' Fetch TCGA metadata
#'
#' @inheritParams create_tcga_pancan_database
#'
#' @return data.frame with TCGA metadata
#' @export
#'
tcga_fetch_metadata <- function(ids = NULL){
  df_pancan_meta <- data.table::fread("https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz")
  df_pancan_meta <- dplyr::select(df_pancan_meta, sampleId = sample, disease = `_primary_disease`, description = sample_type)

  # Filter for only IDs of interest
  if(!is.null(ids)){
    df_pancan_meta <- df_pancan_meta |>
      dplyr::filter(sampleId %in% ids)
  }


  return(df_pancan_meta)
}

#' Fetch MAF
#'
#' Fetch TCGA MC3 MAF
#'
#' @inheritParams create_tcga_pancan_database
#'
#' @return maf data.frame
#' @export
#'
tcga_fetch_pancan_maf <- function(ids = NULL, strict = FALSE){
  df_pancan <- data.table::fread("https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/mc3.v0.2.8.PUBLIC.xena.gz")

  # Filter for IDs we care about
  if (!is.null(ids)) {
    assertions::assert_character(ids)
    assertions::assert_no_duplicates(ids)
    df_pancan <- df_pancan |>
      dplyr::filter(sample %in% ids)


    # Ensure All the IDs we want are actually present in the data
    nsamples <- dplyr::n_distinct(df_pancan$sample)

    if (nsamples != length(ids) & strict)
      cli::cli_abort("Only {nsamples} / {length(ids)} of the TCGA IDs supplied were observed in the data. To continue anyway, rerun with {.arg strict = TRUE}")
    else
      cli::cli_alert_warning("All {nsamples} / {length(ids)} of the TCGA IDs supplied were observed in the data")
  }

  cli::cli_progress_step('Converting to MAF format')
  maf <- df_pancan |>
    dplyr::rename(
      Tumor_Sample_Barcode = sample,
      Chromosome = chr,
      Start_Position = start,
      End_Position = end,
      Reference_Allele = reference,
      Tumor_Seq_Allele2 = alt,
      Hugo_Symbol = gene,
      Consequence = effect,
      HGVSp_Short = Amino_Acid_Change
    )

  return(maf)
}




#' Create Sqlite database
#'
#' @param ref which human reference genome to use
#' @param maf a maf object, path to maf file, or data.frame in MAF format
#' @inheritParams sigminerUtils::sig_analyse_mutations
#' @param outdir outfile directory
#' @param prefix prefix of output filenames
#'
#' @return invisible(NULL)
#'
#' @examples
#' \dontrun{
#' path_maf <- system.file("original_data/mc3.v0.2.8.PUBLIC.xena.gz", package = "TCGAsigDB")
#'
#' # Fix Column Names
#' maf <- data.table::fread(path_maf) |>
#'   dplyr::rename(
#'     Tumor_Sample_Barcode = sample,
#'     Chromosome = chr,
#'     Start_Position = start,
#'     End_Position = end,
#'     Reference_Allele = reference,
#'     Tumor_Seq_Allele2 = alt,
#'     Hugo_Symbol = gene,
#'     Consequence = effect,
#'     HGVSp_Short = Amino_Acid_Change
#'   )
#'
#' create_database(maf, ref = 'hg19', prefix = "TCGA_mc3_pancanatlas")
#'
#' # Or we can pull the data from maftools
#' mafs <- lapply(maftools::tcgaAvailable()[['Study_Abbreviation']], maftools::tcgaLoad)
#' merged_mafs <- maftools::merge_mafs(mafs)
#'
#'
#' create_database(merged_mafs, ref = 'hg19', prefix = "TCGA_mc3_maftools")
#' }
create_database <- function(maf, ref = c("hg38", "hg19"), metadata = NULL, outdir = getwd(), prefix = "cosmic_signatures"){

  # Set up
  if(ref == "hg19"){
    ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
  }else if (ref == "hg38"){
    ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
  }else
    stop("Unknown reference genome", ref)

  if(!is.character(maf) & !is.data.frame(maf))
    assertions::assert_class(maf, class = "MAF")

  if(is.data.frame(maf))
    maf <- sigminer::read_maf_minimal(maf)

  if(is.character(maf))
    maf <- sigminer::read_maf_minimal(data.table::fread(maf))


  cli::cli_progress_step('analysing mutations')

  # Sigminer
  sigminerUtils::sig_analyse_mutations(
    maf = maf,
    ref = ref,
    output_dir = glue::glue("{outdir}/signatures"),
    exposure_type = "absolute",
    n_bootstraps = 100,
    temp_dir = tempdir()
  )


  #sigminerUtils::sig_create_database(sqlite_db = "cosmic_signatures.hg19.sqlite", overwrite = TRUE)
  cli::cli_progress_step('creating database')
  sigminerUtils::sig_create_database(sqlite_db = glue::glue("{outdir}/{prefix}.{ref}.sqlite"), overwrite = TRUE)

  cli::cli_progress_step('Adding to  database')
  sigminerUtils::sig_add_to_database(
    signature_directory = glue::glue("{outdir}/signatures"),
    sqlite_db = as.character(glue::glue("{outdir}/{prefix}.{ref}.sqlite")),
    ref = ref,
    metadata = metadata
  )
}


