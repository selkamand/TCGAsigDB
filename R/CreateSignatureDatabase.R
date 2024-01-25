


#' Create Sqlite database
#'
#' @param path_maf path to a maf file (e.g. pancancer MAF)  or a maf object
#' @param ref which human reference genome to use
#'
#' @return
#'
#' @examples
#' \dontrun{
#' path_maf <- system.file("original_data/mc3.v0.2.8.PUBLIC.xena.gz", package = "TCGAsigDB")
#'
#' # Fix Column Names
#' maf = data.table::fread(path_maf) |>
#  dplyr::rename(
#    Tumor_Sample_Barcode = sample,
#    Chromosome = chr,
#    Start_Position = start,
#    End_Position = end,
#    Reference_Allele = reference,
#    Tumor_Seq_Allele2 = alt,
#    Hugo_Symbol = gene,
#    Consequence = effect,
#    HGVSp_Short = Amino_Acid_Change
#  )

#' create_database(maf, ref = 'hg19', prefix = "TCGA_mc3_pancanatlas")
#' }
create_database <- function(maf, ref = c("hg38", "hg19"), outdir = getwd(), prefix = "cosmic_signatures"){

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
    sigminer::read_maf_minimal(maf = maf)

  if(is.character(maf))
    maf <- sigminer::read_maf_minimal(maf = data.table::fread(maf))


  sample_ids = unique(maf@data[['Tumor_Sample_Barcode']])

  # Sigminer
  sigminerUtils::sig_analyse_mutations(
    maf = maf,
    somatic_ids = sample_ids,
    ref = ref,
    output_dir = "{outdir}/signatures",
    exposure_type = "absolute",
    n_bootstraps = 100,
    temp_dir = tempdir()
  )


  #sigminerUtils::sig_create_database(sqlite_db = "cosmic_signatures.hg19.sqlite", overwrite = TRUE)
  sigminerUtils::sig_create_database(sqlite_db = glue::glue("{outdir}/{prefix}.{ref}.sqlite"), overwrite = TRUE)
  sigminerUtils::sig_add_to_database(signature_directory = "{outdir}/signatures", sqlite_db = as.character(glue::glue("{outdir}/{prefix}.{ref}.sqlite")),ref = ref)
}


