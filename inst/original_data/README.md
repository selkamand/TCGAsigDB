## **TCGA Pan-Cancer Somatic Mutations Dataset**

Page Describing Download (pass-only mc3 MAF variants)

[https://xenabrowser.net/datapages/?cohort=TCGA%20Pan-Cancer%20(PANCAN)&removeHub=http%3A%2F%2F127.0.0.1%3A7222](https://xenabrowser.net/datapages/?cohort=TCGA%20Pan-Cancer%20(PANCAN)&removeHub=http%3A%2F%2F127.0.0.1%3A7222){.uri}

TCGA Unified Ensemble "MC3" mutation calls. This work is described in Scalable Open Science Approach for Mutation Calling of Tumor Exomes Using Multiple Genomic Pipelines. Cell Syst. 2018 Mar 28;6(3):271-281.e7. doi: 10.1016/j.cels.2018.03.002. PubMed PMID: 29596782.

**Downloaded**

```         
aria2c "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/mc3.v0.2.8.PUBLIC.xena.gz
```

**Processing**

Rename so its a

```{r}
data.table::fread("mc3.v0.2.8.PUBLIC.xena.gz") |>
  dplyr::rename(
    Tumor_Sample_Barcode = sample,
    Chromosome = chr,
    Start_Position = start,
    End_Position = end,
    Reference_Allele = reference,
    Tumor_Seq_Allele2 = alt,
    Hugo_Symbol = gene,
    Consequence = effect,
    HGVSp_Short = Amino_Acid_Change, 
    ) 
```
