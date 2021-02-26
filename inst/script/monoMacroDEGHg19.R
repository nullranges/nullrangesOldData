## RNA-seq Differential Gene Expression Analysis for hg19

## Load required libraries
library(data.table)
library(plyranges)
library(tximeta)
library(readr)
library(DESeq2)

## Read in sample sheet as coldata
coldata <-
  fread("inst/extdata/rna/samplesheet.tsv") %>%
  as.data.frame()

## Edit quant paths for hg19
coldata$files <- gsub("rna/quants", "rna/hg19/quants", coldata$files)

## Check that quant paths are correct
file.exists(coldata$files)

## Load Linked Txome
loadLinkedTxome("inst/extdata/rna/hg19/GENCODE.v19.salmon_1.4.0.LinkedTxome.json")

## Import data with tximeta & summarize to gene
se <- tximeta(coldata)
gse <- summarizeToGene(se)

## Add gene symbols
library(org.Hs.eg.db)
gse <- addIds(gse, column = "SYMBOL", gene = T)

## Convert to factors
colData(gse)$Bio_Rep <- as.factor(colData(gse)$Bio_Rep)
colData(gse)$Condition <- as.factor(colData(gse)$Condition)

## Build DESeq object
dds <- DESeqDataSet(gse, design = ~Bio_Rep + Condition)

## Filter out lowly expressed genes (at least 10 counts in at least 4 samples)
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep,]

## Fit model
dds <- DESeq(dds)

## Get DE gene results as GRanges
de_genes <- results(dds,
                    contrast = c("Condition", "PMA", "NON"), 
                    format = "GRanges") %>%
  names_to_column("gene_id")

## Add gene symbols to de_genes
de_genes$gene_symbol <-
  rowData(gse)$SYMBOL[match(de_genes$gene_id, rowData(gse)$gene_id)]

## Change data name
monoMacroDEGHg19 <- de_genes

## Save results to a file
save(monoMacroDEGHg19, file = "data/monoMacroDEGHg19.rda")
