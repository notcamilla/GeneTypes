## ----style, echo=FALSE, message=FALSE-----------------------------------------
library(BiocStyle)
library(knitr)
library(GeneTypes)

## ----loading packages, message=FALSE, warning=FALSE---------------------------
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)

## ----warning message----------------------------------------------------------
# duplicated_tx_names <- duplicated(names(cds_bytx))

## ----suppressWarnings usage---------------------------------------------------
# cds_bytx <- suppressWarnings(cdsBy(txdb, "tx", use.names=TRUE))
# protein_seqs <- sapply(exons_seqs_bytx, function(seq) suppressWarnings(translate(seq)))

## ----protein coding gene example----------------------------------------------
FAM3C <- ProteinCodingGene(id=10447, symbol= 'FAM3C', name='FAM3 metabolism regulating signaling molecule C',
                    description = 'encodes a secreted protein with a GG domain', 
                    structure = GRanges("chr7", IRanges(121348878, 121396364), 
                    strand = '-'), proteinID = 'Q92520')

## ----get structure example----------------------------------------------------
structure(FAM3C)

## ----length product example---------------------------------------------------
lengthProduct(FAM3C)

## ----trna gene example--------------------------------------------------------
TRAAGC11 <- tRNAGene(id=100189346, symbol= 'TRA-AGC1-1', name='tRNA-Ala (anticodon AGC) 1-1',
                    description = 'tRNA-Ala (anticodon AGC) 1-1', 
                    structure = GRanges("chr6", IRanges(28795964, 28796035), strand = "-"),
                    tRNAID='URS000063E4FD_9606', anticodon='AGC')

## ----get name example---------------------------------------------------------
name(TRAAGC11)

## ----length product example two-----------------------------------------------
lengthProduct(TRAAGC11)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

