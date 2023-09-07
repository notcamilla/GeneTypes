#' Get seed region for the MiRNA Gene.
#'
#' This function retrieves the seed region for a MiRNA Gene.
#'
#' @param x The MiRNAGene object.
#' @return A seed region object.
#' @export
setGeneric("getseedseq", function(x) standardGeneric("getseedseq"))

#' @rdname getseedseq
#' @param x The MiRNAGene object.
#' @method getseedseq MiRNAGene
#' @import microRNA
#' @import utils
#'
#' @examples
#' \dontrun{
#' gene3 <- MiRNAGene(id=406988, symbol= 'MIR205', name='MICRORNA 205',
#' description = 'non-coding RNAs that are involved in post-transcriptional
#' regulation of gene expression by affecting both the stability and translation
#' of mRNAs.', structure = GRanges("chr1", IRanges (209432133, 209432242),
#' strand = '+'), final_mature_miRNAID= 'hsa-miR-205')
#' data(hsSeqs)
#' getseedseq(gene3)
#' }
#' @export
setMethod("getseedseq", signature = "MiRNAGene", function(x) {
  data(hsSeqs)
  return(seedRegions(hsSeqs[x@final_mature_miRNAID]))
})


#' Get the mature RNA sequence for the MiRNA Gene class.
#'
#' This function retrieves the mature RNA sequence for a MiRNA Gene.
#'
#' @param x The MiRNAGene object.
#' @return An RNAString object representing the RNA sequence.
#' @export
setGeneric("matureMiRNAseq", function(x) standardGeneric("matureMiRNAseq"))

#' @rdname matureMiRNAseq
#' @param x The MiRNAGene object.
#' @method matureMiRNAseq MiRNAGene
#' @import Biostrings
#'
#' @examples
#' \dontrun{
#' gene3 <- MiRNAGene(id=406988, symbol= 'MIR205', name='MICRORNA 205',
#' description = 'non-coding RNAs that are involved in post-transcriptional
#' regulation of gene expression by affecting both the stability and translation
#' of mRNAs.', structure = GRanges("chr1", IRanges (209432133, 209432242),
#' strand = '+'), final_mature_miRNAID= 'hsa-miR-205')
#' data(hsSeqs)
#' matureMiRNAseq(gene3)
#' }
#' @export
setMethod("matureMiRNAseq", signature = "MiRNAGene", function(x) {
  data(hsSeqs)
  return(RNAString(hsSeqs[x@final_mature_miRNAID]))
})


#' Get the protein sequence for the Protein Coding Gene class.
#'
#' This function retrieves the protein sequence for Protein Coding Gene object.
#'
#' @param x The Protein Coding Gene object.
#' @return A list of protein sequences corresponding to the exonic
#' regions of the gene.
#' @export
setGeneric("proteinseq", function(x) standardGeneric("proteinseq"))

#' @rdname proteinseq
#' @param x The Protein Coding Gene object.
#' @method proteinseq ProteinCodingGene
#' @import GenomicFeatures
#' @import GenomicRanges
#' @import Biostrings
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @export
#' @examples gene1 <- ProteinCodingGene(id=10447, symbol= 'FAM3C', name='FAM3 metabolism
#' regulating signaling molecule C',description = 'encodes a secreted protein
#' with a GG domain', structure = GRanges("chr7",IRanges(121348878, 121396364),
#' strand = '-'), proteinID = 'Q92520')
#' proteinseq(gene1)
setMethod("proteinseq", signature = "ProteinCodingGene", function(x) {
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  genome <- BSgenome.Hsapiens.UCSC.hg38
  gene_txs <- transcriptsBy(txdb, 'gene')[as.character(x@id)]
  tx_names <- mcols(gene_txs[[1]])$tx_name
  cds_bytx <- suppressWarnings(cdsBy(txdb, "tx", use.names=TRUE))
  duplicated_tx_names <- duplicated(names(cds_bytx))
  cds_bytx_unique <- cds_bytx[!duplicated_tx_names]
  gene_cds_bytx <- cds_bytx[intersect(tx_names, names(cds_bytx_unique))]
  exons_seqs_bytx <- getSeq(genome, gene_cds_bytx)
  protein_seqs <- sapply(exons_seqs_bytx, function(seq) suppressWarnings(translate(seq)))
  return(protein_seqs)
})


#' Get the RNA sequence for the LncRNA Gene.
#'
#' This function retrieves the RNA sequence for a LncRNAGene object.
#'
#' @param x The LncRNAGene object.
#' @return An RNAString object representing the RNA sequence of the gene.
#' @export
setGeneric("LncRNAseq", function(x) standardGeneric("LncRNAseq"))

#' @rdname LncRNAseq
#' @param x The LncRNAGene object.
#' @method LncRNAseq LncRNAGene
#' @import Biostrings
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @export
#' @examples gene2 <- LncRNAGene(id=10984, symbol= 'KCNQ1OT1', name='antisense
#' transcript 1', description = 'interacts with chromatin and regulates
#' transcription of multiple target genes through epigenetic modifications',
#' structure = GRanges("chr11", IRanges(2608328, 2699994), strand = '-'),
#' lncRNAID='URS0000759CF4')
#' LncRNAseq(gene2)
setMethod("LncRNAseq", signature = "LncRNAGene", function(x) {
  genome <- BSgenome.Hsapiens.UCSC.hg38
  DNA_seq <- getSeq(genome, x@structure)
  RNA_seq <- RNAString(DNA_seq[[1]])
  return(RNA_seq)
})

#' Get the RNA sequence for the tRNA Gene.
#'
#' This function retrieves the RNA sequence for a tRNAGene object.
#'
#' @param x The tRNAGene object.
#' @return An RNAString object representing the RNA sequence of the gene.
#' @export
setGeneric("tRNAseq", function(x) standardGeneric("tRNAseq"))

#' @rdname tRNAseq
#' @param x The tRNAGene object.
#' @method tRNAseq tRNAGene
#' @import Biostrings
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @export
#' @examples gene5 <- tRNAGene(id=100189346, symbol= 'TRA-AGC1-1',name='tRNA-Ala
#' (anticodon AGC) 1-1', description = 'tRNA-Ala (anticodon AGC) 1-1',
#' structure = GRanges("chr6", IRanges(28795964, 28796035), strand = "-"),
#' tRNAID='URS000063E4FD_9606', anticodon='AGC')
#' tRNAseq(gene5)
setMethod("tRNAseq", signature = "tRNAGene", function(x) {
  genome <- BSgenome.Hsapiens.UCSC.hg38
  DNA_seq <- getSeq(genome, x@structure)
  RNA_seq <- RNAString(DNA_seq[[1]])
  return(RNA_seq)
})

#' Get the RNA sequence for the rRNA Gene
#'
#' This function retrieves the RNA sequence for an rRNAGene object.
#'
#' @param x The rRNAGene object.
#' @return An RNAString object representing the RNA sequence of the gene.
#' @export
setGeneric("rRNAseq", function(x) standardGeneric("rRNAseq"))

#' @rdname rRNAseq
#' @param x The rRNAGene object.
#' @method rRNAseq rRNAGene
#' @import Biostrings
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @export
#' @examples gene4 <- rRNAGene(id=100008588, symbol= 'RNA18SN5', name='RNA,
#' 18S ribosomal N5',description = 'the 45S rDNA repeat unit encodes a 45S rRNA
#' precursor, transcribed by RNA polymerase I', structure = GRanges("chr21",
#' IRanges(8436876, 8438744), strand = "+"), rRNAID='URS0000726FAB')
#' rRNAseq(gene4)
setMethod("rRNAseq", signature = "rRNAGene", function(x) {
  genome <- BSgenome.Hsapiens.UCSC.hg38
  DNA_seq <- getSeq(genome, x@structure)
  RNA_seq <- RNAString(DNA_seq[[1]])
  return(RNA_seq)
})

#' Get the RNA sequence for the SNRNA Gene.
#'
#' This function retrieves the RNA sequence for an SNRNAGene object.
#'
#' @param x The SNRNAGene object.
#' @return An RNAString object representing the RNA sequence of the gene.
#' @export
setGeneric("SNRNAseq", function(x) standardGeneric("SNRNAseq"))

#' @rdname SNRNAseq
#' @param x The SNRNAGene object.
#' @method SNRNAseq SNRNAGene
#' @import Biostrings
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @export
#' @examples gene6 <- SNRNAGene(id=126871, symbol= 'RNU1-1', name='RNU1-1',
#' description = 'RNA, U1 small nuclear 1', structure = GRanges("chr1",
#' IRanges(16514122, 16514285), strand = "-"), SNRNAID='URS000032B6B6')
#' SNRNAseq(gene6)
setMethod("SNRNAseq", signature = "SNRNAGene", function(x) {
  genome <- BSgenome.Hsapiens.UCSC.hg38
  DNA_seq <- getSeq(genome, x@structure)
  RNA_seq <- RNAString(DNA_seq[[1]])
  return(RNA_seq)
})


#' Get the mature mRNA sequence for the Coding Gene
#'
#' This function retrieves the mature mRNA sequence for
#' a Protein Coding Gene object.
#'
#' @param x The Protein Coding Gene object.
#' @return An RNAStringSetList object representing the mature
#' mRNA sequences of the gene.
#' @export
setGeneric("mRNAseq", function(x) standardGeneric("mRNAseq"))

#' @rdname mRNAseq
#' @param x The ProteinCodingGene object.
#' @method mRNAseq ProteinCodingGene
#' @import Biostrings
#' @import GenomicFeatures
#' @import GenomicRanges
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @export
#' @examples gene1 <- ProteinCodingGene(id=10447, symbol= 'FAM3C', name='FAM3
#' metabolism regulating signaling molecule C',description = 'encodes a secreted
#' protein with a GG domain', structure = GRanges("chr7",IRanges(121348878,
#' 121396364), strand = '-'), proteinID = 'Q92520')
#' mRNAseq(gene1)
setMethod("mRNAseq", signature = "ProteinCodingGene", function(x) {
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  genome <- BSgenome.Hsapiens.UCSC.hg38
  gene_txs <- transcriptsBy(txdb, 'gene')[as.character(x@id)]
  tx_names <- mcols(gene_txs[[1]])$tx_name
  cds_bytx <- suppressWarnings(cdsBy(txdb, "tx", use.names=TRUE))
  duplicated_tx_names <- duplicated(names(cds_bytx))
  cds_bytx_unique <- cds_bytx[!duplicated_tx_names]
  gene_cds_bytx <- cds_bytx[intersect(tx_names, names(cds_bytx_unique))]
  exons_seqs_bytx <- getSeq(genome, gene_cds_bytx)
  return(RNAStringSetList(exons_seqs_bytx))
})


#' Calculate the length product for a protein-coding gene.
#'
#' This function calculates the length product for a Coding Gene object,
#' which is the product of the lengths of all its protein sequences.
#'
#' @param x The ProteinCodingGene object.
#' @return The length product of the protein sequences.
#' @export
setGeneric("lengthProduct", function(x) standardGeneric("lengthProduct"))

#' @rdname lengthProduct
#' @param x The ProteinCodingGene object.
#' @method lengthProduct ProteinCodingGene
#' @export
#' @examples gene1 <- ProteinCodingGene(id=10447, symbol= 'FAM3C', name='FAM3
#' metabolism regulating signaling molecule C',description = 'encodes a secreted
#' protein with a GG domain', structure = GRanges("chr7",IRanges(121348878,
#' 121396364), strand = '-'), proteinID = 'Q92520')
#' lengthProduct(gene1)
setMethod("lengthProduct", signature = "ProteinCodingGene", function(x) {
  protein_seqs <- proteinseq(x)
  aa_length <- vapply(protein_seqs, function(seq) sum(nchar(seq)), integer(1))
  return(aa_length)
})


#' Calculate the length product for a long non-coding RNA gene.
#'
#' This function calculates the length product for a LncRNAGene object,
#' which is the sum of the widths of its structure (GRanges) intervals.
#'
#' @param x The LncRNAGene object.
#' @return The length product of the gene structure.
#' @export
#' @examples gene2 <- LncRNAGene(id=10984, symbol= 'KCNQ1OT1', name='antisense
#' transcript 1', description = 'interacts with chromatin and regulates
#' transcription of multiple target genes through epigenetic modifications',
#' structure = GRanges("chr11", IRanges(2608328, 2699994), strand = '-'),
#' lncRNAID='URS0000759CF4')
#' lengthProduct(gene2)
setMethod("lengthProduct", signature = "LncRNAGene", function(x) {
  return(sum(width(x@structure)))
})


#' Calculate the length product for a micro-RNA gene.
#'
#' This function calculates the length product for a MiRNAGene object,
#' which is the number of characters in the RNA sequence of the micro-RNA.
#'
#' @param x The MiRNAGene object.
#' @return The length product of the micro-RNA gene.
#' @import microRNA
#' @import utils
#'
#' @examples
#' \dontrun{
#' gene3 <- MiRNAGene(id=406988, symbol= 'MIR205',name='MICRORNA 205',
#' description = 'non-coding RNAs that are involved in post-transcriptional
#' regulation of gene expression by affecting both the stability and translation
#' of mRNAs.', structure = GRanges("chr1", IRanges (209432133, 209432242),
#' strand = '+'), final_mature_miRNAID= 'hsa-miR-205')
#' data(hsSeqs)
#' lengthProduct(gene3)
#' }
#' @export
setMethod("lengthProduct", signature = "MiRNAGene", function(x) {
  data(hsSeqs)
  return(nchar(hsSeqs[x@final_mature_miRNAID]))
})

#' Calculate the length product for a tRNA gene.
#'
#' This function calculates the length product for a tRNAGene object,
#' which is the sum of the widths of the gene structure.
#'
#' @param x The tRNAGene object.
#' @return The length product of the tRNA gene.
#' @export
#' @examples gene5 <- tRNAGene(id=100189346, symbol= 'TRA-AGC1-1',name='tRNA-Ala
#' (anticodon AGC) 1-1', description = 'tRNA-Ala (anticodon AGC) 1-1',
#' structure = GRanges("chr6", IRanges(28795964, 28796035), strand = "-"),
#' tRNAID='URS000063E4FD_9606', anticodon='AGC')
#' lengthProduct(gene5)
setMethod("lengthProduct", signature = "tRNAGene", function(x) {
  return(sum(width(x@structure)))
})


#' Calculate the length product for an rRNA gene.
#'
#' This function calculates the length product for an rRNAGene object,
#' which is the sum of the widths of the gene structure.
#'
#' @param x The rRNAGene object.
#' @return The length product of the rRNA gene.
#' @export
#' @examples gene4 <- rRNAGene(id=100008588, symbol= 'RNA18SN5', name='RNA,
#' 18S ribosomal N5',description = 'the 45S rDNA repeat unit encodes a 45S rRNA
#' precursor, transcribed by RNA polymerase I', structure = GRanges("chr21",
#' IRanges(8436876, 8438744), strand = "+"), rRNAID='URS0000726FAB')
#' lengthProduct(gene4)
setMethod("lengthProduct", signature = "rRNAGene", function(x) {
  return(sum(width(x@structure)))
})


#' Calculate the length product for an SNRNA gene.
#'
#' This function calculates the length product for an SNRNAGene object,
#' which is the sum of the widths of the gene structure.
#'
#' @param x The SNRNAGene object.
#' @return The length product of the SNRNA gene.
#' @export
#' @examples gene6 <- SNRNAGene(id=126871, symbol= 'RNU1-1', name='RNU1-1',
#' description = 'RNA, U1 small nuclear 1', structure = GRanges("chr1",
#' IRanges(16514122,16514285), strand = "-"), SNRNAID='URS000032B6B6')
#' lengthProduct(gene6)
setMethod("lengthProduct", signature = "SNRNAGene", function(x) {
  return(sum(width(x@structure)))
})

