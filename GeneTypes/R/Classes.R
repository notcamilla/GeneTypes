#' A general class for Genes with various attributes
#'
#' @title Gene
#' @slot id Numeric vector specifying the unique identifier for the gene.
#' @slot symbol Character vector specifying the symbol associated with the gene.
#' @slot name Character vector specifying the name of the gene.
#' @slot description Character vector providing additional description of gene.
#' @slot structure GRanges object representing the genomic structure of gene.
setClass("Gene", representation(id="numeric", symbol="character",
                                name="character", description="character",
                                structure="GRanges", "VIRTUAL") ,
         validity = function(object) {
           if (length(object@id) > 1 || (length(object@symbol) > 1 )){
             return("@id and @symbol must be unique")
           } else if (!is(object@structure, "GRanges") ||
                      is.null(seqnames(object@structure)) ||
                      is.null(start(object@structure)) ||
                      is.null(end(object@structure)) ||
                      is.null(strand(object@structure))) {
             return("@structure must be a valid GRanges object")
           } else { return(TRUE) }
           }
         )


#' A subclass of gene representing protein-coding genes.
#' @title ProteinCodingGene
#' @slot proteinID The Uniprot ID, character.
#' @importFrom methods new
#' @export
#' @examples
#' gene1 <- ProteinCodingGene(id=10447, symbol= 'FAM3C', name='FAM3 metabolism
#' regulating signaling molecule C',description = 'encodes a secreted protein
#' with a GG domain', structure = GRanges("chr7",IRanges(121348878, 121396364),
#' strand = '-'), proteinID = 'Q92520')
ProteinCodingGene <- setClass("ProteinCodingGene", contains="Gene",
                       slots= list(proteinID="character"))

#' A subclass of gene representing long non-coding genes.
#' @title LncRNAGene
#' @slot lncRNAID The RNAcentral ID, character.
#' @importFrom methods new
#' @export
#' @examples
#' gene2 <- LncRNAGene(id=10984, symbol= 'KCNQ1OT1', name='antisense transcript
#' 1', description = 'interacts with chromatin and regulates transcription of
#' multiple target genes through epigenetic modifications', structure =
#' GRanges("chr11", IRanges(2608328, 2699994), strand = '-'),
#' lncRNAID='URS0000759CF4')
LncRNAGene <- setClass("LncRNAGene", contains="Gene",
                       slots=list(lncRNAID='character'))

#' A subclass of gene representing microRNA genes.
#' @title MiRNAGene
#' @slot final_mature_miRNAID The miRBase Sanger ID, character.
#' @importFrom methods new
#' @export
#' @examples
#' gene3 <- MiRNAGene(id=406988, symbol= 'MIR205', name='MICRORNA 205',
#' description = 'non-coding RNAs that are involved in post-transcriptional
#' regulation of gene expression by affecting both the stability and translation
#' of mRNAs.', structure = GRanges("chr1", IRanges (209432133, 209432242),
#' strand = '+'), final_mature_miRNAID= 'hsa-miR-205')
MiRNAGene <- setClass("MiRNAGene", contains="Gene",
                      slots= list(final_mature_miRNAID= 'character'))

#' A subclass of gene representing ribosomal RNA genes.
#' @title rRNAGene
#' @slot rRNAID The RNAcentral ID, character.
#' @importFrom methods new
#' @export
#' @examples
#' gene4 <- rRNAGene(id=100008588, symbol= 'RNA18SN5', name='RNA, 18S ribosomal
#' N5',description = 'the 45S rDNA repeat unit encodes a 45S rRNA precursor,
#' transcribed by RNA polymerase I', structure = GRanges("chr21",
#' IRanges(8436876, 8438744), strand = "+"), rRNAID='URS0000726FAB')
rRNAGene <- setClass("rRNAGene", contains="Gene",
                     slots = list(rRNAID='character'))

#' A subclass of gene representing transfer RNA genes.
#' @title tRNAGene
#' @slot tRNAID The RNAcentral ID, character.
#' @slot tRNAID A character vector specifying the tRNA identifier.
#' @slot anticodon A character vector specifying the anticodon sequence of tRNA.
#' @param object An instance of the tRNAGene class to be validated.
#' @return TRUE if the object is valid, an error message otherwise.
#' @details This validity function checks if the anticodon slot contains only
#' ribo-nucleotides and if it is 3 nucleotides long.
#' If the anticodon fails any of these checks, an appropriate error message is
#' returned; otherwise, TRUE is returned.
#' @importFrom methods new
#' @export
#' @examples
#' gene5 <- tRNAGene(id=100189346, symbol= 'TRA-AGC1-1', name='tRNA-Ala
#' (anticodon AGC) 1-1', description = 'tRNA-Ala (anticodon AGC) 1-1',
#' structure = GRanges("chr6", IRanges(28795964, 28796035), strand = "-"),
#' tRNAID='URS000063E4FD_9606', anticodon='AGC')
tRNAGene <- setClass("tRNAGene", contains="Gene",
                     slots = list(tRNAID='character', anticodon='character'),
                     validity = function(object) {
                       if (!grepl("^[ACGU]+$", object@anticodon)) {
                        return("@anticodon must be contain only ribo-nucleotides")
                       } else if (nchar(object@anticodon) != 3) {
                         return('@anticodom must be 3 ribo-nucleotides long')
                       } else { return(TRUE) }
                       })


#' A subclass of gene representing small-nuclear RNA genes.
#' @title SNRNAGene
#' @slot SNRNAID The RNAcentral ID, character.
#' @importFrom methods new
#' @export
#' @examples
#' gene6 <- SNRNAGene(id=126871, symbol= 'RNU1-1', name='RNU1-1', description =
#' 'RNA, U1 small nuclear 1', structure = GRanges("chr1", IRanges(16514122,
#' 16514285), strand = "-"), SNRNAID='URS000032B6B6')
SNRNAGene <- setClass("SNRNAGene", contains="Gene",
                      slots = list(SNRNAID='character'))


### ACCESSOR AND REPLACEMENT FUNCTIONS ###
#' Accessor function for the symbol of a gene object.
#'
#' @param x An object of class Gene.
#'
#' @return The symbol of the Gene object.
#' @export
setGeneric("symbol", function(x) standardGeneric("symbol"))

#' Accessor function for the symbol of a Gene object.
#'
#' @rdname symbol
#' @param x A Gene object.
#'
#' @return The updated Gene object.
#'
#' @examples gene1 <- ProteinCodingGene(id=10447, symbol= 'FAM3C', name='FAM3
#' metabolism regulating signaling molecule C',description = 'encodes a secreted
#' protein with a GG domain', structure = GRanges("chr7",IRanges(121348878,
#' 121396364), strand = '-'), proteinID = 'Q92520')
#' symbol(gene1)
#' @export
setMethod("symbol", "Gene", function(x) x@symbol)

#' Accessor function for the name of a Gene object.
#'
#' @param x A Gene object.
#' @return The name of the Gene.
#' @export
setGeneric("name", function(x) standardGeneric("name"))


#' Accessor function for the name of a Gene object.
#'
#' @rdname name
#' @param x A Gene object.
#'
#' @return The updated Gene object.
#'
#' @examples gene1 <- ProteinCodingGene(id=10447, symbol= 'FAM3C', name='FAM3
#' metabolism regulating signaling molecule C',description = 'encodes a secreted
#' protein with a GG domain', structure = GRanges("chr7",IRanges(121348878,
#' 121396364), strand = '-'), proteinID = 'Q92520')
#' name(gene1)
#' @export
setMethod("name", "Gene", function(x) x@name)

#' Accessor function for the description of a Gene object.
#'
#' @param x A Gene object
#' @return The description of the Gene
#' @export
setGeneric("description", function(x) standardGeneric("description"))


#' Accessor function for the description of a Gene object.
#'
#' @rdname description
#' @param x A Gene object.
#'
#' @return The updated Gene object.
#'
#' @examples gene1 <- ProteinCodingGene(id=10447, symbol= 'FAM3C', name='FAM3
#' metabolism regulating signaling molecule C',description = 'encodes a secreted
#' protein with a GG domain', structure = GRanges("chr7",IRanges(121348878,
#' 121396364), strand = '-'), proteinID = 'Q92520')
#' description(gene1)
#' @export
setMethod("description", "Gene", function(x) x@description)


#' Accessor function for the structure of a Gene object.
#'
#' @param x A Gene object
#' @return The structure of the Gene
#' @export
setGeneric("structure", function(x) standardGeneric("structure"))


#' Accessor function for the structure of a Gene object..
#'
#' @rdname structure
#' @param x A Gene object.
#'
#' @return The updated Gene object.
#'
#' @examples gene1 <- ProteinCodingGene(id=10447, symbol= 'FAM3C', name='FAM3
#' metabolism regulating signaling molecule C',description = 'encodes a secreted
#' protein with a GG domain', structure = GRanges("chr7",IRanges(121348878,
#' 121396364), strand = '-'), proteinID = 'Q92520')
#' structure(gene1)
#' @export
setMethod("structure", "Gene", function(x) (x@structure))

#' Accessor function for the anticodon of a transfer RNAGene object.
#'
#' @param x A transfer RNAGene object
#' @return The anticodon of the transfer RNAGene
#' @export
setGeneric("anticodon", function(x) standardGeneric("anticodon"))


#' Accessor function for the anticodon of a transfer RNA Gene object.
#'
#' @rdname anticodon
#' @param x A Gene object.
#'
#' @return The updated Gene object.
#'
#' @examples gene5 <- tRNAGene(id=100189346, symbol= 'TRA-AGC1-1',name='tRNA-Ala
#' (anticodon AGC) 1-1', description = 'tRNA-Ala (anticodon AGC) 1-1',
#' structure = GRanges("chr6", IRanges(28795964, 28796035), strand = "-"),
#' tRNAID='URS000063E4FD_9606', anticodon='AGC')
#' anticodon(gene5)
#' @export
setMethod("anticodon", "tRNAGene", function(x) x@anticodon)


#' Accessor function for the proteinID of a Coding Gene object.
#'
#' @param x A Coding Gene object
#' @return The proteinID of the Coding Gene
#' @export
setGeneric("proteinID", function(x) standardGeneric("proteinID"))


#' Accessor function for the proteinID of a Coding Gene object.
#'
#' @rdname proteinID
#' @param x A Gene object.
#'
#' @return The updated Gene object.
#'
#' @examples gene1 <- ProteinCodingGene(id=10447, symbol= 'FAM3C', name='FAM3
#' metabolism regulating signaling molecule C',description = 'encodes a secreted
#' protein with a GG domain', structure = GRanges("chr7",IRanges(121348878,
#' 121396364), strand = '-'), proteinID = 'Q92520')
#' proteinID(gene1)
#' @export
setMethod("proteinID", "ProteinCodingGene", function(x) x@proteinID)


#' Accessor function for the lncRNAID of a long non-coding RNA Gene object.
#'
#' @param x A LncRNAGene object
#' @return The lncRNAID of the LncRNAGene
#' @export
setGeneric("lncRNAID", function(x) standardGeneric("lncRNAID"))


#' Accessor function for the lncRNAID of a long non-coding RNA object.
#'
#' @rdname lncRNAID
#' @param x A Gene object.
#'
#' @return The updated Gene object.
#'
#' @examples gene2 <- LncRNAGene(id=10984, symbol= 'KCNQ1OT1', name='antisense
#' transcript 1', description = 'interacts with chromatin and regulates
#' transcription of multiple target genes through epigenetic modifications',
#' structure = GRanges("chr11", IRanges(2608328, 2699994), strand = '-'),
#' lncRNAID='URS0000759CF4')
#' lncRNAID(gene2)
#' @export
setMethod("lncRNAID", "LncRNAGene", function(x) x@lncRNAID)


#' Accessor function for the final_mature_miRNAID of a micro-RNA Gene object.
#'
#' @param x A micro-RNA Gene object.
#' @return The final_mature_miRNAID slot value of the MiRNAGene object.
#'
#' @export
setGeneric("final_mature_miRNAID", function(x)
  standardGeneric("final_mature_miRNAID"))


#' Accessor function for the final_mature_miRNAID of a micro-RNA object.
#'
#' @rdname final_mature_miRNAID
#' @param x A Gene object.
#'
#' @return The updated Gene object.
#'
#' @examples gene3 <- MiRNAGene(id=406988, symbol= 'MIR205', name='MICRORNA 205',
#' description = 'non-coding RNAs that are involved in post-transcriptional
#' regulation of gene expression by affecting both the stability and translation
#' of mRNAs.', structure = GRanges("chr1", IRanges (209432133, 209432242),
#' strand = '+'), final_mature_miRNAID= 'hsa-miR-205')
#' final_mature_miRNAID(gene3)
#' @export
setMethod("final_mature_miRNAID", "MiRNAGene", function(x)
  x@final_mature_miRNAID)


#' Accessor function for the rRNAID of an ribosomial RNA Gene object.
#'
#' @param x An rRNAGene object.
#' @return The rRNAID slot value of the rRNAGene object.
#' @export
setGeneric("rRNAID", function(x) standardGeneric("rRNAID"))


#' Accessor function for the rRNAID of a ribosomial RNA object.
#'
#' @rdname rRNAID
#' @param x A Gene object.
#'
#' @return The updated Gene object.
#'
#' @examples gene4 <- rRNAGene(id=100008588, symbol= 'RNA18SN5', name='RNA,
#' 18S ribosomal N5',description = 'the 45S rDNA repeat unit encodes a 45S rRNA
#' precursor, transcribed by RNA polymerase I', structure = GRanges("chr21",
#' IRanges(8436876, 8438744), strand = "+"), rRNAID='URS0000726FAB')
#' rRNAID(gene4)
#' @export
setMethod("rRNAID", "rRNAGene", function(x) x@rRNAID)


#' Accessor function for the tRNAID of a transfer RNA Gene object.
#'
#' @param x A tRNAGene object.
#' @return The tRNAID slot value of the transfer RNA Gene object.
#' @export
setGeneric("tRNAID", function(x) standardGeneric("tRNAID"))


#' Accessor function for the tRNAID of a transfer RNA object.
#'
#' @rdname tRNAID
#' @param x A Gene object.
#'
#' @return The updated Gene object.
#'
#' @examples gene5 <- tRNAGene(id=100189346, symbol= 'TRA-AGC1-1',name='tRNA-Ala
#' (anticodon AGC) 1-1', description = 'tRNA-Ala (anticodon AGC) 1-1',
#' structure = GRanges("chr6", IRanges(28795964, 28796035), strand = "-"),
#' tRNAID='URS000063E4FD_9606', anticodon='AGC')
#' tRNAID(gene5)
#' @export
setMethod("tRNAID", "tRNAGene", function(x) x@tRNAID)


#' Accessor function for SNRNAID slot in small-nuclear RNA Gene objects.
#'
#' @param x An small-nuclear RNA Gene object.
#'
#' @return The value of the SNRNAID slot.
#'
#' @export
setGeneric("SNRNAID", function(x) standardGeneric("SNRNAID"))


#' Accessor function for the SNRNAID of a small-nuclear RNA object.
#'
#' @rdname SNRNAID
#' @param x A Gene object.
#'
#' @return The updated Gene object.
#' @examples gene6 <- SNRNAGene(id=126871, symbol= 'RNU1-1', name='RNU1-1',
#' description = 'RNA, U1 small nuclear 1', structure = GRanges("chr1",
#' IRanges(16514122,16514285), strand = "-"), SNRNAID='URS000032B6B6')
#' SNRNAID(gene6)
#' @export
setMethod("SNRNAID", "SNRNAGene", function(x) x@SNRNAID)


