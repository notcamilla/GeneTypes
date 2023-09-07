library(microRNA)
data(hsSeqs)
library(methods)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

### Micro-RNA Gene class functions

# getseedseq
test_that("getseedseq function returns expected output", {
  gene3 <- new("MiRNAGene", id = 3456, symbol = "GENE3",
               name = "Gene 3", description = "Description",
               structure = GRanges(), final_mature_miRNAID = "hsa-miR-1234")
  result <- getseedseq(gene3)
  expect_is(result, "character")
  expected_content <- seedRegions(hsSeqs["hsa-miR-1234"])
  expect_identical(result, expected_content)
})

# matureMiRNAseq
test_that("matureMiRNAseq function returns expected output", {
  gene3 <- new("MiRNAGene", id = 3456, symbol = "GENE3",
               name = "Gene 3", description = "Description",
               structure = GRanges(), final_mature_miRNAID = "hsa-miR-1234")
  result <- matureMiRNAseq(gene3)
  expect_is(result, "RNAString")
  expected_content <- RNAString(hsSeqs["hsa-miR-1234"])
  expect_identical(result, expected_content)
})

### proteinseq
test_that("proteinseq method returns expected output ", {
  gene1 <- new("ProteinCodingGene", id = 1234, symbol = "GENE1",
               name = "Gene 1", description = "Description",
               structure = GRanges("chr1", IRanges(16514122, 16514285),
                                   strand = "-"), proteinID = "PROT1")
  result <- proteinseq(gene1)
  expect_is(result, "list")

})


### LncRNAseq
test_that("LncRNAseq function returns expected output", {
  gene2 <- new("LncRNAGene", id = 2345, symbol = "GENE2",
               name = "Gene 2", description = "Description",
               structure = GRanges("chr11", IRanges(2608328, 2699994),
               strand = '-'), lncRNAID = "lncRNA2")
  result <- LncRNAseq(gene2)
  expect_is(result, "RNAString")
})

### tRNAseq
test_that("tRNAseq function returns expected output", {
  gene5 <- new("tRNAGene", id = 4567, symbol = "GENE5",
               name = "Gene 5", description = "Description",
               structure = GRanges("chr6", IRanges(28795964, 28796035),
               strand = "-"), tRNAID = "tRNA5", anticodon = "ACG")
  result <- tRNAseq(gene5)
  expect_is(result, "RNAString")
})

### rRNAseq
test_that("rRNAseq function returns expected output", {
  gene4 <- new("rRNAGene", id = 18, symbol = "GENE18",
               name = "Gene 18", description = "Description",
               structure = GRanges("chr21", IRanges(8436876, 8438744),
               strand = "+"), rRNAID = "rRNA18S")
  result <- rRNAseq(gene4)
  expect_is(result, "RNAString")
})

### SNRNAseq
test_that("SNRNAseq function returns expected output", {
  gene6 <- new("SNRNAGene", id = 6789, symbol = "GENE6",
               name = "Gene 6", description = "Description",
               structure = GRanges("chr1", IRanges(16514122, 16514285),
               strand = "-"), SNRNAID = "SNRNA6")
  result <- SNRNAseq(gene6)
  expect_is(result, "RNAString")
})

### mRNAseq
test_that("mRNAseq function returns expected output", {
  gene1 <- new("ProteinCodingGene", id = 1234, symbol = "GENE1",
               name = "Gene 1", description = "Description",
               structure = GRanges("chr1", IRanges(16514122, 16514285),
                                   strand = "-"), proteinID = "PROT1")
  result <- mRNAseq(gene1)
  expect_is(result, "RNAStringSetList")
})


### Product length function

# ProteinCodingGene Class
test_that("lengthProduct returns the correct total length for
          ProteinCodingGene", {
  gene1 <- new("ProteinCodingGene", id = 1234, symbol = "GENE1",
           name = "Gene 1", description = "Description",
           structure = GRanges("chr1", IRanges(16514122, 16514285),
           strand = "-"), proteinID = "PROT1")
  total_length <- lengthProduct(gene1)
  protein_seqs <- proteinseq(gene1)
  aa_length <- vapply(protein_seqs, function(seq) sum(nchar(seq)), integer(1))
  expect_is(total_length, "integer")
  expect_equal(total_length, aa_length)
})

# LncRNAGene Class
test_that("lengthProduct returns the correct total length for LncRNAGene", {
  gene2 <- new("LncRNAGene", id = 1, symbol = "GENE2", name = "Gene 2",
              description = "Description",structure = GRanges(seqnames = "chr1",
              ranges = IRanges(start = c(100, 200, 300), end = c(150, 250, 350))),
              lncRNAID = "lncRNA2")
  total_length <- lengthProduct(gene2)
  expect_is(total_length, "integer")
  expect_equal(total_length, sum(width(gene2@structure)))
})


# MiRNAGene Class
test_that("lengthProduct method returns expected length product for MiRNAGene", {
  gene3 <- new("MiRNAGene", id = 3456, symbol = "GENE3",
               name = "Gene 3", description = "Description",
               structure = GRanges(), final_mature_miRNAID = "hsa-miR-1234")
  total_length <- lengthProduct(gene3)
  expect_is(total_length, "integer")
  expect_equal(total_length, nchar(hsSeqs["hsa-miR-1234"]))
})


# rRNAGene Class
test_that("lengthProduct returns the correct total length for rRNAGene", {
  gene4 <- new("rRNAGene", id = 18, symbol = "GENE18",name = "Gene 18",
               description = "Description",structure = GRanges(seqnames = "chr1",
               ranges = IRanges(start = c(100, 200, 300), end = c(150, 250, 350))),
               rRNAID = "rRNA18S")
  total_length <- lengthProduct(gene4)
  expect_is(total_length, "integer")
  expect_equal(total_length, sum(width(gene4@structure)))
})


# tRNAGene Class
test_that("lengthProduct returns the correct total length for tRNAGene", {
  gene5 <- new("tRNAGene", id = 4567, symbol = "GENE5", name = "Gene 5",
               description = "Description",structure = GRanges(seqnames = "chr1",
               ranges = IRanges(start = c(100, 200, 300), end = c(150, 250, 350))),
               tRNAID = "tRNA5", anticodon = 'ACG')
  total_length <- lengthProduct(gene5)
  expect_is(total_length, "integer")
  expect_equal(total_length, sum(width(gene5@structure)))
})


# SNRNAGene Class
test_that("lengthProduct returns the correct total length for SNRNAGene", {
  gene6 <- new("SNRNAGene", id = 6789, symbol = "GENE6", name = "Gene 6",
               description = "Description",structure = GRanges(seqnames = "chr1",
               ranges = IRanges(start = c(100, 200, 300), end = c(150, 250, 350))),
               SNRNAID = "SNRNA6")
  total_length <- lengthProduct(gene6)
  expect_is(total_length, "integer")
  expect_equal(total_length, sum(width(gene6@structure)))
})


