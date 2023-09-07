library(GenomicRanges)
library(methods)

### Class object input check
test_that("Gene validity", {
  # Valid gene object
  gene <- new("ProteinCodingGene", id = 1234, symbol = "GENE1",
              name = "Gene 1", description = "Description",
              structure = GRanges(), proteinID = "PROT1")
  expect_true(validObject(gene))

  # Invalid gene object with non-unique id and symbol
  expect_error(gene <- new("MiRNAGene", id = c(111, 234), symbol = "GENE2",
                           name = "Gene 2", description = "Description",
                           structure = GRanges(),
                           final_mature_miRNAID = "miR111"))

  # Invalid gene object with missing structure information
  expect_error(gene <- new("rRNAGene", id = 18, symbol = "GENE18",
                           name = "Gene 18", description = "Description",
                           structure = GRanges(seqnames = "chr1"),
                           rRNAID = "rRNA18S"))

  # Invalid tRNAGene object with non-ribo-nucleotide anticodon
  expect_error(gene <- new("tRNAGene", id = 4567, symbol = "GENE5",
                          name = "Gene 5", description = "Description",
                          structure = GRanges(), tRNAID = "tRNA5",
                          anticodon = "ABC"))
})



### ProteinCodingGene class
test_that("Is an instance of ProteinCodingGene class", {
  gene1 <- new("ProteinCodingGene", id = 1234, symbol = "GENE1",
              name = "Gene 1", description = "Description",
              structure = GRanges(), proteinID = "PROT1")
  expect_s4_class(gene1, "ProteinCodingGene")
})

### LncRNAGene class
test_that("Is an instance of LncRNAGene class", {
  gene2 <- new("LncRNAGene", id = 2345, symbol = "GENE2",
                    name = "Gene 2", description = "Description",
                    structure = GRanges(), lncRNAID = "lncRNA2")
  expect_s4_class(gene2, "LncRNAGene")
})

# MiRNAGene class
test_that("Is an instance of MiRNAGene class", {
  gene3 <- new("MiRNAGene", id = 3456, symbol = "GENE3",
                   name = "Gene 3", description = "Description",
                   structure = GRanges(), final_mature_miRNAID = "miR333")
  expect_s4_class(gene3, "MiRNAGene")
})

# rRNAGene class
test_that("Is an instance of rRNAGene class", {
  gene4 <- new("rRNAGene", id = 18, symbol = "GENE18",
                  name = "Gene 18", description = "Description",
                  structure = GRanges(), rRNAID = "rRNA18S")
  expect_s4_class(gene4, "rRNAGene")
})

# tRNAGene class tests
test_that("Is an instance of tRNAGene class", {
  gene5 <- new("tRNAGene", id = 4567, symbol = "GENE5",
               name = "Gene 5", description = "Description",
               structure = GRanges(), tRNAID = "tRNA5", anticodon = "ACG")
  expect_s4_class(gene5, "tRNAGene")
})

# SNRNAGene class tests
test_that("Is an instance of SNRNAGene class", {
  gene6 <- new("SNRNAGene", id = 6789, symbol = "GENE6",
                   name = "Gene 6", description = "Description",
                   structure = GRanges(), SNRNAID = "SNRNA6")
  expect_s4_class(gene6, "SNRNAGene")
})



### Accessor functions check

# Symbol
test_that("Symbol accessor function works correctly", {
  gene1 <- new("ProteinCodingGene", id = 1234, symbol = "GENE1",
               name = "Gene 1", description = "Description",
               structure = GRanges(), proteinID = "PROT1")
  symbol_value <- symbol(gene1)
  expect_equal(symbol_value, "GENE1")
})

# Name
test_that("Name accessor function works correctly", {
  gene1 <- new("ProteinCodingGene", id = 1234, symbol = "GENE1",
               name = "Gene 1", description = "Description",
               structure = GRanges(), proteinID = "PROT1")
  name_value <- name(gene1)
  expect_equal(name_value, "Gene 1")
})

# Description
test_that("Description accessor function works correctly", {
  gene1 <- new("ProteinCodingGene", id = 1234, symbol = "GENE1",
               name = "Gene 1", description = "Description",
               structure = GRanges(), proteinID = "PROT1")
  description_value <- description(gene1)
  expect_equal(description_value, "Description")
})

# Structure
test_that("Structure accessor function works correctly", {
  gene1 <- new("ProteinCodingGene", id = 1234, symbol = "GENE1",
               name = "Gene 1", description = "Description",
               structure = GRanges("chr7", IRanges(121348878, 121396364),
                                   strand = '-'), proteinID = "PROT1")
  structure_value <- structure(gene1)
  expected_value <- GRanges("chr7", IRanges(121348878, 121396364),
                                                      strand = '-')
  expect_equal(structure_value, expected_value)
})


# Anticodon
test_that("Anticodon accessor function works correctly", {
  gene5 <- new("tRNAGene", id = 4567, symbol = "GENE5",
              name = "Gene 5", description = "Description",
              structure = GRanges(), tRNAID = "tRNA5",
              anticodon = "GCA")
  anticodon_value <- anticodon(gene5)
  expected_output <- "GCA"
  expect_equal(anticodon_value, expected_output)
})

# ProteinID
test_that("ProteinID accessor function works correctly", {
  gene1 <- new("ProteinCodingGene", id = 1234, symbol = "GENE1",
               name = "Gene 1", description = "Description",
               structure = GRanges(), proteinID = "PROT1")
  proteinID_value <- proteinID(gene1)
  expect_equal(proteinID_value, "PROT1")
})

# lncRNAID
test_that("lncRNAID accessor function works correctly", {
  gene2 <- new("LncRNAGene", id = 2345, symbol = "GENE2",
               name = "Gene 2", description = "Description",
               structure = GRanges(), lncRNAID = "lncRNA2")
  lncRNAID_value <- lncRNAID(gene2)
  expect_equal(lncRNAID_value, "lncRNA2")
})

# final_mature_miRNAID
test_that("final_mature_miRNAID accessor function works correctly", {
  gene3 <- new("MiRNAGene", id = 3456, symbol = "GENE3",
               name = "Gene 3", description = "Description",
               structure = GRanges(), final_mature_miRNAID = "miR333")
  final_mature_miRNAID_value <- final_mature_miRNAID(gene3)
  expect_equal(final_mature_miRNAID_value, "miR333")
})

# rRNAID
test_that("rRNAID accessor function works correctly", {
  gene4 <- new("rRNAGene", id = 18, symbol = "GENE18",
               name = "Gene 18", description = "Description",
               structure = GRanges(), rRNAID = "rRNA18S")
  rRNAID_value <- rRNAID(gene4)
  expect_equal(rRNAID_value, "rRNA18S")
})

# tRNAID
test_that("tRNAID accessor function works correctly", {
  gene5 <- new("tRNAGene", id = 4567, symbol = "GENE5",
               name = "Gene 5", description = "Description",
               structure = GRanges(), tRNAID = "tRNA5", anticodon = "ACG")
  tRNAID_value <- tRNAID(gene5)
  expect_equal(tRNAID_value, "tRNA5")
})

# SNRNAID
test_that("SNRNAID accessor function works correctly", {
  gene6 <- new("SNRNAGene", id = 6789, symbol = "GENE6",
               name = "Gene 6", description = "Description",
               structure = GRanges(), SNRNAID = "SNRNA6")
  SNRNAID_value <- SNRNAID(gene6)
  expect_equal(SNRNAID_value, "SNRNA6")
})
