# GeneTypes
GeneTypes R package for gene objects which represent different gene types. 

Based on a set of S4 classes for different gene types which are inherited from a virtual gene class that represents genes in general.
Each gene should contain information about it’s ID, HUGO symbol, gene name, description, gene structure (chromosome, start, end, strand, exons, and gene product(s).
It has a class-specific function "lengthProduct(gene)" which returns the length of the product of a specific gene object. 
