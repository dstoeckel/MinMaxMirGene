MinMaxMirGene
=============

This is a simple ILP formulation for computing the minimum amount of miRNAs
needed to cover a maximum number of Genes in a biological Pathway.

As an input it needs a file specifying the miRNA to target mappings. The format
is a simple line based, tab-delimited text file:

```
miRNA-ID1	GeneSymbol1
miRNA-ID2	GeneSymbol2
miRNA-ID3	GeneSymbol3
```

As dependencies, an installation of the CPLEX ILP solver is required.

This program is not necessarily restricted to miRNA - Gene interactions, but
can be applied in any scenario, where a number of items needs to be covered by
a minimum amount of "controllers".

This program is licensed under the term of the GNU General Public License,
Version 3 (GPLv3).
