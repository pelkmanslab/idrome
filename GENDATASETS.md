# Generated datasets

There are several "datasets" that are generated from the original datasets detailed [here](RAWDATASETS.md). These are also stored as CSVs and can be used independently from the code.

## "polyprots" viral dataset

This dataset is identical to the original raw viral dataset from Uniprot, but any protein that had CHAIN information indicating that they are cleaved during their maturation process were split into separate proteins. The LCR and IDR analysis was re-run on the full dataset. The columns are unchanged from the raw viral dataset.

## "iddrep" datasets

These are IDR representation datasets, where the separate IDRs are represented as some abstraction (currently as 20D frequency vectors). Proteins containing multiple IDRs have multiple rows with each separate IDR having its own row. 

### Columns

- ACCID: This is Uniprot accession ID of the protein containing this IDR with `_dX` appended to it. The `X` is incrementing numerical value that differentiates IDRs located on the same protein. 
- SEQ: the sequence of the IDR
- LENGTH: the length of the IDR
- A,C,D...Y: the frequency of a given amino acid in this IDR
- PARENT: the Uniprot accession ID of the protein containing this IDR. 

## "sgt" dataset

The stress granule targeting dataset. 
