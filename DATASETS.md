# Raw Datasets 

There are two major datasets used in this research. They can be found in the `human/` and `viral/` folders. All the important data files are stored as CSVs. These CSVs are comma-delimited. The rows are separate proteins while the columns contain data pertinent to each specific protein.

## Common

The following attributes are common the both raw human and viral datasets.

- ACCID: the protein's Uniprot Accession ID
- ORGID: the organism's ID in Uniprot
- PROTNAME: the protein's names and synonyms
- LCR: the protein's low complexity regions as predicted by CAST. The data are stored as follows: `LCR amino acid single letter code :  start position _ end position (inclusive) $ CAST score [positions of residues belonging to LCR];` The data can be extracted easily using the following regex: `([A-Z]):([\d]+).*?_([\d]+).*?\$([\S]*)@(\[.*?\]);`, which can be used like so: `import re; re.findall(r'([A-Z]):([\d]+).*?_([\d]+).*?\$([\S]*)@(\[.*?\]);', "target")` 
- GO: the protein's gene ontology terms separated by `; `
- GOID: the protein's gene ontology IDs separated by `; `
- TAXON: the organism's taxonomic information, stored in Uniprot style
- CHAIN: descriptions of the mature polypeptides following processing
- DOMFT: sequence position-dependent annotation (features)
- DOMCC: sequence position-independent annotation
- REGION: region of interest in the sequence
- COMPBIAS: compositional bias in the sequence as reported by Uniprot
- SEQ: the raw sequence
- KEYWORDS: associated Uniprot keywords
- GENENAME: name(s) of the gene(s) encoding a protein
- IDRXX: intrinsically disordered regions that are at least as long as the value indicated by XX (e.g. IDR10 includes all IDRs that are 10 AAs or longer) as predicted by IUPred. The data are stored as follows: `start : end (inclusive) ;`. The necessary regex is: `([\d]+)_([\d]+)`

## Specific to Viral Dataset

Due to the high prevalence of polyproteins (esp. in positive ssRNA viruses), the viral proteins were split into their mature forms using the information in the CHAIN column. The specific implementation can be viewed in `DATA_PREP.ipynb`
