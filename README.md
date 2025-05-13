# haml-interpretation

HLA Antibody Markup Language is a data formats standard for reporting raw data from HLA antibody assays and expert or algorithmic interpretation.

This repository demonstrates how HAML files that contain only raw MFI data for an SAB panel and the default classification from the vendor software can be modified to add a second per-bead interpretation `bead-classification` and a per-assay `interpretation` block with a list of HLA antigens for reporting to allocation systems.

This script demonstrates what a very simple automated antibody analysis program might do with a basic interpretation scheme using MFI threshold settings to make a new set of bead classifications beyond what the vendor software provided (greater than 2000 MFI for positive, between 1000 and 2000 MFI for borderline/questionable, and less than 1000 MFI for negative). The name of the per-assay `interpretation-software` and per-bead `classification-entity` added by this script are indicated as `Python-Demo` version `0.0.1`. The `interpretation-reason` made this script is labeled as `MFI-cutoff-setting` for all beads.

The interpretation block uses the `bead-classification` assignments to make phenotype list strings (HLA-plstring) for positive, questionable, and negative categories. The PL string specificities are then converted into OPTN antigen categories in antibody list strings (HLA-abstring) for reporting to the OPTN organ allocation system. The file `antigen_conversion_table_with_rules.csv` contains the mappings between WHO allele names and OPTN antigen categories.

Usage example (default input file if not specified is `HAML.xml`)

```
python3 haml_interpretation.py OLSampleCI2025-05-13_haml.xml
```

Outputs XML file with "_interpretation" appended to filename prefix:

```
OLSampleCI2025-05-13_haml_interpreted.xml
```

The input HAML file is augmented


More sophisticated antibody analysis techniques involving pattern analysis or identification of false positives could be substituted for the simple `MFI-cutoff-setting` interpretation reason logic in this demo. 
- Class II antigen specificities are based on the DQB1 and DPB1 alleles instead of the heterodimer in this demo. The DQA1 and DPA1 alleles are ignored.
- Allele-specific antigen categories aren't identified in this demo.
- Candidate HLA typing data is not used as input in this demo, but could in the future to rule out self-antigens.
- In addition to bead-level interpretation reasons, assay-level patterns/features that are identified (e.g. Bw4-pattern, no-clear-cutoff) could be also be documented in the interpretation block when this feature is supported in the HAML schema.

The HAML schema is maintained in this GitHub repository of the Society of Immune Polymorphism:

https://github.com/immunomath/haml

HAML converters from CSV exports of antibody assay vendor are available at:

https://github.com/IHIW/Converters/tree/master/HAMLConverterPy
