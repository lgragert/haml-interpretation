# haml-interpretation

HLA Antibody Markup Language is a data formats standard for reporting raw data from HLA antibody assays and expert or algorithmic interpretation.

This repository demonstrates how HAML files that contain only raw MFI data for an SAB panel and the default classification from the vendor software can be modified to add a second per-bead interpretation `bead-classification` and a per-assay `interpretation` block with a list of HLA antigens for reporting to allocation systems.

This script demonstrates a very simple interpretation scheme using MFI thresholds for bead classifications (greater than 2000 MFI for positive, between 1000 and 2000 MFI for borderline/questionable, and less than 1000 MFI for negative). The name of the per-assay `interpretation-software` and per-bead `classification-entity` added by this script are indicated as `Python-Demo` version `0.0.1`.

The interpretation block uses the `bead-classification` assignments to make phenotype list strings (HLA-plstring) for positive, questionable, and negative categories. The PL string specificities are then converted into OPTN antigen categories in antibody list strings (HLA-abstring) for reporting to the OPTN organ allocation system. The file `antigen_conversion_table_with_rules.csv` contains the mappings between WHO allele names and OPTN antigen categories.

Usage example (default input file if not specified is `HAML.xml`)

```
python3 haml_interpretation.py OLSampleCI2025-05-13_haml.xml
```

Outputs XML file with "_interpretation" appended to filename prefix:

```
OLSampleCI2025-05-13_haml_interpreted.xml
```


More sophisticated antibody analysis techniques involving pattern analysis or identification of false positives could be substituted for the simple interpretation logic in this demo. 
- Class II antigen specificities are based on only the DQB1 and DPB1 categories instead of the heterodimer in this demo.
- Allele-specific antigen categories aren't identified in this demo.
- Candidate HLA typing data is not used as input in this demo, but could in the future

The HAML schema is maintained in this GitHub repository of the Society of Immune Polymorphism:

https://github.com/immunomath/haml

HAML converters from CSV exports of antibody assay vendor are available at:

https://github.com/IHIW/Converters/tree/master/HAMLConverterPy
