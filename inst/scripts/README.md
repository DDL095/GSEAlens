# Data Generation for inst/extdata

This directory contains scripts and documentation describing how the data in `inst/extdata/` was generated.

## pathway_annotations.csv

**Source**: Manually curated pathway annotation file for demonstration and custom pathway visualization purposes.

**Columns**:
- `ID`: Pathway identifier (e.g., HALLMARK_HYPOXIA, REACTOME_EXAMPLE)
- `Brief_Description_CN`: Brief description in Simplified Chinese
- `Brief_Description_EN`: Brief description in English
- `Abstract`: Detailed pathway abstract text
- `Custom_Note`: User-defined custom notes for the pathway

**Licensing**: This file is provided as example data for GSEAlens package functionality demonstration.

**Generation Method**:
1. Pathway IDs were selected from commonly used gene set collections (MSigDB Hallmark, Reactome, KEGG)
2. Descriptions were curated based on official pathway documentation
3. The file can be extended by users to add custom pathway annotations for their own gene sets

**Usage in GSEAlens**:
This annotation file can be loaded via the `create_addition_data()` function to provide custom pathway descriptions in the Shiny application.
