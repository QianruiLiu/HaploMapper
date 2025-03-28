# HaploMapper

## Overview
**HaploMapper** is an interactive web application designed for temporal and geographic visualization of mitochondrial DNA (mtDNA) and Y-chromosome haplogroup frequencies from ancient DNA datasets. It enables detailed exploration of haplogroup distributions, migration patterns, and human evolutionary processes by providing dynamic visualizations of both basal haplogroups and their primary subclades.

## Methodology
The analysis workflow in HaploMapper is divided into three main steps:

### Step 1: Data Filtering (`filter_annotation.py`)
- **Purpose:** Filters the annotation data based on a provided list of ancient or modern sample identifiers. Only samples marked "PASS" in quality assessment are retained.
- **Inputs:**
  - `Annotations.xlsx` (annotation data)
  - `Ancient_samples.txt` or `Modern_samples.txt` (sample IDs for filtering)
- **Output:**
  - `annotation_filtered.tsv` (filtered annotation file)

### Step 2: Haplogroup Frequency Calculation (`new_haplogroup_frequency_tables.py`)
- **Purpose:** Calculates frequencies of basal haplogroups (single-letter) and their primary subclades (e.g., A1, B2) from the filtered annotation data.
- **Input:**
  - `annotation_filtered.tsv`
- **Outputs:**
  - `Y_haplogroup_frequencies.tsv` (Y-chromosome frequency table)
  - `mtDNA_haplogroup_frequencies.tsv` (mtDNA frequency table)

### Step 3: Interactive Geographic Visualization (`perfect_geography.py`)
- **Purpose:** Generates an interactive map visualization, displaying double-ring pie charts for each population. The inner ring represents basal haplogroups, while the outer ring represents subclades. Includes interactive filtering by time periods and haplogroup types (mtDNA/Y-chromosome).
- **Inputs:**
  - `Y_haplogroup_frequencies.tsv`
  - `mtDNA_haplogroup_frequencies.tsv`
- **Output:**
  - `haplogroup_map.html` (interactive HTML visualization)

## Files included
```
HaploMapper/
├── data/
│   ├── raw/
│   │   ├── Annotations.xlsx
│   │   └── Ancient_samples.txt
│   └── processed/
│       ├── annotation_filtered.tsv
│       ├── mtDNA_haplogroup_frequencies.tsv
│       └── Y_haplogroup_frequencies.tsv
├── src/
│   ├── step1_filtering/
│   │   └── filter_annotation.py
│   ├── step2_frequency/
│   │   ├── new_haplogroup_frequency_tables.py
│   │   └── haplogroup_frequency_tables.py (deprecated version)
│   └── step3_geography/
│       ├── geography.py (version 1 - basal haplogroups only)
│       ├── newgeography.py (version 2 - includes subclades, alignment issues)
│       ├── lately_geography.py (version 3 - improved alignment)
│       └── perfect_geography.py (version 4 - final version, perfect alignment)
├── outputs/
│   └── haplogroup_map.html
└── requirements.txt
```
**Note: The version 1,2,3 of geographic visualizaion part can be respectively found in branches geography-v1, v2, v3, the old version of haplogroup frequency calculation part can be found in branch geography-v1.***

## Installation
Clone the repository and install required Python dependencies:

The following commands are done in bash:

`git clone git@github.com:QianruiLiu/HaploMapper.git`

`cd HaploMapper`

`pip install -r requirements.txt`

## Dependencies
**Python (recommended 3.12.9)**

**pandas**

**numpy**

**folium**

**openpyxl**

## Running the Code
### Step 1: Filtering
`python src/step1_filtering/filter_annotation.py -a data/raw/Annotations.xlsx -s data/raw/Ancient_samples.txt`

### Step 2: Frequency Tables
`python src/step2_frequencyTable/new_haplogroup_frequency_tables.py --input data/processed/annotation_filtered.tsv`

### Step 3: Interactive Map
`python src/step3_geography/perfect_geography.py --y_input data/processed/Y_haplogroup_frequencies.tsv --mt_input data/processed/mtDNA_haplogroup_frequencies.tsv`

## Known Bugs

Currently, no critical bugs have been identified. Please report any issues via the GitHub issue tracker.

## Frequently Asked Questions (FAQs)
### What if my input file format doesn't match the provided examples?

HaploMapper expects specific column headers. Make sure your files match the expected structure detailed above or adjust your files accordingly.

### Why do I have problems in generating double-ring pie chart?

Make sure you get internet connected and have accession to chart.js: https://cdn.jsdelivr.net/npm/chart.js.

## Contact

### For questions or issues, please contact Qianrui Liu:
Email: qianrui.liu.1767@student.lu.se
