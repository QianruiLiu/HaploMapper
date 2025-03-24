# Geography Visualization (Version 1)

## Overview
This initial version of the geography visualization script (`geography.py`) provides a basic map presentation of haplogroup distributions. It displays only the main haplogroups without analyzing or visualizing any subtypes.

## Usage
- **Input Files:**  
  - `mtDNA_haplogroup_frequencies.tsv`  
  - `Y_haplogroup_frequencies.tsv`
- **Output:** A simple HTML file displaying the geographic distribution of main haplogroups.

## How to Run
Execute the script in your terminal:
```bash
python geography.py --y_input Y_haplogroup_frequencies.tsv --mt_input mtDNA_haplogroup_frequencies.tsv
