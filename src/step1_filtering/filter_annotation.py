#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functionality:
    Retain the rows in an annotation file (e.g., annotation.xlsx) whose "Genetic ID" column
    matches values found in the second column of a given samples file (ancient or modern).

Usage Examples:
1. Filter with an ancient samples file:
   python filter_annotation_by_samples.py --annotation annotation.xlsx --samples Ancient_samples.txt
   
2. Filter with a modern samples file:
   python filter_annotation_by_samples.py --annotation annotation.xlsx --samples Modern_samples.txt -o modern_filtered.tsv

3. Use short flags:
   python filter_annotation_by_samples.py -a annotation.xlsx -s Ancient_samples.txt [-o output.tsv]
"""

import argparse
import pandas as pd
import sys
import os

def main():
    parser = argparse.ArgumentParser(
        description="Filter the annotation file based on a list of sample names (ancient or modern)."
    )
    parser.add_argument(
        "-a", "--annotation",
        required=True,
        help="Path to the annotation.xlsx file"
    )
    parser.add_argument(
        "-s", "--samples",
        required=True,
        help="Path to the file containing the sample names to retain (e.g., ancient or modern)."
    )
    parser.add_argument(
        "-o", "--output",
        default="annotation_filtered.tsv",
        help="Output file name (default: annotation_filtered.tsv)"
    )
    args = parser.parse_args()
    
    # Check if input files exist
    if not os.path.exists(args.annotation):
        sys.exit(f"Error: Annotation file '{args.annotation}' not found.")
    if not os.path.exists(args.samples):
        sys.exit(f"Error: Samples file '{args.samples}' not found.")
    
    try:
        # 1. Read the samples file, get the second column (list of sample names)
        try:
            df_samples = pd.read_csv(
                args.samples,
                sep=r"\s+",
                header=None,
                names=["SampleIndex", "SampleName"],
                usecols=[0, 1]
            )
        except pd.errors.EmptyDataError:
            sys.exit(f"Error: Sample file '{args.samples}' is empty.")
        except pd.errors.ParserError:
            sys.exit(f"Error: Could not parse '{args.samples}'. Please check the file format.")
        
        # Convert to a set for efficient membership checks
        sample_name_set = set(df_samples["SampleName"].unique())
        if not sample_name_set:
            sys.exit(f"Error: No valid sample names found in '{args.samples}'.")
        print(f"Found {len(sample_name_set)} distinct sample names in {args.samples}.")

        # 2. Read annotation.xlsx
        try:
            df_annotation = pd.read_excel(args.annotation)
        except ValueError as e:
            sys.exit(f"Error reading Excel file: {str(e)}")
        except Exception as e:
            sys.exit(f"Error reading annotation file: {str(e)}")
            
        print(f"Loaded {len(df_annotation)} rows from {args.annotation}.")

        # Ensure that "Genetic ID" is in the columns
        if "Genetic ID" not in df_annotation.columns:
            sys.exit("Error: Column 'Genetic ID' not found in the annotation file. Please check the headers.")

        # 3. Filter rows: keep only those whose "Genetic ID" is in sample_name_set
        df_filtered = df_annotation[df_annotation["Genetic ID"].isin(sample_name_set)]
        
        if df_filtered.empty:
            print("Warning: No matching 'Genetic ID' values found. The output file will be empty.")
        else:
            print(f"Retained {len(df_filtered)} rows after filtering.")

        # 4. Further filter rows by "ASSESSMENT" == "PASS"
        if "ASSESSMENT" not in df_filtered.columns:
            sys.exit("Error: Column 'ASSESSMENT' not found in the annotation file. Cannot filter by 'PASS'.")
        
        df_filtered = df_filtered[df_filtered["ASSESSMENT"] == "PASS"]
        if df_filtered.empty:
            print("Warning: No rows with 'PASS' assessment found. The output file will be empty.")
        else:
            print(f"Retained {len(df_filtered)} rows after filtering by ASSESSMENT == 'PASS'.")

        # 5. Save the filtered results to a new file (tab-delimited by default)
        try:
            df_filtered.to_csv(args.output, index=False, sep='\t')
            print(f"Filtered results have been saved to: {args.output}")
        except PermissionError:
            sys.exit(f"Error: Permission denied when writing to '{args.output}'.")
        except IOError as e:
            sys.exit(f"Error writing output file: {str(e)}")
            
    except Exception as e:
        sys.exit(f"An unexpected error occurred: {str(e)}")

if __name__ == "__main__":
    main()