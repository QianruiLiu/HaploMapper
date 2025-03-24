#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functionality:
  1) Reads a filtered annotation file (e.g., annotation_filtered.tsv).
  2) Produces two frequency tables:
       - Y-chromosome haplogroups (basal letter)
       - mtDNA haplogroups (basal letter)
     Each table has columns exactly as shown in the example image:
       Ancient pop name | Country | Age | Lat | Long | A | B | C | ... | Total

  3) (ADDED) For each basal haplogroup (e.g., "A"), also parse and compute 
     the frequency of any 'subclade' in the form "A1", "A2", etc. 
     (i.e., if the second character is a digit). These subclade frequencies 
     are calculated relative to the count of the basal haplogroup only (not the overall total).

Notes:
  - Y-chromosome and mtDNA frequency tables are output into two separate files.
  - This script assumes that the following columns exist in the input:
      "Group ID"
      "Political Entity"  -> used for "Country" and for forming "Ancient pop name"
      "Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]" -> for age
      "Lat." -> "Lat"
      "Long." -> "Long"
      "Y haplogroup (manual curation in ISOGG format)" -> for Y-chromosome
      "mtDNA haplogroup if >2x or published" -> for mtDNA

Usage example:
  python new_haplogroup_frequency_tables.py --input annotation_filtered.tsv [--sep "\t"]
"""

import argparse
import pandas as pd
import numpy as np
import os
import sys

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Generate basal and subclade haplogroup frequency tables (Y & mtDNA).")
    parser.add_argument(
        "--input",
        required=True,
        help="Path to the filtered file (CSV/TSV) that includes columns for haplogroups and metadata.")
    parser.add_argument(
        "--sep",
        default="\t",
        help="Delimiter for the input file. Default is tab ('\\t').")
    parser.add_argument(
        "--y_output",
        default="Y_haplogroup_frequencies.tsv",
        help="Output file name for the Y-chromosome haplogroup table.")
    parser.add_argument(
        "--mt_output",
        default="mtDNA_haplogroup_frequencies.tsv",
        help="Output file name for the mtDNA haplogroup table.")
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input):
        sys.exit(f"Error: Input file '{args.input}' does not exist.")
    
    try:
        # 1. Read the filtered file
        df = pd.read_csv(args.input, sep=args.sep, dtype=str)
        print(f"Loaded {len(df)} rows from {args.input}.")
    except pd.errors.EmptyDataError:
        sys.exit(f"Error: Input file '{args.input}' is empty.")
    except pd.errors.ParserError:
        sys.exit(f"Error: Could not parse '{args.input}'. Please check the file format and delimiter.")
    except Exception as e:
        sys.exit(f"Error reading input file: {str(e)}")
    
    # Fill NaN values with empty strings
    df = df.fillna("")
    
    # 2. Define relevant column names from the input file
    col_group_id = "Group ID"
    col_country = "Political Entity"
    col_age = "Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]"
    col_lat = "Lat."
    col_long = "Long."
    col_y_hg = "Y haplogroup (manual curation in ISOGG format)"
    col_mt_hg = "mtDNA haplogroup if >2x or published"
    
    # Check for required columns
    required_columns = [col_country, col_age, col_lat, col_long]
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        sys.exit(f"Error: The following required columns are missing: {', '.join(missing_columns)}")
    
    try:
        # Preprocess age and generate aggregated group labels 
        # Convert age to numeric (in BP)
        df['age_numeric'] = pd.to_numeric(df[col_age], errors='coerce')
        
        # Check if all age values are NaN
        if df['age_numeric'].isna().all():
            sys.exit(f"Error: No valid age values found in column '{col_age}'.")
            
        # Compute the lower bound for each 1000 BP bin
        df['age_bin_lower'] = (np.floor(df['age_numeric'] / 1000) * 1000).astype('Int64')
        # Create a binned age string, e.g., "0-1000 BP"
        df["binned_age"] = df['age_bin_lower'].astype(str) + "-" + (df['age_bin_lower'] + 1000).astype(str) + " BP"
        # Generate new "Ancient pop name": combine Political Entity and binned age
        df["Ancient pop name"] = df[col_country] + " " + df["binned_age"]
        
        # Convert Latitude and Longitude to numeric for later aggregation
        df["lat_numeric"] = pd.to_numeric(df[col_lat], errors='coerce')
        df["long_numeric"] = pd.to_numeric(df[col_long], errors='coerce')
    except Exception as e:
        sys.exit(f"Error during data preprocessing: {str(e)}")
    
    # Verify that necessary columns exist in the DataFrame
    needed_cols = ["Ancient pop name", col_country, "binned_age"]
    for c in needed_cols:
        if c not in df.columns:
            raise ValueError(f"Missing column '{c}' in your input file. Check your headers.")
    
    # Define grouping keys: Ancient pop name, Political Entity, and binned_age
    group_keys = ["Ancient pop name", col_country, "binned_age"]
    
    # 3. Build the Y-chromosome haplogroup frequency table
    if col_y_hg not in df.columns:
        print(f"Warning: Y haplogroup column '{col_y_hg}' not found in input file. Skipping Y frequency table.")
    elif all(df[col_y_hg].str.strip() == ""):
        print(f"Warning: No Y haplogroup info found under '{col_y_hg}'. Skipping Y frequency table.")
    else:
        try:
            make_haplo_table(
                df=df,
                haplo_col=col_y_hg,
                group_cols=group_keys,
                output_file=args.y_output,
                label="Y-chr"
            )
        except Exception as e:
            print(f"Error creating Y haplogroup table: {str(e)}")
    
    # 4. Build the mtDNA haplogroup frequency table
    if col_mt_hg not in df.columns:
        print(f"Warning: mtDNA haplogroup column '{col_mt_hg}' not found in input file. Skipping mtDNA frequency table.")
    elif all(df[col_mt_hg].str.strip() == ""):
        print(f"Warning: No mtDNA haplogroup info found under '{col_mt_hg}'. Skipping mtDNA frequency table.")
    else:
        try:
            make_haplo_table(
                df=df,
                haplo_col=col_mt_hg,
                group_cols=group_keys,
                output_file=args.mt_output,
                label="mtDNA"
            )
        except Exception as e:
            print(f"Error creating mtDNA haplogroup table: {str(e)}")

def make_haplo_table(df, haplo_col, group_cols, output_file, label):
    """
    Groups the DataFrame by group_cols and computes the frequency (in percent)
    of the basal haplogroup letter (first character in haplo_col).
    It also parses any subclade (if the second character is a digit, e.g., "A1", "B2"),
    and computes its frequency relative to the parent's raw count.
    
    Then, it calculates the calendar year using: 1950 - mean(age_numeric),
    formatting the result with "CE" or "BCE" suffix.
    The final table includes columns:
      Ancient pop name | Country | Age | Lat | Long | A | A1 | A2 | B | B1 | B2 | ... | Total
    where A, B, etc. are basal haplogroups, and A1, A2, etc. are subclade columns appended after each basal group.
    
    Note:
    - The "Total" is recalculated as the sum of the raw counts of the basal haplogroups and their subclades.
    - All aggregations (including coordinate averages) are based only on the filtered data.
    """
    print(f"Building {label} haplogroup frequency table...")

    # 1) Filter out invalid entries in haplo_col (e.g., empty strings, "..", "n/a", etc.)
    invalid_patterns = ["", "..", "n/a", "N/A", "Neanderthal"]
    df_valid = df[~df[haplo_col].str.strip().isin(invalid_patterns)].copy()
    
    # Check if there are any valid entries left after filtering
    if len(df_valid) == 0:
        raise ValueError(f"No valid haplogroup entries found in column '{haplo_col}' after filtering.")
    
    # Apply more complex filtering with error handling
    def is_valid_haplo(x):
        try:
            x = str(x).strip()
            return len(x) == 1 or (len(x) >= 2 and x[1].isdigit())
        except:
            return False
            
    df_valid = df_valid[df_valid[haplo_col].apply(is_valid_haplo)]
    
    # Check again after the second filter
    if len(df_valid) == 0:
        raise ValueError(f"No valid haplogroup entries matching the basal or subclade pattern in column '{haplo_col}'.")

    # 2) Extract the basal haplogroup letter and convert to uppercase
    df_valid["haplo_basal"] = df_valid[haplo_col].str[0].str.upper()

    def parse_subclade(h):
        try:
            h = str(h).strip()
            if len(h) < 2:
                return None
            if h[1].isdigit():
                return h[:2].upper()
            return None
        except:
            return None

    df_valid["haplo_sub"] = df_valid[haplo_col].apply(parse_subclade)

    try:
        # 3) Group by the specified keys and basal haplogroup; count occurrences (raw count)
        grouped_basal = df_valid.groupby(group_cols + ["haplo_basal"]).size().reset_index(name="Count_basal")
        
        # Check if grouping produced any results
        if len(grouped_basal) == 0:
            raise ValueError(f"Grouping by {group_cols} and haplogroups produced no results.")
            
        # 4) Pivot the table: rows = group_cols, columns = haplo_basal, values = Count_basal
        pivot_basal = grouped_basal.pivot_table(
            index=group_cols,
            columns="haplo_basal",
            values="Count_basal",
            fill_value=0
        )
    except Exception as e:
        raise ValueError(f"Error during grouping or pivoting: {str(e)}")
        
    # Save a copy of raw counts for later calculations (including the basal columns and total count)
    pivot_basal_counts = pivot_basal.copy()
    pivot_basal_counts["Total_count"] = pivot_basal_counts.sum(axis=1)

    # 5) Compute the percentage for basal haplogroups (for output); "Total" is the sum of basal counts
    pivot_basal["Total"] = pivot_basal.sum(axis=1)
    
    # Check for division by zero
    if (pivot_basal["Total"] == 0).any():
        print(f"Warning: Some rows have a total count of 0. Percentages for these rows will be set to 0.")
        
    for col in pivot_basal.columns:
        if col != "Total":
            # Avoid division by zero
            pivot_basal[col] = np.where(
                pivot_basal["Total"] > 0,
                pivot_basal[col] / pivot_basal["Total"] * 100,
                0.0
            )

    pivot_basal.reset_index(inplace=True)
    pivot_basal_counts.reset_index(inplace=True)

    # 6) Count the occurrences for subclades (raw count)
    try:
        # Filter out rows where haplo_sub is NaN
        sub_valid = df_valid.dropna(subset=["haplo_sub"])
        
        # Check if there are any subclades
        if len(sub_valid) == 0:
            print(f"Warning: No valid subclade haplogroups found for {label}.")
            grouped_sub = pd.DataFrame(columns=group_cols + ["haplo_basal", "haplo_sub", "Count_sub"])
        else:
            grouped_sub = (
                sub_valid
                .groupby(group_cols + ["haplo_basal", "haplo_sub"])
                .size()
                .reset_index(name="Count_sub")
            )
            
        if len(grouped_sub) > 0:
            pivot_sub = grouped_sub.pivot_table(
                index=group_cols,
                columns="haplo_sub",
                values="Count_sub",
                fill_value=0
            ).reset_index()
        else:
            # Create an empty pivot table with the correct structure
            pivot_sub = pd.DataFrame(columns=group_cols)
    except Exception as e:
        print(f"Warning: Error processing subclades: {str(e)}. Continuing without subclade data.")
        pivot_sub = pd.DataFrame(columns=group_cols)

    # 7) Merge the basal percentage table and subclade count table (subclade counts are still raw)
    try:
        merged_df = pd.merge(pivot_basal, pivot_sub, on=group_cols, how="left", suffixes=("", "_sub"))

        # 8) Merge the raw basal counts (with a "_count" suffix) into the table
        merged_df = pd.merge(merged_df, pivot_basal_counts, on=group_cols, how="left", suffixes=("", "_count"))
    except Exception as e:
        raise ValueError(f"Error merging tables: {str(e)}")

    # 9) Recalculate "Total": total = raw basal total count ("Total_count") + sum of all subclade raw counts
    subclade_cols = sorted(set(pivot_sub.columns) - set(group_cols)) if hasattr(pivot_sub, 'columns') else []
    merged_df["Total"] = merged_df["Total_count"]
    
    # Only add subclade counts if they exist
    if subclade_cols:
        merged_df["Total"] += merged_df[subclade_cols].sum(axis=1)

    # 10) For each subclade, compute its frequency: (subclade count / parent's raw count) * 100
    for sc in subclade_cols:
        try:
            parent_letter = sc[0].upper()
            merged_df[sc] = pd.to_numeric(merged_df[sc], errors="coerce").fillna(0)
            parent_count = pd.to_numeric(merged_df[parent_letter + "_count"], errors="coerce").fillna(0)
            
            # Avoid division by zero
            merged_df[sc] = np.where(
                parent_count > 0,
                merged_df[sc] / parent_count * 100,
                0.0
            )
        except Exception as e:
            print(f"Warning: Error processing subclade {sc}: {str(e)}. Setting values to 0.")
            merged_df[sc] = 0.0

    # 11) Merge coordinate data using the filtered dataset (df_valid)
    try:
        coord_df = df_valid.groupby(group_cols)[["lat_numeric", "long_numeric", "age_numeric"]].mean().reset_index()
        merged_df = pd.merge(merged_df, coord_df, on=group_cols, how="left")
    except Exception as e:
        print(f"Warning: Error processing coordinate data: {str(e)}. Coordinates may be missing or incorrect.")

    # 12) Rename the second grouping key column to "Country"
    merged_df.rename(columns={group_cols[1]: "Country"}, inplace=True)

    # 13) Calculate calendar years using age_numeric (1950 - age) and format with "CE"/"BCE" suffix
    try:
        computed_age = (1950 - merged_df["age_numeric"]).round().astype(int)
        def format_calendar_year(year):
            try:
                if year > 0:
                    return f"{year} CE"
                elif year < 0:
                    return f"{abs(year)} BCE"
                else:
                    return "0 CE"
            except:
                return "Unknown"
                
        merged_df["Age"] = computed_age.apply(format_calendar_year)
    except Exception as e:
        print(f"Warning: Error calculating calendar years: {str(e)}. Setting 'Age' to 'Unknown'.")
        merged_df["Age"] = "Unknown"

    # 14) Drop columns that are no longer needed: "binned_age" and "age_numeric"
    try:
        merged_df.drop(columns=[group_cols[2], "age_numeric"], inplace=True)
    except KeyError as e:
        print(f"Warning: Could not drop column: {str(e)}. Continuing.")

    # 15) Rename coordinate columns to "Lat" and "Long"
    merged_df.rename(columns={"lat_numeric": "Lat", "long_numeric": "Long"}, inplace=True)

    # 16) Construct the final output column order: metadata + basal haplogroups and their subclades + Total
    base_cols = ["Ancient pop name", "Country", "Age", "Lat", "Long"]
    basal_letters = [c for c in pivot_basal.columns if c not in group_cols and c != "Total"]
    final_order_of_haplos = []
    for b in sorted(basal_letters):
        final_order_of_haplos.append(b)
        related_subs = sorted([sc for sc in subclade_cols if sc.startswith(b)])
        final_order_of_haplos.extend(related_subs)
    final_order_of_haplos.append("Total")
    
    # Make sure all requested columns exist in the dataframe
    final_cols = [c for c in base_cols + final_order_of_haplos if c in merged_df.columns]
    
    # Check if we have the minimum required columns
    required_result_cols = ["Ancient pop name", "Country", "Total"]
    missing_required = [c for c in required_result_cols if c not in final_cols]
    if missing_required:
        raise ValueError(f"Missing required columns in final output: {missing_required}")
        
    merged_df = merged_df[final_cols]

    # 17) Format each haplogroup column: convert to percentage with two decimal places and append "%" (except for Total which is a raw count)
    for hc in final_order_of_haplos:
        if hc == "Total":
            try:
                merged_df[hc] = merged_df[hc].astype(int)
            except Exception as e:
                print(f"Warning: Error formatting Total column: {str(e)}. Values may be incorrect.")
        else:
            if hc in merged_df.columns:
                try:
                    merged_df[hc] = merged_df[hc].astype(float).apply(lambda x: f"{x:.2f}%" if pd.notna(x) else "")
                except Exception as e:
                    print(f"Warning: Error formatting column {hc}: {str(e)}. Values may be incorrect.")

    # 18) Output the final merged DataFrame to a TSV file
    try:
        # Check if output directory exists, create if it doesn't
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        merged_df.to_csv(output_file, sep="\t", index=False)
        print(f"{label} haplogroup table -> {output_file} (rows: {len(merged_df)})")
    except Exception as e:
        sys.exit(f"Error writing output file '{output_file}': {str(e)}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.exit(f"Fatal error: {str(e)}")