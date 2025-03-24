#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Functionality:
  1) Read a filtered annotation file (e.g., annotation_filtered.tsv).
  2) Produce two frequency tables:
       - Y-chr haplogroups (basal letter)
       - mtDNA haplogroups (basal letter)
     Each table has columns exactly as shown in the example image:
       Ancient pop name | Country | Age | Lat | Long | A | B | C | ... | Total

Notes:
  - We separate Y-chr and mtDNA into two output files.
  - This script assumes the following columns exist in the input:
      "Group ID" -> originally provided but not used for naming groups now
      "Political Entity" -> used for "Country" and to form the Ancient pop name
      "Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]" -> used for age information
      "Lat." -> "Lat"
      "Long." -> "Long"
      "Y haplogroup (manual curation in ISOGG format)" -> for Y-chr
      "mtDNA haplogroup if >2x or published" -> for mtDNA

Usage example:
  python build_haplogroup_tables.py --input annotation_filtered.tsv --sep "\t"
"""

import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Generate basal haplogroup frequency tables (Y & mtDNA).")
    parser.add_argument(
        "--input",
        required=True,
        help="Path to the filtered file (CSV/TSV) that includes columns for haplogroups and metadata."
    )
    parser.add_argument(
        "--sep",
        default="\t",
        help="Delimiter for the input file. Default is tab ('\\t')."
    )
    parser.add_argument(
        "--y_output",
        default="Y_haplogroup_frequencies.tsv",
        help="Output file name for the Y-chr haplogroup table."
    )
    parser.add_argument(
        "--mt_output",
        default="mtDNA_haplogroup_frequencies.tsv",
        help="Output file name for the mtDNA haplogroup table."
    )
    args = parser.parse_args()
    
    # === 1. Read the filtered file ===
    df = pd.read_csv(args.input, sep=args.sep, dtype=str).fillna("")
    print(f"Loaded {len(df)} rows from {args.input}.")
    
    # === 2. Define relevant column names ===
    col_group_id = "Group ID"
    col_country  = "Political Entity"
    col_age      = "Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]"
    col_lat      = "Lat."
    col_long     = "Long."
    col_y_hg     = "Y haplogroup (manual curation in ISOGG format)"
    col_mt_hg    = "mtDNA haplogroup if >2x or published"
    
    # --- Preprocess age and generate aggregated group labels ---
    # Convert age to numeric (BP)
    df['age_numeric'] = pd.to_numeric(df[col_age], errors='coerce')
    # Compute lower bound for each 1000 BP bin
    df['age_bin_lower'] = (np.floor(df['age_numeric'] / 1000) * 1000).astype('Int64')
    # Create a binned age string, e.g., "0-1000 BP"
    df["binned_age"] = df['age_bin_lower'].astype(str) + "-" + (df['age_bin_lower'] + 1000).astype(str) + " BP"
    # Generate new Ancient pop name: Political Entity + " " + binned_age
    df["Ancient pop name"] = df[col_country] + " " + df["binned_age"]
    
    # Convert Lat and Long to numeric (for later aggregation)
    df["lat_numeric"] = pd.to_numeric(df[col_lat], errors='coerce')
    df["long_numeric"] = pd.to_numeric(df[col_long], errors='coerce')
    
    # Check required columns (using aggregated columns)
    needed_cols = [ "Ancient pop name", col_country, "binned_age" ]
    for c in needed_cols:
        if c not in df.columns:
            raise ValueError(f"Missing column '{c}' in your input file. Check your headers.")
    
    # --- Define grouping keys ---
    # Group by: Ancient pop name, Political Entity, and binned_age.
    group_keys = ["Ancient pop name", col_country, "binned_age"]
    
    # === 3. Build Y-chromosome haplogroup frequency table ===
    if col_y_hg not in df.columns or all(df[col_y_hg].str.strip() == ""):
        print(f"No Y haplogroup info found under '{col_y_hg}'. Skipping Y frequency table.")
    else:
        make_haplo_table(
            df=df,
            haplo_col=col_y_hg,
            group_cols=group_keys,
            output_file=args.y_output,
            label="Y-chr"
        )
    
    # === 4. Build mtDNA haplogroup frequency table ===
    if col_mt_hg not in df.columns or all(df[col_mt_hg].str.strip() == ""):
        print(f"No mtDNA haplogroup info found under '{col_mt_hg}'. Skipping mtDNA frequency table.")
    else:
        make_haplo_table(
            df=df,
            haplo_col=col_mt_hg,
            group_cols=group_keys,
            output_file=args.mt_output,
            label="mtDNA"
        )

def make_haplo_table(df, haplo_col, group_cols, output_file, label):
    """
    Group the given DataFrame by group_cols and compute the frequency (in percent)
    of the basal haplogroup letter (the first letter in haplo_col).
    Also aggregate latitude, longitude, and age_numeric (mean) for each group.
    Then calculate the calendar year using: 1950 - mean(age_numeric)
    and add a suffix ("CE" for positive years, "BCE" for negative years).
    The final table contains:
      Ancient pop name | Country | Age | Lat | Long | haplo columns... | Total
    """
    print(f"Building {label} haplogroup frequency table...")
    
    # 1) Filter out invalid entries in haplo_col (e.g., "", "..", "n/a")
    invalid_patterns = ["", "..", "n/a", "N/A"]
    df_valid = df[~df[haplo_col].str.strip().isin(invalid_patterns)].copy()
    
    # 2) Extract the basal haplogroup letter and convert to uppercase
    df_valid["haplo_basal"] = df_valid[haplo_col].str[0].str.upper()
    
    # 3) Group by the specified keys and haplo_basal, then count occurrences
    grouped = df_valid.groupby(group_cols + ["haplo_basal"]).size().reset_index(name="Count")
    
    # 4) Pivot: rows = group_cols, columns = haplo_basal, values = Count
    pivot_df = grouped.pivot_table(
        index=group_cols,
        columns="haplo_basal",
        values="Count",
        fill_value=0
    )
    
    # 5) Calculate total count per group
    pivot_df["Total"] = pivot_df.sum(axis=1)
    
    # 6) Convert haplo counts to percentages (except Total)
    for col in pivot_df.columns:
        if col != "Total":
            pivot_df[col] = pivot_df[col] / pivot_df["Total"] * 100
    
    # 7) Reset index to make group_cols regular columns
    pivot_df.reset_index(inplace=True)
    
    # 8) Aggregate lat, long, and age_numeric (mean) for each group
    coord_df = df.groupby(group_cols)[["lat_numeric", "long_numeric", "age_numeric"]].mean().reset_index()
    
    # 9) Merge the haplogroup pivot table with the aggregated coordinates
    merged_df = pd.merge(pivot_df, coord_df, on=group_cols, how="left")
    
    # 10) Rename columns: use group_cols[1] as Country
    merged_df.rename(columns={
        group_cols[1]: "Country"
    }, inplace=True)
    
    # 11) Compute calendar year using the mean age_numeric: 1950 - mean(age_numeric)
    computed_age = (1950 - merged_df["age_numeric"]).round().astype(int)
    # 根据计算结果添加后缀：大于0加"CE"，小于0加"BCE"，等于0则显示"0 CE"
    def format_calendar_year(year):
        if year > 0:
            return f"{year} CE"
        elif year < 0:
            return f"{abs(year)} BCE"
        else:
            return "0 CE"
    merged_df["Age"] = computed_age.apply(format_calendar_year)
    
    # 12) 删除不再需要的辅助列
    merged_df.drop(columns=[group_cols[2], "age_numeric"], inplace=True)
    
    # 13) Rename coordinate columns
    merged_df.rename(columns={
        "lat_numeric": "Lat",
        "long_numeric": "Long"
    }, inplace=True)
    
    # 14) Define final column order:
    # "Ancient pop name", "Country", "Age", "Lat", "Long", haplogroup columns..., "Total"
    base_cols = ["Ancient pop name", "Country", "Age", "Lat", "Long"]
    all_cols = list(merged_df.columns)
    haplo_cols = [col for col in all_cols if col not in base_cols + ["Total"]]
    haplo_cols_sorted = sorted(haplo_cols)
    final_cols = base_cols + haplo_cols_sorted + ["Total"]
    merged_df = merged_df[final_cols]
    
    # 15) Format haplo percentage columns to two decimals with "%" sign
    for col in haplo_cols_sorted:
        merged_df[col] = merged_df[col].astype(float).apply(lambda x: f"{x:.2f}%" if pd.notna(x) else "")
    
    # 16) Convert Total to integer
    merged_df["Total"] = merged_df["Total"].astype(int)
    
    # 17) Output TSV file
    merged_df.to_csv(output_file, sep="\t", index=False)
    print(f"{label} haplogroup table -> {output_file} (rows: {len(merged_df)})")

if __name__ == "__main__":
    main()



