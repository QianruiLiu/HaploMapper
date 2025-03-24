#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import folium
from folium import IFrame

def process_dataframe(path, sep):
    """
    Read CSV file and convert numeric columns
    """
    df = pd.read_csv(path, sep=sep, dtype=str).fillna("")
    print(f"Loaded {len(df)} rows from {path}.")

    meta_cols = ["Ancient pop name", "Country", "Age", "Lat", "Long"]
    all_cols = list(df.columns)
    haplo_cols = [c for c in all_cols if c not in meta_cols and c != "Total"]

    for c in haplo_cols:
        # Remove percentage signs and convert to numeric values
        df[c] = df[c].str.replace("%", "", regex=False)
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)

    if "Total" in df.columns:
        df["Total"] = pd.to_numeric(df["Total"], errors="coerce").fillna(0).astype(int)

    df["Lat"] = pd.to_numeric(df["Lat"], errors="coerce")
    df["Long"] = pd.to_numeric(df["Long"], errors="coerce")
    
    return df, haplo_cols

def create_popup_html(marker_id, pop_name, country, age, total, labels, values, haplo_type):
    """
    Construct an HTML fragment containing a <canvas> element and Chart.js script for displaying pie charts.
    Parameter haplo_type indicates the data type (Y-chromosome or mtDNA)
    """
    labels_js = ", ".join([f"'{lbl}'" for lbl in labels])
    values_js = ", ".join([f"{val:.2f}" for val in values])
    
    # Using Chart.js CDN, change to local reference if offline use is needed
    chart_js_cdn = "https://cdn.jsdelivr.net/npm/chart.js"

    html = f"""
        <h4>{pop_name} ({country})</h4>
        <p>Age: {age}, Total: {total}</p>
        <p>{haplo_type}</p>
        <canvas id="{marker_id}" width="250" height="200"></canvas>
        <script src="{chart_js_cdn}"></script>
        <script>
          var ctx = document.getElementById('{marker_id}').getContext('2d');
          var myPieChart = new Chart(ctx, {{
            type: 'pie',
            data: {{
              labels: [{labels_js}],
              datasets: [{{
                data: [{values_js}],
                backgroundColor: [
                  'rgba(220, 20, 60, 0.6)',
                  'rgba(30, 144, 255, 0.6)',
                  'rgba(34, 139, 34, 0.6)',
                  'rgba(255, 165, 0, 0.6)',
                  'rgba(128, 0, 128, 0.6)',
                  'rgba(255, 192, 203, 0.6)',
                  'rgba(0, 250, 154, 0.6)',
                  'rgba(70, 130, 180, 0.6)',
                  'rgba(106, 90, 205, 0.6)',
                  'rgba(188, 143, 143, 0.6)',
                  'rgba(95, 158, 160, 0.6)',
                  'rgba(128, 128, 128, 0.6)'
                ]
              }}]
            }},
            options: {{
              responsive: false,
              plugins: {{
                legend: {{ position: 'right' }}
              }}
            }}
          }});
        </script>
    """
    return html

def add_markers_to_group(df, haplo_cols, feature_group, id_prefix, haplo_type):
    """
    Iterate through data rows and add markers to the corresponding feature group for each row.
    id_prefix is used to generate unique canvas element ids
    haplo_type is used to display data type information in the popup
    """
    for idx, row in df.iterrows():
        pop_name = row["Ancient pop name"]
        country  = row["Country"]
        age      = row["Age"]
        lat      = row["Lat"]
        lon      = row["Long"]
        total    = row.get("Total", 1)
        
        freq_data = []
        for c in haplo_cols:
            val = row[c]
            if val > 0:
                freq_data.append((c, val))
        
        if not freq_data:
            continue
        
        labels = [x[0] for x in freq_data]
        values = [float(x[1]) for x in freq_data]
        
        marker_id = f"{id_prefix}_marker_{idx}"
        popup_html = create_popup_html(marker_id, pop_name, country, age, total, labels, values, haplo_type)
        
        # Skip if latitude or longitude is invalid
        if pd.isna(lat) or pd.isna(lon):
            continue
        
        # Use IFrame to embed HTML in the Popup
        iframe = IFrame(html=popup_html, width=300, height=300)
        popup = folium.Popup(iframe, max_width=450)
        
        folium.Marker(
            location=[lat, lon],
            tooltip=f"{pop_name} ({country})",
            popup=popup
        ).add_to(feature_group)

def main():
    parser = argparse.ArgumentParser(
        description="Plot pie charts of Y-chromosome and mtDNA haplogroups on a Folium map.")
    parser.add_argument("--y_input", required=True,
                        help="Y-chromosome haplogroup frequency CSV file path")
    parser.add_argument("--mt_input", required=True,
                        help="mtDNA haplogroup frequency CSV file path")
    parser.add_argument("--sep", default="\t",
                        help="CSV file separator (default: tab)")
    parser.add_argument("--output_html", default="haplogroup_map.html",
                        help="Output HTML filename")
    args = parser.parse_args()

    # Read two data files
    df_y, haplo_cols_y = process_dataframe(args.y_input, args.sep)
    df_mt, haplo_cols_mt = process_dataframe(args.mt_input, args.sep)
    
    # Calculate the center point (average latitude and longitude) of all data
    all_lat = pd.concat([df_y["Lat"], df_mt["Lat"]])
    all_lon = pd.concat([df_y["Long"], df_mt["Long"]])
    center_lat = all_lat.mean() if not all_lat.empty else 20
    center_lon = all_lon.mean() if not all_lon.empty else 0
    
    # Initialize Folium map
    folium_map = folium.Map(location=[center_lat, center_lon], zoom_start=3)
    
    # Create two layers: Y-chromosome and mtDNA
    fg_y = folium.FeatureGroup(name="Y-chromosome")
    fg_mt = folium.FeatureGroup(name="mtDNA")
    
    add_markers_to_group(df_y, haplo_cols_y, fg_y, "y", "Y-chromosome haplogroup")
    add_markers_to_group(df_mt, haplo_cols_mt, fg_mt, "mt", "mtDNA haplogroup")
    
    # Add layers to the map
    folium_map.add_child(fg_y)
    folium_map.add_child(fg_mt)
    
    # Add layer control panel for interactive layer switching
    folium.LayerControl().add_to(folium_map)
    
    folium_map.save(args.output_html)
    print(f"Map has been saved to {args.output_html}")

if __name__ == "__main__":
    main()