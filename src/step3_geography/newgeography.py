#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import random
import colorsys
import pandas as pd
import folium
from folium import IFrame

def generate_color_palette(labels):
    """
    Assign a color to each label in the given list.
    Use uniform sampling in the HSV color space, then convert to RGBA.
    """
    n = len(labels)
    if n == 0:
        return {}
    sorted_labels = sorted(labels)
    color_map = {}
    for i, lbl in enumerate(sorted_labels):
        h = i / float(n)
        s = 0.65 + 0.35 * random.random()
        v = 0.8 + 0.2 * random.random()
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        color_str = f"rgba({int(r*255)}, {int(g*255)}, {int(b*255)}, 0.75)"
        color_map[lbl] = color_str
    return color_map

def process_dataframe(path, sep):
    """
    Read the file and remove the '%' from haplo columns, converting them to float.
    """
    df = pd.read_csv(path, sep=sep, dtype=str).fillna("")
    print(f"Loaded {len(df)} rows from {path}.")

    meta_cols = ["Ancient pop name", "Country", "Age", "Lat", "Long"]
    all_cols = list(df.columns)
    haplo_cols = [c for c in all_cols if c not in meta_cols and c != "Total"]

    # Remove "%" and convert to numeric values
    for c in haplo_cols:
        df[c] = df[c].str.replace("%", "", regex=False)
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)

    if "Total" in df.columns:
        df["Total"] = pd.to_numeric(df["Total"], errors="coerce").fillna(0).astype(int)

    # Convert latitude and longitude to numeric values
    df["Lat"] = pd.to_numeric(df["Lat"], errors="coerce")
    df["Long"] = pd.to_numeric(df["Long"], errors="coerce")
    
    return df, haplo_cols

def is_basal_col(col):
    """Determine if the column represents a basal haplogroup (e.g. 'A', 'B', 'C', etc.)"""
    return len(col) == 1 and col.isalpha()

def is_subclade_col(col):
    """Determine if the column represents a subclade (e.g. 'A1', 'B2', etc.: letter followed by a digit)"""
    if len(col) < 2:
        return False
    return col[0].isalpha() and col[1].isdigit()

def build_two_ring_data(row, haplo_cols):
    """
    Given a row (representing a population) with frequency values,
    return:
      label_list  : All basal and subclade labels present in the population
                    (first all basal, then all subclade, only those with nonzero values)
      data_main   : Percentage of basal haplogroups in the population (inner ring)
      data_sub    : Percentage of subclades in the population (outer ring), taken directly from the file
    """
    basal_data = []     # (col, pct, abs_count)
    subclade_data = []  # (col, pct, abs_count)

    total_ind = row.get("Total", 0)
    if not isinstance(total_ind, (int, float)):
        total_ind = 0

    # Store the absolute count for each basal haplogroup (used for calculating subclade percentages later)
    store_basal_count = {}

    # Interpret row[col] as the overall percentage (0-100) and compute the absolute count:
    for col in haplo_cols:
        val_pct = float(row[col])
        if val_pct <= 0:
            continue
        absolute_count = val_pct / 100.0 * total_ind if total_ind > 0 else val_pct

        if is_basal_col(col):
            basal_data.append((col, val_pct, absolute_count))
            store_basal_count[col] = absolute_count
        elif is_subclade_col(col):
            subclade_data.append((col, val_pct, absolute_count))
        else:
            pass

    # Optional sorting
    basal_data.sort(key=lambda x: x[1], reverse=True)
    subclade_data.sort(key=lambda x: x[1], reverse=True)

    # Basal data (inner ring)
    label_list = [b[0] for b in basal_data]
    data_main  = [b[1] for b in basal_data]  # Basal percentages from the file

    # Pre-fill outer ring data with zeros
    data_sub = [0] * len(basal_data)

    # Simple function to get the parent letter for a subclade
    def get_parent_letter(sub):
        return sub[0]  # For example, A1 -> A

    # Collect subclade data (directly use the percentage values stored in the input file)
    for (col_s, pct_s, abs_s) in subclade_data:
        # Directly use the subclade percentage as outer ring data
        label_list.append(col_s)
        data_main.append(0)      # Subclades have 0 in the inner ring
        data_sub.append(pct_s)   # Outer ring uses the percentage from the file

    return label_list, data_main, data_sub

def create_popup_html(marker_id, pop_name, country, age, total,
                      label_list, data_main, data_sub, color_list,
                      haplo_type):
    """
    Construct HTML to display a double-layer doughnut chart:
      dataset[0] => subclades (outer ring)
      dataset[1] => basal haplogroups (inner ring)

    Use the DataLabels plugin to differentiate the display:
      - For basal haplogroups (datasetIndex=1), display only the label.
      - For subclades (datasetIndex=0), display "label + percentage".
    """
    chart_js_cdn = "https://cdn.jsdelivr.net/npm/chart.js"
    datalabels_cdn = "https://cdn.jsdelivr.net/npm/chartjs-plugin-datalabels@2.0.0"

    labels_js = ", ".join([f"'{lbl}'" for lbl in label_list])
    main_js = ", ".join([f"{v:.2f}" for v in data_main])
    sub_js = ", ".join([f"{v:.2f}" for v in data_sub])
    colors_js = ", ".join([f"'{c}'" for c in color_list])

    # Calculate the sum of the outer ring values
    sum_sub = sum(data_sub)

    # dataset 0 => subclades (outer ring)
    ds_sub = f"""
      {{
        label: 'Subclades (Outside)',
        data: [{sub_js}],
        backgroundColor: [{colors_js}],
        borderWidth: 1
      }}
    """
    # dataset 1 => basal haplogroups (inner ring)
    ds_main = f"""
      {{
        label: 'Basal (Inside)',
        data: [{main_js}],
        backgroundColor: [{colors_js}],
        borderWidth: 1
      }}
    """

    # If the outer ring sums to 0, only draw the basal ring
    if sum_sub > 0:
        datasets_code = f"[{ds_sub}, {ds_main}]"
    else:
        datasets_code = f"[{ds_main}]"

    html = f"""
    <h4>{pop_name} ({country})</h4>
    <p>Age: {age}, Total: {total}</p>
    <p>{haplo_type}</p>
    <canvas id="{marker_id}" width="600" height="550"></canvas>

    <!-- Include Chart.js and DataLabels plugin -->
    <script src="{chart_js_cdn}"></script>
    <script src="{datalabels_cdn}"></script>
    <script>
      Chart.register(ChartDataLabels);

      var ctx = document.getElementById('{marker_id}').getContext('2d');
      var myChart = new Chart(ctx, {{
        type: 'doughnut',
        data: {{
          labels: [{labels_js}],
          datasets: {datasets_code}
        }},
        options: {{
          responsive: false,
          plugins: {{
            legend: {{ position: 'right' }},
            datalabels: {{
              display: function(context) {{
                var val = context.dataset.data[context.dataIndex];
                return val > 0; // Only display labels for values > 0
              }},
              formatter: function(value, context) {{
                // datasetIndex=0 => subclades (outer ring), datasetIndex=1 => basal haplogroups (inner ring)
                var dsIndex = context.datasetIndex;
                var lbl = context.chart.data.labels[context.dataIndex];
                if(dsIndex === 0) {{
                  // For subclades, display "label + percentage"
                  return lbl + ' ' + value.toFixed(1) + '%';
                }} else {{
                  // For basal haplogroups, display only the label
                  return lbl;
                }}
              }},
              color: '#000',
              font: {{
                weight: 'bold',
                size: 12
              }}
            }}
          }}
        }}
      }});
    </script>
    """
    return html

def add_markers_to_group(df, haplo_cols, feature_group, id_prefix, haplo_type, color_map):
    """
    Build and add markers.
      - Inner ring: basal haplogroups (datasetIndex=1, display label only)
      - Outer ring: subclades (datasetIndex=0, display label + percentage)
    """
    for idx, row in df.iterrows():
        pop_name = row["Ancient pop name"]
        country  = row["Country"]
        age = row["Age"]
        lat = row["Lat"]
        lon = row["Long"]
        total = row.get("Total", 0)

        if pd.isna(lat) or pd.isna(lon):
            continue

        # Calculate label_list, data_main, and data_sub
        label_list, data_main, data_sub = build_two_ring_data(row, haplo_cols)
        if not label_list:
            continue

        # Retrieve colors in the same order as labels
        color_list = [color_map[lbl] for lbl in label_list]

        marker_id = f"{id_prefix}_marker_{idx}"
        popup_html = create_popup_html(
            marker_id=marker_id,
            pop_name=pop_name,
            country=country,
            age=age,
            total=total,
            label_list=label_list,
            data_main=data_main,      # Basal (inner ring)
            data_sub=data_sub,        # Subclades (outer ring)
            color_list=color_list,
            haplo_type=haplo_type
        )

        iframe = IFrame(html=popup_html, width=750, height=800)
        popup = folium.Popup(iframe, max_width=850)

        folium.Marker(
            location=[lat, lon],
            tooltip=f"{pop_name} ({country})",
            popup=popup
        ).add_to(feature_group)

def main():
    parser = argparse.ArgumentParser(
        description="Plot a multi-layer doughnut chart on a Folium map for basal and subclade haplogroups: inner ring shows basal labels only, outer ring shows subclade percentages."
    )
    parser.add_argument("--y_input", required=True,
                        help="Y-chromosome haplogroup frequency file (CSV/TSV)")
    parser.add_argument("--mt_input", required=True,
                        help="mtDNA haplogroup frequency file (CSV/TSV)")
    parser.add_argument("--sep", default="\t",
                        help="Delimiter for the input file (default is \\t). For CSV, use --sep ','")
    parser.add_argument("--output_html", default="haplogroup_map.html",
                        help="Output HTML filename")
    args = parser.parse_args()

    # Read the two files
    df_y, haplo_cols_y = process_dataframe(args.y_input, args.sep)
    df_mt, haplo_cols_mt = process_dataframe(args.mt_input, args.sep)

    # Collect all labels (basal + subclade) for unified color mapping
    all_labels = sorted(set(haplo_cols_y + haplo_cols_mt))
    color_map = generate_color_palette(all_labels)
    
    # Calculate the center of the map
    all_lat = pd.concat([df_y["Lat"], df_mt["Lat"]])
    all_lon = pd.concat([df_y["Long"], df_mt["Long"]])
    center_lat = all_lat.mean() if not all_lat.empty else 20
    center_lon = all_lon.mean() if not all_lon.empty else 0
    
    # Initialize the map
    folium_map = folium.Map(location=[center_lat, center_lon], zoom_start=3)
    
    # Create feature groups
    fg_y  = folium.FeatureGroup(name="Y-chr")
    fg_mt = folium.FeatureGroup(name="mtDNA")

    # Add markers to each feature group
    add_markers_to_group(df_y,  haplo_cols_y,  fg_y,  "y",  "Y-chr haplogroup",  color_map)
    add_markers_to_group(df_mt, haplo_cols_mt, fg_mt, "mt", "mtDNA haplogroup", color_map)
    
    folium_map.add_child(fg_y)
    folium_map.add_child(fg_mt)
    folium.LayerControl().add_to(folium_map)
    
    folium_map.save(args.output_html)
    print(f"The map has been saved to {args.output_html}")

if __name__ == "__main__":
    main()
