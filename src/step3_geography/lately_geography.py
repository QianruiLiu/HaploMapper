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
    Generate unified label arrays and two-layer data based on input row:
      label_list : Ordered, with parent groups followed by their subclades
      data_main  : Inner ring data, parent groups use actual percentages, subclades filled with 0
      data_sub   : Outer ring data, parent groups filled with 0, subclades use normalized values (normalized sum equals parent group percentage)
      display_sub: Percentages for display in the outer ring, parent groups are None, subclades are (subclade original value/total of parent group's subclades)*100
    """
    basal_values = {}
    subclade_values = {}
    
    for col in haplo_cols:
        try:
            val_pct = float(row[col])
        except:
            continue
        if val_pct <= 0:
            continue
        
        if is_basal_col(col):
            basal_values[col] = val_pct
        elif is_subclade_col(col):
            parent = col[0]  # For example, the parent group of A1 is A
            subclade_values.setdefault(parent, []).append((col, val_pct))
    
    basal_sorted = sorted(basal_values.items(), key=lambda x: x[1], reverse=True)
    
    label_list = []
    data_main = []
    data_sub = []
    display_sub = []  # For displaying subclade percentages (as a percentage of parent group)
    
    for basal, basal_val in basal_sorted:
        # Add parent group: inner ring shows actual value, outer ring filled with 0; display data doesn't show percentage
        label_list.append(basal)
        data_main.append(basal_val)
        data_sub.append(0)
        display_sub.append(None)
        
        if basal in subclade_values:
            subclades = subclade_values[basal]
            sum_sub = sum([val for (_, val) in subclades])
            for sc, val in subclades:
                label_list.append(sc)
                data_main.append(0)  # Fill subclade position with 0 in inner ring
                if sum_sub > 0:
                    rel_val = (val / sum_sub) * basal_val
                    disp_val = (val / sum_sub) * 100  # Calculate subclade percentage relative to parent group
                else:
                    rel_val = 0
                    disp_val = 0
                data_sub.append(rel_val)
                display_sub.append(disp_val)
    return label_list, data_main, data_sub, display_sub


def create_popup_html(marker_id, pop_name, country, age, total,
                      label_list, data_main, data_sub, display_sub, color_list,
                      haplo_type):
    """
    Generate double-layer pie chart HTML, with inner and outer rings aligned, 
    and outer ring labels showing subclade percentages relative to parent group.
    """
    chart_js_cdn = "https://cdn.jsdelivr.net/npm/chart.js"
    datalabels_cdn = "https://cdn.jsdelivr.net/npm/chartjs-plugin-datalabels@2.0.0"

    labels_js = ", ".join([f"'{lbl}'" for lbl in label_list])
    main_js = ", ".join([f"{v:.2f}" for v in data_main])
    sub_js = ", ".join([f"{v:.2f}" for v in data_sub])
    colors_js = ", ".join([f"'{c}'" for c in color_list])
    display_js = ", ".join([f"{(v if v is not None else 'null')}" for v in display_sub])

    # dataset 0: outer ring (subclades); dataset 1: inner ring (parent groups)
    ds_sub = f"""
      {{
        label: 'Subclades (Outside)',
        data: [{sub_js}],
        backgroundColor: [{colors_js}],
        borderAlign: 'inner',
        borderWidth: 1
      }}
    """
    ds_main = f"""
      {{
        label: 'Basal (Inside)',
        data: [{main_js}],
        backgroundColor: [{colors_js}],
        borderAlign: 'inner',
        borderWidth: 1
      }}
    """

    datasets_code = f"[{ds_sub}, {ds_main}]"
    
    html = f"""
    <h4>{pop_name} ({country})</h4>
    <p>Age: {age}, Total: {total}</p>
    <p>{haplo_type}</p>
    <canvas id="{marker_id}" width="600" height="550"></canvas>

    <!-- Import Chart.js and DataLabels plugin -->
    <script src="{chart_js_cdn}"></script>
    <script src="{datalabels_cdn}"></script>
    <script>
      Chart.register(ChartDataLabels);
      // displaySub array stores percentages for display, parent group portions are null
      var displaySub = [{display_js}];

      var ctx = document.getElementById('{marker_id}').getContext('2d');
      var myChart = new Chart(ctx, {{
        type: 'doughnut',
        data: {{
          labels: [{labels_js}],
          datasets: {datasets_code}
        }},
        options: {{
          responsive: false,
          cutout: '50%',
          rotation: -0.5 * Math.PI,
          plugins: {{
            legend: {{ position: 'right' }},
            datalabels: {{
              display: function(context) {{
                var val = context.dataset.data[context.dataIndex];
                return val > 0;
              }},
              formatter: function(value, context) {{
                var dsIndex = context.datasetIndex;
                var lbl = context.chart.data.labels[context.dataIndex];
                if(dsIndex === 0) {{
                  // For outer ring subclades, display the percentage stored in displaySub (if exists)
                  var disp = displaySub[context.dataIndex];
                  return lbl + (disp !== null ? ' ' + disp.toFixed(1) + '%' : '');
                }} else {{
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
    for idx, row in df.iterrows():
        pop_name = row["Ancient pop name"]
        country  = row["Country"]
        age = row["Age"]
        lat = row["Lat"]
        lon = row["Long"]
        total = row.get("Total", 0)

        if pd.isna(lat) or pd.isna(lon):
            continue

        # Note: This now returns 4 arrays
        label_list, data_main, data_sub, display_sub = build_two_ring_data(row, haplo_cols)
        if not label_list:
            continue

        color_list = [color_map[lbl] for lbl in label_list]

        marker_id = f"{id_prefix}_marker_{idx}"
        popup_html = create_popup_html(
            marker_id=marker_id,
            pop_name=pop_name,
            country=country,
            age=age,
            total=total,
            label_list=label_list,
            data_main=data_main,
            data_sub=data_sub,
            display_sub=display_sub,
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

    # Collect all labels (parent groups and subclades) for unified color mapping
    all_labels = sorted(set(haplo_cols_y + haplo_cols_mt))
    color_map = generate_color_palette(all_labels)
    
    # Calculate map center point
    all_lat = pd.concat([df_y["Lat"], df_mt["Lat"]])
    all_lon = pd.concat([df_y["Long"], df_mt["Long"]])
    center_lat = all_lat.mean() if not all_lat.empty else 20
    center_lon = all_lon.mean() if not all_lon.empty else 0
    
    # Initialize map
    folium_map = folium.Map(location=[center_lat, center_lon], zoom_start=3)
    
    # Create layers
    fg_y  = folium.FeatureGroup(name="Y-chr")
    fg_mt = folium.FeatureGroup(name="mtDNA")

    # Add markers to each layer
    add_markers_to_group(df_y,  haplo_cols_y,  fg_y,  "y",  "Y-chr haplogroup",  color_map)
    add_markers_to_group(df_mt, haplo_cols_mt, fg_mt, "mt", "mtDNA haplogroup", color_map)
    
    folium_map.add_child(fg_y)
    folium_map.add_child(fg_mt)
    folium.LayerControl().add_to(folium_map)
    
    folium_map.save(args.output_html)
    print(f"The map has been saved to {args.output_html}")

if __name__ == "__main__":
    main()