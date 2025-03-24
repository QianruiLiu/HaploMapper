#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functionality:
  1) Reads frequency tables for Y-chromosome and mtDNA haplogroups (e.g., Y_haplogroup_frequencies.tsv and mtDNA_haplogroup_frequencies.tsv).
  2) Processes the input files to extract population metadata and haplogroup frequency data.
  3) Generates an interactive geographic map using Folium, displaying markers for ancient populations.
  4) Each marker features an interactive popup containing a double-ring pie chart that illustrates:
       - Basal haplogroup frequencies (inner ring)
       - Subclade frequencies (outer ring), with subclade percentages normalized relative to their corresponding basal group.
  5) Adds interactive filtering controls for BP (Before Present) intervals and haplogroup types (Y-chr/mtDNA) directly on the map.

Notes:
  - This script assumes that the input frequency tables include the following columns:
       "Ancient pop name", "Country", "Age", "Lat", "Long", and haplogroup frequency columns.
  - The interactive popups use Chart.js and the Chart.js datalabels plugin (loaded via CDN) for dynamic chart rendering.
  - Markers are organized into two separate feature groups for Y-chromosome and mtDNA data.
  - A control panel is added to allow users to filter markers based on chronological BP intervals and haplogroup type.

Usage example:
  python perfect_geography.py --y_input Y_haplogroup_frequencies.tsv --mt_input mtDNA_haplogroup_frequencies.tsv [--sep "\t"] [--output_html haplogroup_map.html]
"""

import argparse
import random
import re
import colorsys
import pandas as pd
import folium
from folium import IFrame
from collections import defaultdict
import sys
import os

#############################
# 1) Preserve original functions: color generation, loading tables, parent group/subclade identification
#############################

def generate_color_palette(labels):
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
    try:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Input file not found: {path}")
            
        df = pd.read_csv(path, sep=sep, dtype=str).fillna("")
        print(f"Loaded {len(df)} rows from {path}.")

        meta_cols = ["Ancient pop name", "Country", "Age", "Lat", "Long"]
        all_cols = list(df.columns)
        
        # Check if required columns exist
        missing_cols = [col for col in meta_cols if col not in all_cols]
        if missing_cols:
            raise ValueError(f"Missing required columns in {path}: {', '.join(missing_cols)}")
            
        haplo_cols = [c for c in all_cols if c not in meta_cols and c != "Total"]

        for c in haplo_cols:
            df[c] = df[c].str.replace("%", "", regex=False)
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)

        if "Total" in df.columns:
            df["Total"] = pd.to_numeric(df["Total"], errors="coerce").fillna(0).astype(int)

        df["Lat"] = pd.to_numeric(df["Lat"], errors="coerce")
        df["Long"] = pd.to_numeric(df["Long"], errors="coerce")
        
        # Check if we have valid coordinates
        if df["Lat"].isna().all() or df["Long"].isna().all():
            print(f"Warning: No valid coordinates found in {path}")
            
        return df, haplo_cols
    except pd.errors.EmptyDataError:
        raise ValueError(f"The file {path} is empty or contains no valid data")
    except pd.errors.ParserError:
        raise ValueError(f"Unable to parse {path}. Please check the file format and separator (current: '{sep}')")
    except Exception as e:
        raise Exception(f"Error processing {path}: {str(e)}")

def is_basal_col(col):
    return len(col) == 1 and col.isalpha()

def is_subclade_col(col):
    if len(col) < 2:
        return False
    return col[0].isalpha() and col[1].isdigit()

#############################
# 2) Modify the inner and outer ring functions to solve alignment issues
#############################
def build_two_ring_data(row, haplo_cols):
    """
    Generate unified label arrays and two-level data based on input row:
      label_list : Ordered, with parent groups immediately followed by their subclades
      data_main  : Inner ring data, parent groups use actual percentages, subclades are 0
      data_sub   : Outer ring data, parent groups are 0, subclades are normalized values
      display_sub: For displaying outer ring percentages, parent groups are None, subclades show relative percentages
      is_placeholder: Marks which items are placeholders
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
            parent = col[0]  # Parent group of A1 is A
            subclade_values.setdefault(parent, []).append((col, val_pct))

    basal_sorted = sorted(basal_values.items(), key=lambda x: x[1], reverse=True)

    label_list = []
    data_main = []
    data_sub = []
    display_sub = []
    is_placeholder = []

    for basal, basal_val in basal_sorted:
        # Add parent group: inner ring shows actual value, outer ring is 0
        label_list.append(basal)
        data_main.append(basal_val)
        data_sub.append(0)
        display_sub.append(None)
        is_placeholder.append(False)

        # Check if there are subclades
        if basal in subclade_values and subclade_values[basal]:
            subclades = subclade_values[basal]
            sum_sub = sum([val for (_, val) in subclades])
            for sc, val in sorted(subclades, key=lambda x: x[1], reverse=True):
                label_list.append(sc)
                data_main.append(0)  # Inner ring subclade position is 0
                if sum_sub > 0:
                    rel_val = (val / sum_sub) * basal_val
                    disp_val = (val / sum_sub) * 100
                else:
                    rel_val = 0
                    disp_val = 0
                data_sub.append(rel_val)
                display_sub.append(disp_val)
                is_placeholder.append(False)
        else:
            # No subclades, add an "empty" placeholder subclade to maintain inner/outer ring alignment
            placeholder_label = f"{basal}_placeholder"
            label_list.append(placeholder_label)
            data_main.append(0)
            data_sub.append(basal_val)
            display_sub.append(None)
            is_placeholder.append(True)

    # Handle subclades without corresponding parent groups
    leftover_parents = set(
        [subclade[0][0] for subclade in subclade_values.items() if subclade[0] not in basal_values]
    )
    for parent in leftover_parents:
        # For subclades without parent groups, first add a virtual parent group
        virtual_parent_label = f"{parent}"
        label_list.append(virtual_parent_label)
        data_main.append(0)
        data_sub.append(0)
        display_sub.append(None)
        is_placeholder.append(False)

        subclades = subclade_values[parent]
        sum_sub = sum([val for (_, val) in subclades])

        for sc, val in sorted(subclades, key=lambda x: x[1], reverse=True):
            label_list.append(sc)
            data_main.append(0)
            rel_val = val
            disp_val = (val / sum_sub) * 100 if sum_sub > 0 else 0
            data_sub.append(rel_val)
            display_sub.append(disp_val)
            is_placeholder.append(False)

    return label_list, data_main, data_sub, display_sub, is_placeholder

#############################
# 3) Update create_popup_html() function to handle transparent placeholders
#############################
def create_popup_html(chart_id, pop_name, country, age, total,
                      label_list, data_main, data_sub, display_sub, color_list,
                      haplo_type, is_placeholder):
    chart_js_cdn = "https://cdn.jsdelivr.net/npm/chart.js"
    datalabels_cdn = "https://cdn.jsdelivr.net/npm/chartjs-plugin-datalabels@2.0.0"

    labels_js = ", ".join([f"'{lbl}'" for lbl in label_list])
    main_js = ", ".join([f"{v:.2f}" for v in data_main])
    sub_js = ", ".join([f"{v:.2f}" for v in data_sub])
    colors_js = ", ".join([f"'{c}'" for c in color_list])
    display_js = ", ".join([f"{(v if v is not None else 'null')}" for v in display_sub])
    placeholder_js = ", ".join([f"{'true' if p else 'false'}" for p in is_placeholder])

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

    # Note: Here all curly braces in JS are written as '{{' or '}}'
    html = f"""
    <h4>{pop_name} ({country})</h4>
    <p>Age: {age}, Total: {total}</p>
    <p>{haplo_type}</p>
    <canvas id="{chart_id}" width="600" height="550"></canvas>

    <script src="{chart_js_cdn}"></script>
    <script src="{datalabels_cdn}"></script>
    <script>
      Chart.register(ChartDataLabels);
      // displaySub array stores percentages for display, parent groups are null
      var displaySub = [{display_js}];
      // isPlaceholder array marks which labels are placeholders
      var isPlaceholder = [{placeholder_js}];
      
      var ctx = document.getElementById('{chart_id}').getContext('2d');
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
            legend: {{ 
              position: 'right',
              labels: {{
                filter: function(legendItem, chartData) {{
                  // Filter out placeholder labels, don't display in legend
                  var index = legendItem.index;
                  return !isPlaceholder[index];
                }}
              }}
            }},
            tooltip: {{
              callbacks: {{
                label: function(context) {{
                  var index = context.dataIndex;
                  // Don't show tooltips for placeholders
                  if (isPlaceholder[index]) {{
                    return null;
                  }}
                  var label = context.label || '';
                  var value = context.raw || 0;
                  return label + ': ' + value.toFixed(2) + '%';
                }}
              }}
            }},
            datalabels: {{
              display: function(context) {{
                var val = context.dataset.data[context.dataIndex];
                var idx = context.dataIndex;
                // Don't show labels for placeholders
                return val > 0 && !isPlaceholder[idx];
              }},
              formatter: function(value, context) {{
                var dsIndex = context.datasetIndex;
                var lbl = context.chart.data.labels[context.dataIndex];
                var idx = context.dataIndex;
                
                // Skip placeholders
                if (isPlaceholder[idx]) {{
                  return null;
                }}
                
                if(dsIndex === 0) {{
                  // For outer ring subclades, show percentages stored in displaySub (if exists)
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

#############################
# 4) BP parsing function: extract numbers from row["Age"]
#############################
def parse_historical_age(age_str):
    if not isinstance(age_str, str):
        return 0
    
    is_bce = "BCE" in age_str or "BC" in age_str
    match = re.search(r"(\d+)", age_str)
    if not match:
        return 0
    
    year = int(match.group(1))
    if is_bce:
        return -year
    else:
        return year

#############################
# 5) Main function
#############################
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--y_input", required=True)
    parser.add_argument("--mt_input", required=True)
    parser.add_argument("--sep", default="\t")
    parser.add_argument("--output_html", default="haplogroup_map.html")
    args = parser.parse_args()

    try:
        df_y, haplo_cols_y = process_dataframe(args.y_input, args.sep)
        df_mt, haplo_cols_mt = process_dataframe(args.mt_input, args.sep)

        # Check if we have any data
        if len(df_y) == 0 and len(df_mt) == 0:
            raise ValueError("No data found in input files")

        # Get all labels
        all_haplo_cols = set(haplo_cols_y + haplo_cols_mt)
        if not all_haplo_cols:
            raise ValueError("No haplogroup columns found in the input files")
            
        placeholder_labels = [f"{col}_placeholder" for col in all_haplo_cols if is_basal_col(col)]
        all_labels = sorted(set(list(all_haplo_cols) + placeholder_labels))
        color_map = generate_color_palette(all_labels)

        all_lat = pd.concat([df_y["Lat"], df_mt["Lat"]])
        all_lon = pd.concat([df_y["Long"], df_mt["Long"]])
        
        # Handle case with no valid coordinates
        if all_lat.isna().all() or all_lon.isna().all():
            print("Warning: No valid coordinates found in either input file. Using default map center.")
            center_lat, center_lon = 20, 0
        else:
            center_lat = all_lat.mean() if not all_lat.empty else 20
            center_lon = all_lon.mean() if not all_lon.empty else 0

        folium_map = folium.Map(location=[center_lat, center_lon], zoom_start=3)

        fg_y = folium.FeatureGroup(name="Y-chr", overlay=True, control=True)
        fg_mt = folium.FeatureGroup(name="mtDNA", overlay=True, control=True)
        folium_map.add_child(fg_y)
        folium_map.add_child(fg_mt)

        folium.LayerControl(position='topright').add_to(folium_map)

        marker_count = 0
        markers_js = []
        y_markers = []
        mt_markers = []

        def add_marker(row, haplo_cols, feature_group, haplo_type):
            nonlocal marker_count
            try:
                pop_name = row["Ancient pop name"]
                country = row["Country"]
                age_str = row["Age"]
                lat = row["Lat"]
                lon = row["Long"]
                total = row.get("Total", 0)

                if pd.isna(lat) or pd.isna(lon):
                    return

                bp_val = parse_historical_age(age_str)

                label_list, data_main, data_sub, display_sub, is_placeholder = build_two_ring_data(row, haplo_cols)
                if not label_list:
                    return

                color_list = []
                for i, lbl in enumerate(label_list):
                    if is_placeholder[i]:
                        color_list.append('rgba(0, 0, 0, 0)')
                    else:
                        color_list.append(color_map.get(lbl, 'rgba(128, 128, 128, 0.75)'))  # Fallback color if missing

                chart_id = f"chart_{marker_count}"
                popup_html = create_popup_html(
                    chart_id=chart_id,
                    pop_name=pop_name,
                    country=country,
                    age=age_str,
                    total=total,
                    label_list=label_list,
                    data_main=data_main,
                    data_sub=data_sub,
                    display_sub=display_sub,
                    color_list=color_list,
                    haplo_type=haplo_type,
                    is_placeholder=is_placeholder
                )

                iframe = IFrame(html=popup_html, width=750, height=800)
                popup = folium.Popup(iframe, max_width=850)

                m = folium.Marker(
                    location=[lat, lon],
                    tooltip=f"{pop_name} ({country})",
                    popup=popup,
                )
                m.add_to(feature_group)

                varname = m.get_name()  # e.g. "marker_ABCDEFGH"
                # Note: here curly braces need to be {{ }}
                markers_js.append(f"markersData.push({{ bp: {bp_val}, markerVar: '{varname}', type: '{haplo_type}' }});")
                
                if haplo_type == "Y-chr haplogroup":
                    y_markers.append(varname)
                else:
                    mt_markers.append(varname)

                marker_count += 1
            except Exception as e:
                print(f"Warning: Failed to add marker for {row.get('Ancient pop name', 'unknown')}: {str(e)}")

        # Add markers to both FeatureGroups
        for idx, row in df_y.iterrows():
            add_marker(row, haplo_cols_y, fg_y, "Y-chr haplogroup")

        for idx, row in df_mt.iterrows():
            add_marker(row, haplo_cols_mt, fg_mt, "mtDNA haplogroup")

        if marker_count == 0:
            raise ValueError("No valid markers could be created from the input data")

        map_var = folium_map.get_name()

        # BP filtering functionality, using comma-separated year ranges, needs double curly braces
        bp_filter_html = f"""
        <div id="bpFilterControl" style="position: absolute; top: 10px; left: 10px; z-index:9999; background: white; padding: 10px; border-radius: 4px; box-shadow: 0 1px 5px rgba(0,0,0,0.4);">
          <label for="bpDropdown"><strong>Select Bp Interval:</strong></label><br>
          <select id="bpDropdown" style="width: 100%; margin-top: 5px; padding: 4px;">
            <option value="all">All Ages</option>
            <option value="-100000,-7000">Before 7000 BCE</option>
            <option value="-7000,-6000">7000-6000 BCE</option>
            <option value="-6000,-5000">6000-5000 BCE</option>
            <option value="-5000,-4000">5000-4000 BCE</option>
            <option value="-4000,-3000">4000-3000 BCE</option>
            <option value="-3000,-2500">3000-2500 BCE</option>
            <option value="-2500,-2000">2500-2000 BCE</option>
            <option value="-2000,-1500">2000-1500 BCE</option>
            <option value="-1500,-1000">1500-1000 BCE</option>
            <option value="-1000,-500">1000-500 BCE</option>
            <option value="-500,0">500-0 BCE</option>
            <option value="0,500">0-500 CE</option>
            <option value="500,1000">500-1000 CE</option>
            <option value="1000,1500">1000-1500 CE</option>
            <option value="1500,2000">1500-2000 CE</option>
            <option value="2000,100000">After 2000 CE</option>
          </select>
          <div style="margin-top: 8px;">
            <input type="checkbox" id="showYChr" checked><label for="showYChr"> Y-chr</label>
            <input type="checkbox" id="showMtDNA" checked style="margin-left: 10px;"><label for="showMtDNA"> mtDNA</label>
          </div>
        </div>

        <script>
        // Store: {{ bp: <int>, markerVar: 'marker_xxxx', type: 'Y-chr/mtDNA' }}
        var markersData = [];
        var yMarkers = [{', '.join([f"'{m}'" for m in y_markers])}];
        var mtMarkers = [{', '.join([f"'{m}'" for m in mt_markers])}];
        var currentBpLow = null;
        var currentBpHigh = null;
        var showY = true;
        var showMt = true;

        // This function displays/hides markers based on current filter conditions
        function applyAllFilters() {{
          for(var i=0; i<markersData.length; i++) {{
            var info = markersData[i];
            var bpVal = info.bp;
            var mkVar = info.markerVar;
            var markerType = info.type;
            var markerObj = window[mkVar];
            if(!markerObj) continue;

            // Check if marker type should be displayed
            var showByType = true;
            if(markerType === "Y-chr haplogroup" && !showY) {{
              showByType = false;
            }} else if(markerType === "mtDNA haplogroup" && !showMt) {{
              showByType = false;
            }}

            // Check if BP value is within selected range
            var showByBp = true;
            if(currentBpLow !== null && currentBpHigh !== null) {{
              if(bpVal < currentBpLow || bpVal > currentBpHigh) {{
                showByBp = false;
              }}
            }}

            // Decide whether to show marker
            if(showByType && showByBp) {{
              if(!{map_var}.hasLayer(markerObj)) {{
                {map_var}.addLayer(markerObj);
              }}
            }} else {{
              if({map_var}.hasLayer(markerObj)) {{
                {map_var}.removeLayer(markerObj);
              }}
            }}
          }}
        }}

        // Initialize controls
        setTimeout(function() {{
          // BP dropdown change event handler
          var bpSel = document.getElementById('bpDropdown');
          bpSel.onchange = function() {{
            var val = this.value;
            if(val === 'all') {{
              currentBpLow = null;
              currentBpHigh = null;
            }} else {{
              var parts = val.split(',');
              currentBpLow = parseFloat(parts[0]);
              currentBpHigh = parseFloat(parts[1]);
            }}
            applyAllFilters();
          }};

          // Y-chromosome checkbox change event
          var yChk = document.getElementById('showYChr');
          yChk.onchange = function() {{
            showY = this.checked;
            applyAllFilters();
          }};

          // mtDNA checkbox change event
          var mtChk = document.getElementById('showMtDNA');
          mtChk.onchange = function() {{
            showMt = this.checked;
            applyAllFilters();
          }};

          // Ensure layer controls work with filtering
          var layerControls = document.querySelectorAll('.leaflet-control-layers-selector');
          layerControls.forEach(function(control) {{
            control.addEventListener('change', function() {{
              // When layer state changes, we allow layer control to override filter results
              // But all filter conditions will be reapplied at the next filtering
            }});
          }});
        }}, 500);
        </script>
        """

        markers_definition_js = "\n".join(markers_js)
        combined_control = bp_filter_html + f"<script>\n{markers_definition_js}\n</script>"

        from folium import Element
        folium_map.get_root().html.add_child(Element(combined_control))

        # Check if output directory exists
        output_dir = os.path.dirname(os.path.abspath(args.output_html))
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        folium_map.save(args.output_html)
        print(f"Saved to {args.output_html}")
        
    except FileNotFoundError as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {str(e)}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()