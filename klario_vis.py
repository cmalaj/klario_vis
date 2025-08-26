import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from io import StringIO
import re
import copy

st.set_page_config(layout="wide")
st.title("ClarioSTAR Visualisation Portal v. 2.0")

# Generate 96 distinct colours from the rainbow colormap
rainbow_cmap = cm.get_cmap("gist_rainbow", 96)
well_order = [f"{row}{col}".replace("0", "") for row in "ABCDEFGH" for col in range(1, 13)]
well_colours = {well: mcolors.to_hex(rainbow_cmap(i)) for i, well in enumerate(well_order)}

uploaded_files = st.file_uploader("Upload one or more ClarioSTAR CSV files", type=["csv", "CSV"], accept_multiple_files=True)

def time_to_minutes(t):
    h, m, s = map(int, t.split(":"))
    return h * 60 + m + s / 60

def parse_growth_file(file):
    import pandas as pd
    import numpy as np
    import re
    from io import StringIO

    # Read all lines
    lines = file.getvalue().decode("utf-8").splitlines()

    # Find header line where the actual OD data starts
    header_line_idx = None
    for i, line in enumerate(lines):
        if line.startswith("Well,Content,"):
            header_line_idx = i
            break
    if header_line_idx is None:
        raise ValueError("Could not find OD data header.")

    # Extract time header
    time_header_line = lines[header_line_idx + 1]
    time_labels = time_header_line.split(",")[2:]  # Skip "Well" and "Content"

    # Convert time labels to minutes and hours
    def time_str_to_minutes(s):
        h, m = 0, 0
        match = re.match(r"(\d+)\s*h(?:\s*(\d+)\s*min)?", s.strip())
        if match:
            h = int(match.group(1))
            m = int(match.group(2)) if match.group(2) else 0
        return h * 60 + m

    time_mins = [time_str_to_minutes(t) for t in time_labels]
    time_hours = [t / 60 for t in time_mins]

    # Read OD data into DataFrame
    data_block = "\n".join(lines[header_line_idx + 2:])
    df = pd.read_csv(StringIO(data_block), header=None)

    # Drop rows with empty wells or labels
    df = df.dropna(subset=[0, 1])

    # ✅ Normalise well labels (e.g. A01 → A1)
    df[0] = df[0].astype(str).str.upper()
    df[0] = df[0].str.replace(r"^([A-H])0?(\d)$", r"\1\2", regex=True)
    df = df[df[0].str.match(r"^[A-H]\d{1,2}$")]

    # Extract well data
    wells = df[0].values
    labels = df[1].values
    od_data = df.iloc[:, 2:].replace("OVRFLW", np.nan).astype(float).values

    return wells, labels, od_data, time_mins, time_hours
def generate_preset_layout(strain, phages):
    rows = list("ABCDEFGH")
    cols = [str(c) for c in range(1, 13)]
    layout_df = pd.DataFrame("", index=rows, columns=cols)

    tech_reps = ["T1", "T2"]
    batches = ["B1", "B2", "B3"]

    for row_idx, row_letter in enumerate(rows):
        phage_id = phages[row_idx // 2]
        tech_rep = tech_reps[row_idx % 2]

        # Columns 1–9 (MOI combinations)
        well_values = []
        for moi in ["MOI1", "MOI0.5", "MOI0.1"]:
            for batch in batches:
                label = f"{phage_id}_{moi}-{strain}_{batch}-{tech_rep}"
                well_values.append(label)

        # Columns 10–12 (specials)
        if row_idx == 0:
            well_values += [phage_id, "EMPTY", f"{strain}_B1"]
        elif row_idx == 1:
            well_values += [phage_id, "EMPTY", f"{strain}_B1"]
        elif row_idx == 2:
            well_values += [phage_id, "EMPTY", "EMPTY"]
        elif row_idx == 3:
            well_values += [phage_id, "EMPTY", f"{strain}_B2"]
        elif row_idx == 4:
            well_values += [phage_id, "EMPTY", f"{strain}_B2"]
        elif row_idx == 5:
            well_values += [phage_id, "EMPTY", "EMPTY"]
        elif row_idx == 6:
            well_values += [phage_id, "EMPTY", f"{strain}_B3"]
        elif row_idx == 7:
            well_values += [phage_id, "EMPTY", f"{strain}_B3"]

        layout_df.loc[row_letter, :] = well_values

    return layout_df

if uploaded_files:
    all_data = []
    all_layouts = {}  # Dict to store well label maps per plate
    all_summary = []

    for i, file in enumerate(uploaded_files):
        wells, labels, od_data, time_mins, time_hours = parse_growth_file(file)
        original_label_map = dict(zip(wells, labels))

        df = pd.DataFrame(od_data.T, index=time_mins, columns=wells)
        df.index.name = "Time (min)"
        df = df[~df.index.duplicated(keep="first")]  # Drop duplicate timepoints if any
        df["Plate"] = f"Plate {i + 1}"  # ✅ only metadata column
        plate_name = f"Plate {i + 1}"
        filename = file.name
        default_title = filename or plate_name
        custom_title = st.text_input(
            f"Custom Title for {plate_name}",
            value=default_title,
            key=f"title_{plate_name}"
        )


        st.markdown(f"---\n### {plate_name} Layout Settings")

        layout_mode = st.radio(
            f"Layout Mode for {plate_name}",
            ["Use preset layout", "Start with empty layout"],
            horizontal=True,
            key=f"layout_mode_{plate_name}"
        )

        host_strain = st.text_input(
            f"{plate_name} - Bacterial Host Strain", value="PAO1", key=f"strain_{plate_name}"
        )

        phage_input = st.text_input(
            f"{plate_name} - Phage(s) (comma-separated)", value="P1,P2,P3,P4", key=f"phages_{plate_name}"
        )

        phages = [p.strip() for p in phage_input.split(",") if p.strip()]
        well_label_map = {}

        if layout_mode == "Use preset layout" and len(phages) == 4:
            layout_df = generate_preset_layout(host_strain, phages)
            for row in layout_df.index:
                for col in layout_df.columns:
                    well = f"{row}{col}"
                    label = layout_df.loc[row, col]
                    well_label_map[well] = label

            # Optional: show layout table
            st.markdown(f"**{plate_name} - Auto-generated Layout Preview**")
            st.dataframe(layout_df, use_container_width=True)
        elif layout_mode == "Use preset layout":
            st.warning(f"{plate_name}: You must enter exactly 4 phages for the preset layout.")

        # Store layout
        all_layouts[plate_name] = {
            "custom_title": custom_title,
            "well_map": {**original_label_map, **well_label_map}
        }

        if df.empty:
            st.warning(f"The file **{file.name}** could not be processed (empty or invalid data). Skipping.")
            continue

        all_data.append(df)

        numeric_cols = df.columns.drop(["Plate"], errors="ignore")
        numeric_cols = [col for col in numeric_cols if not col.startswith("T°")]
        numeric_cols = [col for col in df.columns if col not in ["Well", "Label", "Plate"]]

        summary = pd.DataFrame({
            "Well": numeric_cols,
            "Mean": df[numeric_cols].mean(),
            "SD": df[numeric_cols].std()
        }).reset_index(drop=True)
        summary["Plate"] = f"Plate {i + 1}"
        all_summary.append(summary)

    # Sidebar: Time-series well selection controls
    st.sidebar.header("Time-Series Controls")

    all_rows = list("ABCDEFGH")
    all_cols = list(range(1, 13))

    # Row selector
    st.sidebar.subheader("Rows")
    select_all_rows = st.sidebar.checkbox("Select all rows", value=True)
    if select_all_rows:
        selected_rows = all_rows
    else:
        selected_rows = st.sidebar.multiselect("Choose rows (A–H):", all_rows, default=all_rows, key="row_select")

    # Column selector
    st.sidebar.subheader("Columns")
    select_all_cols = st.sidebar.checkbox("Select all columns", value=True)
    if select_all_cols:
        selected_cols = all_cols
    else:
        selected_cols = st.sidebar.multiselect("Choose columns (1–12):", all_cols, default=all_cols, key="col_select")

    # Per-plate visualisation
    for idx, df in enumerate(all_data):
        plate = df["Plate"].iloc[0]
        st.subheader(f"{plate} - Time Series")

        # ✅ Retrieve stored title
        plate_info = all_layouts.get(plate, {})
        layout_map = plate_info.get("well_map", {})
        custom_title = plate_info.get("custom_title", plate)

        # Custom well labels for this plate only
        custom_labels = {}
        layout_map = all_layouts.get(plate, {}).get("well_map", {})
        with st.sidebar.expander(f"Custom Labels for {plate}"):
            for row in selected_rows:
                for col_num in selected_cols:
                    well_id = f"{row}{col_num}"
                    default_label = layout_map.get(well_id, well_id)
                    label_key = f"{plate}_{well_id}_label"
                    label = st.text_input(f"{plate} - Label for {well_id}", value=default_label, key=label_key)
                    custom_labels[well_id] = label

        # Time unit toggle
        time_unit = st.radio(
            f"{plate} – X-axis time unit",
            options=["Minutes", "Hours"],
            horizontal=True,
            key=f"time_unit_{plate}"
        )

        # Axis range override UI
        with st.expander(f"Adjust axis ranges for {plate}"):
            col1, col2 = st.columns(2)
            with col1:
                x_min_raw = df.index.min()
                x_max_raw = df.index.max()
                x_min_default = x_min_raw if time_unit == "Minutes" else x_min_raw / 60
                x_max_default = x_max_raw if time_unit == "Minutes" else x_max_raw / 60

                x_min = st.number_input(f"{plate} X min ({time_unit})", value=float(x_min_default), step=0.1, key=f"{plate}_xmin")
                x_max = st.number_input(f"{plate} X max ({time_unit})", value=float(x_max_default), step=0.1, key=f"{plate}_xmax")
            with col2:
                y_min = st.number_input(f"{plate} Y min (OD600)", value=float(df.drop(columns='Plate', errors='ignore').min().min()), step=0.1, key=f"{plate}_ymin")
                y_max = st.number_input(f"{plate} Y max (OD600)", value=float(df.drop(columns='Plate', errors='ignore').max().max()), step=0.1, key=f"{plate}_ymax")


        # Blank correction UI
        with st.expander(f"Blank Correction for {plate}"):
            apply_blank = st.checkbox(f"Apply blank correction for {plate}", key=f"{plate}_blank_toggle")
            blank_well = st.selectbox(
                f"Select blank well for {plate}",
                options=[f"{r}{c}" for r in "ABCDEFGH" for c in range(1, 13)],
                index=95,  # Default to H12
                key=f"{plate}_blank_select"
            )
        if apply_blank and blank_well in df.columns:
            df_corrected = df.copy()
            blank_values = df[blank_well]
            well_cols = df.filter(regex=r"^[A-H]\d{1,2}$").columns
            df_corrected[well_cols] = df[well_cols].subtract(blank_values, axis=0)
            df = df_corrected


        group_replicates = st.checkbox(
            f"Group technical replicates for {plate}?",
            value=False,
            key=f"group_reps_{plate}"
        )
        # Build plot
        fig = go.Figure()

        if group_replicates:
            import re

            # Group by normalized label (e.g., strip '-T1', '-T2')
            label_to_wells = {}

            for col in df.columns:
                if col in ["Plate"] or col.startswith("T°"):
                    continue
                match = re.match(r"([A-H])(\d{1,2})", col)
                if not match:
                    continue
                row, col_num = match.groups()
                col_num = int(col_num)
                well_id = f"{row}{col_num}"
                if row not in selected_rows or col_num not in selected_cols:
                    continue

                # Get the label from layout or fallback
                label = custom_labels.get(well_id) or layout_map.get(well_id, well_id)

                # ✅ Strip trailing '-T1', '-T2', etc. for grouping
                group_label = re.sub(r"-T\d$", "", label)

                label_to_wells.setdefault(group_label, []).append(col)

            for group_label, replicate_cols in label_to_wells.items():
                if not replicate_cols:
                    continue
                # Use first well to determine color
                colour = well_colours.get(replicate_cols[0], "#CCCCCC")
                x_vals = df.index if time_unit == "Minutes" else df.index / 60
                values = df[replicate_cols].values  # shape: time x replicates
                mean_vals = np.nanmean(values, axis=1)
                std_vals = np.nanstd(values, axis=1)

                # Convert matplotlib RGBA to valid Plotly rgba string
                rgba = mcolors.to_rgba(colour, alpha=0.2)
                fillcolor = f"rgba({int(rgba[0]*255)}, {int(rgba[1]*255)}, {int(rgba[2]*255)}, {rgba[3]})"

                # Mean line
                fig.add_trace(go.Scatter(
                    x=x_vals,
                    y=mean_vals,
                    mode='lines',
                    name=group_label,
                    line=dict(color=colour, width=2),
                    legendgroup=group_label,
                    showlegend=True
                ))

                # SD ribbon
                fig.add_trace(go.Scatter(
                    x=np.concatenate([x_vals, x_vals[::-1]]),
                    y=np.concatenate([mean_vals + std_vals, (mean_vals - std_vals)[::-1]]),
                    fill='toself',
                    fillcolor=fillcolor,
                    line=dict(color='rgba(255,255,255,0)'),
                    hoverinfo="skip",
                    showlegend=False,
                    legendgroup=group_label
                ))

        else:
            # Fall back to individual wells
            for col in df.columns:
                if col in ["Plate"] or col.startswith("T°"):
                    continue
                match = re.match(r"([A-H])(\d{1,2})", col)
                if not match:
                    continue
                row, col_num = match.groups()
                col_num = int(col_num)
                well_id = f"{row}{col_num}"
                if row not in selected_rows or col_num not in selected_cols:
                    continue
                label = custom_labels.get(well_id) or layout_map.get(well_id, well_id)
                colour = well_colours.get(well_id, "#CCCCCC")
                x_vals = df.index if time_unit == "Minutes" else df.index / 60

                fig.add_trace(go.Scatter(
                    x=x_vals,
                    y=df[col],
                    name=label,
                    mode='lines',
                    line=dict(color=colour)
                ))

        # Final plot layout + render
        fig.update_layout(
            title=custom_title,
            xaxis_title=f"Time ({time_unit})",
            yaxis_title="OD600",
            legend_title="Well Label",
            margin=dict(l=50, r=50, t=50, b=50),
            xaxis=dict(range=[x_min, x_max]),
            yaxis=dict(range=[y_min, y_max])
        )

        st.plotly_chart(fig, use_container_width=True)


st.markdown("---")
st.header("Across-File Comparison")

show_comparison = st.checkbox("Enable Comparison Plot", value=False)

if show_comparison:
    compare_by = st.radio(
        "Compare by:",
        options=["Well location", "Custom label", "Manual selection"],
        horizontal=True
    )

    group_replicates = False
    if compare_by == "Custom label":
        group_replicates = st.checkbox(
            "Group technical replicates (same base label)?",
            value=True,
            help="Strips -T1/-T2 etc. and aggregates the matching profiles"
        )
    label_pool = set()
    well_pool = set()

    # Optional: Toggle for x-axis units
    time_unit_compare = st.radio(
        "X-axis time unit for comparison plot",
        options=["Minutes", "Hours"],
        horizontal=True,
        key="time_unit_compare"
    )

    for df in all_data:
        plate = df["Plate"].iloc[0]
        layout_map = all_layouts.get(plate, {}).get("well_map", {})
        for col in df.columns:
            if re.match(r"^[A-H]\d{1,2}$", col):
                label = st.session_state.get(f"{plate}_{col}_label", col)
                label_pool.add(label)
                well_pool.add(col)

    if compare_by == "Custom label":
        if group_replicates:
            base_labels = {re.sub(r"-T\d$", "", lbl) for lbl in label_pool}
            selection = st.multiselect(
                "Select sample labels to compare across plates",
                options=sorted(base_labels),
                key="compare_label_selector"
            )
        else:
            selection = st.multiselect(
                "Select sample labels to compare across plates",
                options=sorted(label_pool),
                key="compare_label_selector"
            )
    elif compare_by == "Well location":
        selection = st.multiselect(
            "Select wells to compare across plates",
            options=sorted(well_pool),
            key="compare_well_selector"
        )
    elif compare_by == "Manual selection":
        manual_selection = []
        st.markdown("**Manually select specific wells or labels from each file:**")
        for df in all_data:
            plate = df["Plate"].iloc[0]
            layout_map = all_layouts.get(plate, {}).get("well_map", {})

            st.markdown(f"**{plate}**")
            cols_in_plate = [col for col in df.columns if re.match(r"^[A-H]\d{1,2}$", col)]
            well_options = []
            for col in cols_in_plate:
                label = st.session_state.get(f"{plate}_{col}_label", layout_map.get(col, col))
                well_options.append(f"{plate}::{col}::{label}")

            selected = st.multiselect(
                f"Select wells or labels from {plate}",
                options=well_options,
                key=f"manual_select_{plate}"
            )
            manual_selection.extend(selected)
    
    with st.expander("Adjust axis ranges for comparison plot"):
        col1, col2 = st.columns(2)
        with col1:
            all_timepoints = np.concatenate([df.index.values for df in all_data])
            x_min_raw = all_timepoints.min()
            x_max_raw = all_timepoints.max()
            x_min_default = x_min_raw if time_unit_compare == "Minutes" else x_min_raw / 60
            x_max_default = x_max_raw if time_unit_compare == "Minutes" else x_max_raw / 60

            x_min = st.number_input("X min", value=float(x_min_default), step=0.1, key="comp_xmin")
            x_max = st.number_input("X max", value=float(x_max_default), step=0.1, key="comp_xmax")
        with col2:
            # Precompute reasonable y-limits from all selected wells
            all_y_values = []
            for df in all_data:
                numeric_cols = df.columns.drop(["Plate"], errors="ignore")
                for col in numeric_cols:
                    if re.match(r"^[A-H]\d{1,2}$", col):
                        all_y_values.append(df[col].values)
            if all_y_values:
                all_y_array = np.concatenate(all_y_values)
                y_min_default = float(np.nanmin(all_y_array))
                y_max_default = float(np.nanmax(all_y_array))
            else:
                y_min_default = 0.0
                y_max_default = 1.5

            y_min = st.number_input("Y min (OD600)", value=y_min_default, step=0.1, key="comp_ymin")
            y_max = st.number_input("Y max (OD600)", value=y_max_default, step=0.1, key="comp_ymax")
    
    dash_styles = ["solid", "dash", "dot", "dashdot", "longdash", "longdashdot"]
    fig = go.Figure()

    for plate_idx, df in enumerate(all_data):
        plate = df["Plate"].iloc[0]
        layout_map = all_layouts.get(plate, {}).get("well_map", {})

        matched = {}

        for col in df.columns:
            if not re.match(r"^[A-H]\d{1,2}$", col):
                continue

            label = st.session_state.get(f"{plate}_{col}_label", layout_map.get(col, col))

            if compare_by == "Well location" and col in selection:
                matched.setdefault(f"{plate}::{col}", []).append(col)

            elif compare_by == "Custom label":
                base_label = re.sub(r"-T\d$", "", label) if group_replicates else label
                if base_label in selection:
                    matched.setdefault(f"{plate}::{base_label}", []).append(col)

            elif compare_by == "Manual selection":
                for item in manual_selection:
                    p, w, l = item.split("::")
                    if p == plate and (w == col or l == label):
                        matched.setdefault(f"{plate}::{l}", []).append(col)

        for name, cols in matched.items():
            if not cols:
                continue
            colour = well_colours.get(cols[0], "#888888")
            x_vals = df.index if time_unit_compare == "Minutes" else df.index / 60
            values = df[cols].values
            mean_vals = np.nanmean(values, axis=1)
            std_vals = np.nanstd(values, axis=1)

            rgba = mcolors.to_rgba(colour, alpha=0.2)
            fillcolor = f"rgba({int(rgba[0]*255)}, {int(rgba[1]*255)}, {int(rgba[2]*255)}, {rgba[3]})"

            fig.add_trace(go.Scatter(
                x=x_vals,
                y=mean_vals,
                mode='lines',
                name=f"{plate} – {name}",
                line=dict(color=colour, dash=dash_styles[plate_idx % len(dash_styles)]),
                legendgroup=name
            ))
            fig.add_trace(go.Scatter(
                x=np.concatenate([x_vals, x_vals[::-1]]),
                y=np.concatenate([mean_vals + std_vals, (mean_vals - std_vals)[::-1]]),
                fill='toself',
                fillcolor=fillcolor,
                line=dict(color='rgba(255,255,255,0)'),
                showlegend=False,
                hoverinfo="skip",
                legendgroup=name
            ))

    fig.update_layout(
        title="Comparison Plot",
        xaxis_title=f"Time ({time_unit_compare})",
        yaxis_title="OD600",
        legend_title="Sample",
        margin=dict(l=50, r=50, t=50, b=50),
        xaxis=dict(range=[x_min, x_max]),
        yaxis=dict(range=[y_min, y_max])
    )
    st.plotly_chart(fig, use_container_width=True)

    # ----------------------
# Bacterial Growth Threshold Analysis
# ----------------------
st.markdown("---")
st.header("Growth Threshold Analysis")

thresholds = [1000, 100, 250, 500, 2000, 3000, 5000, 7500, 10000]

for idx, df in enumerate(all_data):
    plate = df["Plate"].iloc[0]
    st.subheader(f"{plate} – AUC Grid Analysis")

    # Get time axis
    time_vals = df.index.values
    time_hours = time_vals / 60

    # Identify candidate wells (standard 96-well IDs)
    candidate_wells = [col for col in df.columns if re.fullmatch(r"[A-H]1[0-2]?|[A-H][1-9]", col)]

    preferred_default_wells = ["A12", "B12", "D12", "E12", "G12", "H12"]
    default = [w for w in preferred_default_wells if w in candidate_wells]

    selected_positions = st.multiselect(
        f"Choose control well positions (e.g. A12, B11) for {plate}",
        options=candidate_wells,
        default=default,
        key=f"ctrl_wells_{plate}"
    )

    if not selected_positions:
        st.warning(f"No control well positions selected for {plate}. Skipping.")
        continue

    mean_vals = df[selected_positions].mean(axis=1)
    baseline = mean_vals.iloc[0]

    st.subheader(f"{plate} – Select Growth Threshold for AUC Integration")
    threshold_to_use = st.selectbox(
        f"Select threshold (× growth) for {plate}",
        thresholds,
        index=0,
        key=f"threshold_selector_{plate}"
    )

    # Threshold logic
    thresh_val = baseline * threshold_to_use
    cross_idx = np.argmax(mean_vals.values >= thresh_val)
    cross_time = time_hours[cross_idx] if cross_idx < len(time_hours) else None

    # --- Plot mean growth curve with CI ---
    fig = go.Figure()
    std_vals = df[selected_positions].std(axis=1)
    upper_bound = mean_vals + std_vals
    lower_bound = mean_vals - std_vals

    fig.add_trace(go.Scatter(
        x=time_hours,
        y=mean_vals,
        name="Mean Control Growth",
        line=dict(color="blue", width=3),
    ))
    fig.add_trace(go.Scatter(
        x=time_hours,
        y=lower_bound,
        line=dict(color="rgba(0,0,255,0)"),
        showlegend=False,
        hoverinfo='skip',
    ))
    fig.add_trace(go.Scatter(
        x=time_hours,
        y=upper_bound,
        fill='tonexty',
        fillcolor='rgba(0,0,255,0.2)',
        line=dict(color="rgba(0,0,255,0)"),
        name="±1 SD",
        hoverinfo='skip',
    ))

    if cross_time is not None:
        fig.add_shape(
            type="line",
            x0=cross_time, x1=cross_time,
            y0=0, y1=mean_vals.max() * 1.1,
            line=dict(dash="dash", color="red")
        )
        fig.add_trace(go.Scatter(
            x=[cross_time], y=[thresh_val],
            mode="markers+text",
            marker=dict(color="red", size=6),
            text=[f"{threshold_to_use}×"],
            textposition="top center",
            showlegend=False
        ))

    # --- Compute ΔAUC grid ---
    if cross_time is not None:
        valid_mask = time_hours <= cross_time
        control_auc = np.trapz(mean_vals[valid_mask], time_hours[valid_mask])

        delta_auc_grid = pd.DataFrame(index=list("ABCDEFGH"), columns=[str(i) for i in range(1, 13)])
        for well in candidate_wells:
            row, col = well[0], well[1:]
            if well not in df.columns or row not in delta_auc_grid.index or col not in delta_auc_grid.columns:
                continue
            curve = df[well]
            well_auc = np.trapz(curve[valid_mask], time_hours[valid_mask])
            delta_auc = well_auc - control_auc
            delta_auc_grid.loc[row, col] = delta_auc

        delta_auc_grid = delta_auc_grid.apply(pd.to_numeric)
        norm = mcolors.TwoSlopeNorm(vcenter=0, vmin=delta_auc_grid.min().min(), vmax=delta_auc_grid.max().max())
        cmap = cm.get_cmap("coolwarm_r")

        st.subheader(f"{plate} – ΔAUC Well Grid (up to {threshold_to_use}×)")

        for row in list("ABCDEFGH"):
            cols = st.columns(12, gap="small")
            for i, col_num in enumerate(range(1, 13)):
                well_id = f"{row}{col_num}"
                delta = delta_auc_grid.loc[row, str(col_num)]
                if pd.isna(delta):
                    cols[i].markdown(f"<div style='height:36px;'></div>", unsafe_allow_html=True)
                    continue
                colour = mcolors.to_hex(cmap(norm(delta)))
                delta_formatted = f"{delta:+.3f}"
                cols[i].markdown(
                    f"""
                    <div style="
                        background-color: {colour};
                        border: 1px solid #444;
                        border-radius: 3px;
                        height: 36px;
                        line-height: 1.2;
                        text-align: center;
                        font-size: 11px;
                        font-weight: bold;
                        color: black;
                        padding-top: 2px;
                    ">
                        {well_id}<br><span style='font-weight: normal'>{delta_formatted}</span>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )

    st.plotly_chart(fig, use_container_width=True)