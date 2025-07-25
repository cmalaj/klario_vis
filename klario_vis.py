
import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import re

st.set_page_config(layout="wide")
st.title("ClarioSTAR Growth Curve Visualisation Portal (Enhanced)")

# Generate distinct colours for 96 wells
rainbow_cmap = cm.get_cmap("gist_rainbow", 96)
well_order = [f"{row}{col:02}" for row in "ABCDEFGH" for col in range(1, 13)]
well_colours = {well: mcolors.to_hex(rainbow_cmap(i)) for i, well in enumerate(well_order)}

uploaded_files = st.file_uploader("Upload up to 4 CSV files at any one time.", type="csv", accept_multiple_files=True)

# Helper: Convert time like '1h 15 min' -> minutes
def convert_to_minutes(label):
    match = re.match(r'(?:(\d+)\s*h)?\s*(?:(\d+)\s*min)?', label.strip())
    if match:
        hours = int(match.group(1)) if match.group(1) else 0
        minutes = int(match.group(2)) if match.group(2) else 0
        return hours * 60 + minutes
    return None

# Helper: Convert minutes -> hh:mm
def minutes_to_hhmm(mins):
    h = int(mins // 60)
    m = int(mins % 60)
    return f"{h}:{m:02}"

# Layout generator
def generate_preset_layout(strain, phages):
    rows = list("ABCDEFGH")
    cols = [str(c).zfill(2) for c in range(1, 13)]
    layout_df = pd.DataFrame("", index=rows, columns=cols)

    tech_reps = ["T1", "T2"]
    batches = ["B1", "B2", "B3"]

    for row_idx, row_letter in enumerate(rows):
        phage_id = phages[row_idx // 2]
        tech_rep = tech_reps[row_idx % 2]

        well_values = []
        for moi in ["MOI1", "MOI0.5", "MOI0.1"]:
            for batch in batches:
                label = f"{phage_id}_{moi}-{strain}_{batch}-{tech_rep}"
                well_values.append(label)

        # Columns 10–12
        extras = [
            [phage_id, "BROTH", f"{strain}_B1"],
            [phage_id, "VEHICLE", f"{strain}_B1"],
            [phage_id, "PAO1", "EMPTY"],
            [phage_id, "EMPTY", f"{strain}_B2"],
            [phage_id, "BROTH", f"{strain}_B2"],
            [phage_id, "VEHICLE", "EMPTY"],
            [phage_id, "PAO1", f"{strain}_B3"],
            [phage_id, "EMPTY", f"{strain}_B3"]
        ]
        well_values += extras[row_idx]
        layout_df.loc[row_letter, :] = well_values

    return layout_df

# Parse CSV
def parse_growth_csv(file, plate_num):
    lines = file.getvalue().decode('utf-8').splitlines()
    header_idx = next(i for i, line in enumerate(lines) if line.startswith("Well,Content"))
    time_header_line = lines[header_idx + 1].split(',')[2:]
    timepoints = [convert_to_minutes(t) for t in time_header_line]
    data_lines = lines[header_idx + 2:]
    raw_data = [row.split(',')[2:] for row in data_lines]
    well_ids = [row.split(',')[0] for row in data_lines]
    df_raw = pd.DataFrame(raw_data, index=well_ids, columns=timepoints).apply(pd.to_numeric, errors='coerce')
    df = df_raw.transpose()
    df.index.name = "Time"
    df["Plate"] = f"Plate {plate_num}"
    return df

# Sidebar controls
st.sidebar.header("Global Plot Controls")
time_unit = st.sidebar.radio("X-axis units", options=["minutes", "hours", "hh:mm"], index=0)

all_rows = list("ABCDEFGH")
all_cols = list(range(1, 13))
select_all_rows = st.sidebar.checkbox("Select all rows", value=True)
selected_rows = all_rows if select_all_rows else st.sidebar.multiselect("Choose rows (A–H):", all_rows, default=all_rows)
select_all_cols = st.sidebar.checkbox("Select all columns", value=True)
selected_cols = all_cols if select_all_cols else st.sidebar.multiselect("Choose columns (1–12):", all_cols, default=all_cols)

# Remove global layout toggle — we now control layout per file directly in main UI

# Main app logic
if uploaded_files:
    all_data, all_summary, all_labels = [], [], []

    for i, file in enumerate(uploaded_files):
        df = parse_growth_csv(file, i + 1)
        if df.empty:
            st.warning(f"{file.name} could not be parsed.")
            continue

        plate = df["Plate"].iloc[0]
        st.subheader(f"{plate} - Time Series")
        df_plot = df.drop(columns="Plate")

        # Auto layout
        layout_map = {}
        if use_autolabels and len(phages) >= 4:
            layout_df = generate_preset_layout(strain, phages)
            layout_map = {f"{r}{int(c):02}": layout_df.loc[r, c] for r in layout_df.index for c in layout_df.columns}
        else:
            for row in selected_rows:
                for col in selected_cols:
                    well = f"{row}{col:02}"
                    label_key = f"label_{plate}_{well}"
                    layout_map[well] = st.sidebar.text_input(f"{plate} - Label for {well}", value=well, key=label_key)

        # Axis settings
        with st.expander(f"Adjust axis for {plate}"):
            x1, x2 = st.columns(2)
            with x1:
                x_min = st.number_input(f"{plate} X min", value=float(df_plot.index.min()), step=1.0)
                x_max = st.number_input(f"{plate} X max", value=float(df_plot.index.max()), step=1.0)
            with x2:
                y_min = st.number_input(f"{plate} Y min", value=float(df_plot.min().min()), step=0.1)
                y_max = st.number_input(f"{plate} Y max", value=float(df_plot.max().max()), step=0.1)


        # Per-file layout controls
        use_autolabels = st.checkbox(f"Use auto layout for {file.name}?", key=f"autolabel_{plate}")
        if use_autolabels:
            strain = st.text_input(f"{plate} - Strain name", value="PAO1", key=f"strain_{plate}")
            phage_input = st.text_input(f"{plate} - Phage IDs (comma-separated)", value="P1,P2,P3,P4", key=f"phages_{plate}")
            phages = [p.strip() for p in phage_input.split(",") if p.strip()]
            if len(phages) >= 4:
                layout_df = generate_preset_layout(strain, phages)
                layout_map = {f"{r}{int(c):02}": layout_df.loc[r, c] for r in layout_df.index for c in layout_df.columns}
        else:
            for row in selected_rows:
                for col in selected_cols:
                    well = f"{row}{col:02}"
                    label_key = f"label_{plate}_{well}"
                    layout_map[well] = st.text_input(f"{plate} - Label for {well}", value=well, key=label_key)

        # Plot individual wells
        fig = go.Figure()
        for col in df_plot.columns:
            match = re.match(r"([A-H])(\d{1,2})", col)
            if not match:
                continue
            row, col_num = match.groups()
            col_num = int(col_num)
            if row not in selected_rows or col_num not in selected_cols:
                continue
            colour = well_colours.get(f"{row}{col_num:02}", "#888888")
            label = layout_map.get(f"{row}{col_num:02}", f"{row}{col_num:02}")
            time_vals = df_plot.index
            if time_unit == "hours":
                time_vals = df_plot.index / 60
            elif time_unit == "hh:mm":
                time_vals = [minutes_to_hhmm(t) for t in df_plot.index]

            fig.add_trace(go.Scatter(x=time_vals, y=df_plot[col], name=label, line=dict(color=colour)))

        fig.update_layout(
            xaxis_title=f"Time ({time_unit})",
            yaxis_title="OD600",
            legend_title="Well Label",
            title=file.name,
            xaxis=dict(range=[x_min, x_max]),
            yaxis=dict(range=[y_min, y_max])
        )
        st.plotly_chart(fig, use_container_width=True)

        df_melted = df_plot.reset_index().melt(id_vars=["Time"], var_name="Well", value_name="OD600")
        df_melted["Label"] = df_melted["Well"].map(layout_map)
        all_data.append(df_melted)
        all_summary.append(df_plot.agg(["mean", "std"], axis=1))
        all_labels.append(layout_map)

    # Comparison plot
    st.subheader("Comparison Plot (mean ± SD)")
    common_labels = set.intersection(*[set(m["Label"].unique()) for m in all_data])
    selected_label = st.selectbox("Select label to compare:", sorted(common_labels))

    fig_compare = go.Figure()
    for i, df_m in enumerate(all_data):
        sub = df_m[df_m["Label"] == selected_label]
        if sub.empty:
            continue
        df_grp = sub.groupby("Time")["OD600"].agg(["mean", "std"]).reset_index()
        x_vals = df_grp["Time"]
        if time_unit == "hours":
            x_vals = x_vals / 60
        elif time_unit == "hh:mm":
            x_vals = [minutes_to_hhmm(t) for t in df_grp["Time"]]

        fig_compare.add_trace(go.Scatter(
            x=x_vals,
            y=df_grp["mean"],
            mode="lines",
            name=f"Plate {i+1}",
            line=dict(width=2)
        ))
        fig_compare.add_trace(go.Scatter(
            x=x_vals.tolist() + x_vals[::-1].tolist(),
            y=(df_grp["mean"] + df_grp["std"]).tolist() + (df_grp["mean"] - df_grp["std"])[::-1].tolist(),
            fill="toself",
            name=f"Plate {i+1} ± SD",
            line=dict(width=0),
            opacity=0.3,
            showlegend=False
        ))

    fig_compare.update_layout(
        title=f"Comparison of {selected_label}",
        xaxis_title=f"Time ({time_unit})",
        yaxis_title="OD600"
    )
    st.plotly_chart(fig_compare, use_container_width=True)