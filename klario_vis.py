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
st.title("ClarioSTAR Growth Curve Visualisation Portal (CSV-compatible)")

# Generate distinct colours for 96 wells
rainbow_cmap = cm.get_cmap("gist_rainbow", 96)
well_order = [f"{row}{col:02}" for row in "ABCDEFGH" for col in range(1, 13)]
well_colours = {well: mcolors.to_hex(rainbow_cmap(i)) for i, well in enumerate(well_order)}

uploaded_files = st.file_uploader("Upload up to 4 CSV files at any one time.", type="csv", accept_multiple_files=True)

# Time conversion helper
def convert_to_minutes(label):
    match = re.match(r'(?:(\d+)\s*h)?\s*(?:(\d+)\s*min)?', label.strip())
    if match:
        hours = int(match.group(1)) if match.group(1) else 0
        minutes = int(match.group(2)) if match.group(2) else 0
        return hours * 60 + minutes
    return None

# Parser for your CSV file structure
def parse_growth_csv(file, plate_num):
    lines = file.getvalue().decode('utf-8').splitlines()
    header_idx = next(i for i, line in enumerate(lines) if line.startswith("Well,Content"))
    time_header_line = lines[header_idx + 1].split(',')[2:]  # skip 'Well' and 'Content'

    timepoints = [convert_to_minutes(t) for t in time_header_line]
    data_lines = lines[header_idx + 2:]
    raw_data = [row.split(',')[2:] for row in data_lines]
    well_ids = [row.split(',')[0] for row in data_lines]

    df_raw = pd.DataFrame(raw_data, index=well_ids, columns=timepoints).apply(pd.to_numeric, errors='coerce')
    df = df_raw.transpose()
    df.index.name = "Time"
    df["Plate"] = f"Plate {plate_num}"
    return df

# Sidebar for row/column selection (always show)
st.sidebar.header("Time-Series Controls")
all_rows = list("ABCDEFGH")
all_cols = list(range(1, 13))

# ROW selection
st.sidebar.subheader("Rows")
select_all_rows = st.sidebar.checkbox("Select all rows", value=True)
if select_all_rows:
    selected_rows = all_rows
else:
    selected_rows = st.sidebar.multiselect("Choose rows (A–H):", all_rows, default=all_rows)

# COLUMN selection
st.sidebar.subheader("Columns")
select_all_cols = st.sidebar.checkbox("Select all columns", value=True)
if select_all_cols:
    selected_cols = all_cols
else:
    selected_cols = st.sidebar.multiselect("Choose columns (1–12):", all_cols, default=all_cols)

# Now proceed if files are uploaded
if uploaded_files:
    all_data = []
    all_summary = []

    for i, file in enumerate(uploaded_files):
        df = parse_growth_csv(file, i + 1)
        if df.empty:
            st.warning(f"The file **{file.name}** could not be processed. Skipping.")
            continue

        all_data.append(df)

        numeric_cols = df.columns.drop("Plate", errors="ignore")
        summary = pd.DataFrame({
            "Well": numeric_cols,
            "Mean": df[numeric_cols].mean(),
            "SD": df[numeric_cols].std()
        }).reset_index(drop=True)
        summary["Plate"] = f"Plate {i + 1}"
        all_summary.append(summary)

    for idx, df in enumerate(all_data):
        plate = df["Plate"].iloc[0]
        st.subheader(f"{plate} - Time Series")

        # Custom labels specific to this plate
        custom_labels = {}
        with st.sidebar.expander(f"Custom Labels for {plate}"):
            for row in selected_rows:
                for col_num in selected_cols:
                    well_id = f"{row}{col_num:02}"
                    label_key = f"label_{plate}_{well_id}"
                    custom_label = st.text_input(f"{plate} - Label for {well_id}", value=well_id, key=label_key)
                    custom_labels[well_id] = custom_label

        # Axis controls
        with st.expander(f"Adjust axis ranges for {plate}"):
            col1, col2 = st.columns(2)
            with col1:
                x_min = st.number_input(f"{plate} X min (minutes)", value=float(df.index.min()), step=1.0, key=f"{plate}_xmin")
                x_max = st.number_input(f"{plate} X max (minutes)", value=float(df.index.max()), step=1.0, key=f"{plate}_xmax")
            with col2:
                y_min = st.number_input(f"{plate} Y min (OD600)", value=float(df.drop(columns='Plate').min().min()), step=0.1, key=f"{plate}_ymin")
                y_max = st.number_input(f"{plate} Y max (OD600)", value=float(df.drop(columns='Plate').max().max()), step=0.1, key=f"{plate}_ymax")

        fig = go.Figure()

        for col in df.columns:
            if col == "Plate":
                continue
            match = re.match(r"([A-H])(\d{1,2})", col)
            if not match:
                continue
            row, col_num = match.groups()
            col_num = int(col_num)
            well_id = f"{row}{col_num:02}"
            if row not in selected_rows or col_num not in selected_cols:
                continue
            colour = well_colours.get(well_id, "#CCCCCC")
            label = custom_labels.get(well_id, well_id)
            fig.add_trace(go.Scatter(
                x=df.index,
                y=df[col],
                name=label,
                mode='lines',
                line=dict(color=colour)
            ))

        fig.update_layout(
            xaxis_title="Time (minutes)",
            yaxis_title="OD600",
            legend_title="Well Label",
            margin=dict(l=50, r=50, t=50, b=50),
            xaxis=dict(range=[x_min, x_max]),
            yaxis=dict(range=[y_min, y_max])
        )

        st.plotly_chart(fig, use_container_width=True)

    # Summary heatmaps
    metrics = ["Mean", "SD"]
    fig, axes = plt.subplots(len(metrics), len(all_summary), figsize=(5 * len(all_summary), 5 * len(metrics)))

    if len(all_summary) == 1:
        axes = np.array([[axes[0]], [axes[1]]])

    for j, metric in enumerate(metrics):
        for i, summary in enumerate(all_summary):
            plate = summary["Plate"].iloc[0]
            sub = summary[["Well", metric]]

            well_ids = sub["Well"].dropna().unique()
            rows = sorted(set([w[0] for w in well_ids if re.match(r"^[A-Z]\d+$", w)]))
            cols = sorted(set([int(re.search(r"\d+$", w).group()) for w in well_ids if re.match(r"^[A-Z]\d+$", w)]))

            heatmap = pd.DataFrame(index=rows, columns=cols, dtype=float)

            for _, row in sub.iterrows():
                match = re.match(r"([A-Z])(\d{1,2})", row["Well"])
                if match:
                    r, c = match.groups()
                    if r in heatmap.index and int(c) in heatmap.columns:
                        heatmap.loc[r, int(c)] = row[metric]

            heatmap.columns = heatmap.columns.astype(int)

            sns.heatmap(
                heatmap,
                ax=axes[j][i],
                cmap="rainbow_r",
                annot=False,
                cbar=True
            )
            axes[j][i].set_title(f"{plate} - {metric}")
            axes[j][i].set_xlabel("Column")
            axes[j][i].set_ylabel("Row")

    plt.tight_layout()
    st.subheader("Plate Summary Heatmaps")
    st.pyplot(fig)