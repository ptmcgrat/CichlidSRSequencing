import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
import os

# 1. Define File Paths 
temp_dir = '/home/mcgrath-lab/nkumar317/backup_pca/temp_full_backup/lg10_Deep_Benthic_Inversion/'
sscore_file = temp_dir + "lg10_Deep_Benthic_Inversion_new_projection_no_ld.sscore"
eigenval_file = temp_dir + "lg10_Deep_Benthic_Inversion_subset_yes_pca_generation_no_ld_pruning_generated_step_3.eigenval"
metadata_file = "SampleDatabase_v2.xlsx" 

# 2. Calculate Variance Explained from Eigenval file
df_eigenval = pd.read_csv(eigenval_file, header=None)
total_variance = df_eigenval[0].sum()
pc1_variance = round((df_eigenval.iloc[0, 0] / total_variance) * 100, 1)
pc2_variance = round((df_eigenval.iloc[1, 0] / total_variance) * 100, 1)

# 3. Load ONLY the sscore data 
df_proj = pd.read_csv(sscore_file, sep='\t')
df_proj = df_proj.rename(columns={'#IID': 'SampleID'})

# 4. Load Metadata and Merge
df_meta = pd.read_excel(metadata_file, sheet_name='SampleLevel')
df_merged = pd.merge(df_proj, df_meta, on=['SampleID'], how='inner')

# 5. Map Colors and Sizes 
malinksy_color_map = {
    'Mbuna': '#A020F0', 'AC': '#A2CD5A', 'Shallow_Benthic': '#FF6347', 
    'Deep_Benthic': '#4876FF', 'Rhamphochromis': '#8B4513', 
    'Diplotaxodon': '#FFA54F', 'Utaka': '#006400'
}
bionano_shape_map = {'No': 'circle', 'Yes': 'x'}

df_merged['Color'] = df_merged['Ecogroup_PTM'].map(malinksy_color_map)
df_merged.loc[df_merged['BionanoData'] == 'Yes', 'Color'] = 'black'

# INCREASED SIZES HERE: Base points are now 7, Bionano background X's are 9
df_merged['HTML_Size'] = df_merged['BionanoData'].apply(lambda x: 9 if x == 'Yes' else 7)
df_merged['Opacity'] = df_merged['BionanoData'].apply(lambda x: 0.6 if x == 'Yes' else 0.9)

# 6. Generate Base Plot 
fig = px.scatter(
    df_merged, 
    x='PC1_AVG', 
    y='PC2_AVG', 
    hover_data=['SampleID', 'Ecogroup_PTM', 'Organism', 'ProjectID_PTM'],
    title="LG10 Inversion PCA",
    labels={
        'PC1_AVG': f'PC1 {pc1_variance}%',
        'PC2_AVG': f'PC2 {pc2_variance}%'
    },
    width=1100, 
    height=800
)

# Reverse the Y-axis
fig.update_yaxes(autorange="reversed")

# Apply marker updates
fig.update_traces(marker=dict(
    color=df_merged['Color'],
    symbol=df_merged['BionanoData'].map(bionano_shape_map),
    size=df_merged['HTML_Size'],
    opacity=df_merged['Opacity'],
    line_width=0
))

# 7. Isolate Target Bionano Samples for Highlighting
df_bionano_targets = df_merged[
    (df_merged['BionanoData'] == 'Yes') & 
    (df_merged['Organism'].str.contains('Mchenga|Aulonocara', case=False, na=False))
]

# 8. Add the Highlight Layer for the differential haplotypes
fig.add_trace(
    go.Scatter(
        x=df_bionano_targets['PC1_AVG'],
        y=df_bionano_targets['PC2_AVG'],
        mode='markers',
        marker=dict(
            symbol='x',
            color='black',
            size=18, # INCREASED SIZE HERE to maintain visual hierarchy 
            line=dict(width=3, color='black') 
        ),
        name='Target Bionano (Aulonocara/Mchenga)',
        text=df_bionano_targets['SampleID'] + "<br>" + df_bionano_targets['Organism'],
        hoverinfo='text'
    )
)

# 9. Export HTML and PDF
output_dir = "/home/mcgrath-lab/nkumar317/backup_pca/pca_out/"
os.makedirs(output_dir, exist_ok=True)

output_html = os.path.join(output_dir, "LG10_Inversion_Bionano_Highlighted_PCA.html")
output_pdf = os.path.join(output_dir, "LG10_Inversion_Bionano_Highlighted_PCA.pdf")

# Write HTML
fig.write_html(output_html)
print(f"Success! Interactive HTML generated and saved as {output_html}")

# Write PDF (Vector format for Illustrator)
fig.write_image(output_pdf)
print(f"Success! Vector PDF generated and saved as {output_pdf}")