import pandas as pd # sometimes pandas will need to be conda remove'd and pip uninstall'd to install python >=3.7 to get this script to work. pip install pandas again afterwards.
import plotly.express as px # Do not conda install plotly . Use this: pip install plotly==5.11.0
import plotly.graph_objs as go # Do not conda install plotly . Use this: pip install plotly==5.11.0
import pdb, pathlib


# The parameters I want to be able to control are color, size, 

class Plotter():
    def __init__(self, output_dir, linkage_groups):
        self.output_dir = output_dir
        self.linkage_group_list = linkage_groups
        pathlib.Path(self.output_dir).mkdir(parents=True, exist_ok=True) # build the file path with pathlib.Path
        color_map = {'Mbuna': 'purple', 'AC': 'limegreen', 'Shallow_Benthic': 'red', 'Deep_Benthic': 'blue', 'Rhamphochromis': 'brown', 'Diplotaxodon': 'orange', 'Utaka': 'darkgreen', 'Riverine': 'pink'}
        # The below color hexcodes match what's in the Malinksy paper, extracted from Illustrator using the PDF of the publication 
        malinksy_color_map = {'Mbuna': '#A020F0', 'AC': '#A2CD5A', 'Shallow_Benthic': '#FF6347', 'Deep_Benthic': '#4876FF', 'Rhamphochromis': '#8B4513', 'Diplotaxodon': '#FFA54F', 'Utaka': '#006400'}
        # project_ID_shape_map = {'MalinskyData': 'square', 'Streelman_McGrathData': 'diamond', 'BrainDiversity_s1': 'star', 'MC_males': 'circle', 'MC_females': 'circle-open'} # removed for now to exclude shapes when generating data for Patrick's grant.
        bionano_shape_map = {'No': 'circle', 'Yes': 'x'}

        for lg in self.linkage_group_list:
            print('GENERATING PCA FOR ' + lg)
            # calculate percent variance explained by pc1 and 2. Round to 2 decimals
            variance_df = pd.read_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.eigenval', header=None)
            pc1_variance = (variance_df.loc[0][0] / variance_df.sum())[0]*100
            pc2_variance = (variance_df.loc[1][0] / variance_df.sum())[0]*100
            pc1_variance = round(pc1_variance, 2)
            pc2_variance = round(pc2_variance, 2)

            eigen_df = pd.read_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_new_projection.sscore', sep='\t')
            eigen_df = eigen_df.rename(columns = {'#IID':'SampleID'})
            
            df_merged = pd.merge(eigen_df, self.df, on=['SampleID'])
            df_merged['Color'] = df_merged['Ecogroup_PTM'].map(malinksy_color_map)  # Map Ecogroup_PTM to the malinsky_color_map
            df_merged.loc[df_merged['BionanoData'] == 'Yes', 'Color'] = 'black'  # Override to black for BionanoData 'Yes'
            df_merged['Size'] = df_merged['BionanoData'].apply(lambda x: 15 if x == 'Yes' else 8)
            if lg.startswith('NC'):
                plot_title = list(self.linkage_group_map.keys())[list(self.linkage_group_map.values()).index(lg)]
            else:
                plot_title = lg
            # Do the plotting magic using the new 'Color' column
            
            fig = px.scatter(df_merged, x='PC1_AVG', y='PC2_AVG',
                            labels={
                                'PC1_AVG': 'PC1 ' + str(pc1_variance) + '%',
                                'PC2_AVG': 'PC2 ' + str(pc2_variance) + '%'
                            },
                            # symbol_map=bionano_shape_map, # if you have a shape_map, uncomment and add that here
                            title=plot_title, hover_data=['SampleID', 'Ecogroup_PTM', 'Organism', 'ProjectID_PTM'])

            # Set various properties like size and color and shape here
            # fig.update_traces(marker=dict(size=8)) # No need for selectors as all points will be circles. You can add a selector here to specify which points get the size applied

            # Update traces for Bionano samples (diamonds) to size 8
            # fig.update_traces(marker=dict(size=9), selector=dict(marker_symbol='x'))

            # Update traces for non-Bionano samples (circles) to size 5
            # fig.update_traces(marker=dict(size=5), selector=dict(marker_symbol='circle'))

            # Apply custom colors and marker shapes
            fig.update_traces(marker=dict(color=df_merged['Color'])) # maps color based on the Color column in df_merged
            fig.update_traces(marker=dict(symbol=df_merged['BionanoData'].map(bionano_shape_map))) # maps the X's for BionanoData samples 
            fig.update_traces(marker=dict(size=df_merged['Size'])) # maps the size we want for the points based on the "Size" column.
            fig.write_html(self.plotly_out + lg + '_PCA.html')
        
        
        
        pass