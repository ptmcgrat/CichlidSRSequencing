#!/Users/kmnike/anaconda3/envs/mcgrath/bin/python3
import pandas as pd, subprocess as sp, os
import argparse

parser = argparse.ArgumentParser(usage = 'This script will take in an excel or csv datasheet and filtering options to output a list of filtered samples to use in downstream analyses with VCF data')
parser.add_argument('metadata', help = 'csv or xlsx datasheet with sampleIDs and metadata information')
parser.add_argument('-f', '--filters', help = 'one or multiple eco group names for filtering the data', choices = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhampochromis', 'Diplotaxodon', 'Riverine', 'AC', 'Non_Riverine', 'All'], nargs = '*', default = 'All')
parser.add_argument('-o', '--output', help = 'name of output sample file', default = 'filtered_samples.csv')
args = parser.parse_args()

if args.filters == 'All':
    eg = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhampochromis', 'Diplotaxodon', 'Riverine', 'AC']

elif args.filters[0] == 'Non_Riverine':
    eg = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhampochromis', 'Diplotaxodon']

else:
    eg = args.filters

df = pd.read_csv(args.metadata)
df_filtered = df[df['Ecogroup'].isin(eg)]
df_filtered['metadata_id'] = df_filtered['SampleID'] + "_" + df_filtered['Ecogroup']
df_filtered[['SampleID', 'metadata_id']].to_csv('filtered_samples.csv')
only_samples = df_filtered['SampleID']
only_samples.to_csv('only_sample_names.csv', index=False, header=False)
