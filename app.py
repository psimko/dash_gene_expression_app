#Imports
import os
import pandas as pd
import numpy as np
from numpy import inf
import matplotlib.pyplot as plt
import pickle
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import flask
import dash
from dash import dcc, html, Input, Output, State
import plotly.colors
import matplotlib.colors as mcolors
import boto3
from botocore.exceptions import ClientError
import gc
import psutil

# Call the s3 bucket
client = boto3.client('s3')

# Set the values
bucket_name = 'gene-app-mouse-data'
#key = '/data'

################################################################################################################

def print_memory_usage(tag=""):
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / 1024 / 1024  # In MB
    print(f"[Memory] {tag} - {mem:.2f} MB")

def read_file(bucket, key_value):
    try:
        s3 = boto3.client('s3')
        obj = s3.get_object(Bucket=bucket, Key=key_value)
        df = pd.read_csv(obj['Body'])
        return df
    except ClientError as ex:
        if ex.response['Error']['Code'] == 'NoSuchKey':
            print(f"[ERROR] No such key: {key_value}")
        else:
            print(f"[ERROR] AWS ClientError: {ex}")
    except Exception as e:
        print(f"[ERROR] Unexpected error reading CSV: {e}")
    
    return None  

def read_pickle(bucket, key_value):
    try:
        s3 = boto3.client('s3')
        obj = s3.get_object(Bucket=bucket, Key=key_value)
        data = pickle.load(obj['Body'])
        return data
    except ClientError as ex:
        if ex.response['Error']['Code'] == 'NoSuchKey':
            print("Key doesn't match. Please check the key value entered.")
    except Exception as e:
        print("Error during unpickling:", e)
    
# Write the data
# avg_expression_genesAll_div_df = read_file(bucket_name, '/data/avg_expression_div_genesAll_notNormalized_df.csv')
# avg_expression_genesAll_class_df = read_file(bucket_name, '/data/avg_expression_class_genesAll_notNormalized_df.csv') 
# avg_expression_genesAll_subclass_df =  read_file(bucket_name, '/data/avg_expression_subclass_genesAll_notNormalized_df.csv')    
# avg_expression_genesAll_supertype_df = read_file(bucket_name, '/data/avg_expression_supertypes_genesAll_notNormalized_df.csv')       
    
################################################################################################################ 
    

# Pick standard genes from expression_df and the average expression files

# If you want to test only genes from a list use this
gene_list_otherVendors = ['Acta2', 'Adgre1', 'Abcc8', 'Abcc9', 'Actl6b', 'Aif1', 'Akap5', 'Aldh5a1', 'App', 'Aqp1', 'Aqp4', 'Arg1', 'Bcl11b', 'Calca',
            'Ccnd1', 'Cd247', 'Cd3e', 'Cd4', 'Cd5', 'Cd68', 'Cd86', 'Cd8a', 'Cdh1', 'Chat', 'Cnp', 'Cntnap1', 'Cntnap2', 'Col4a3/5/2/1',
            'Creb1', 'Cspg4', 'Ctnnb1', 'Dbh', 'Dcx', 'Ddx5', 'Dlg2', 'Eea1', 'Eea5', 'Egr1', 'Emcn', 'Epm2a', 'Ewsr1', 'Fn1', 'Foxa2', 'Gad1', 'Gad2',
            'Gad2/1', 'Gap43', 'Gfap', 'Gria2', 'Grin1', 'Grm2', 'Gsk3a', 'Gsk3a/b', 'Gucy1b1', 'Hcls1', 'Hopx', 'Htr2b', 'Htr7', 'Il10', 'Ins',
            'Itgam', 'Itgax', 'Khdrbs1', 'Lamp1', 'Lyve1', 'Mag', 'Maoa', 'Maob', 'Map2', 'Mapk3', 'Mapk8/9/10', 'Mapt', 'Mbp', 'Mki67', 'Mog',
            'Mrc1', 'Myb', 'Ncam1', 'Nefh', 'Nefl', 'Nefm', 'Nefm/h', 'Nfasc', 'Nfatc1', 'Nos1', 'Nos3', 'Npy', 'Nr3c2', 'Nrp1', 'Ntrk3', 'Ocrl', 'Oxt',
            'P2rx4', 'P2ry12', 'Pax6', 'Pax7', 'Pdgfrb', 'Pecam1', 'Plp1', 'Ppp1r1b', 'Prkca/b/g', 'Pvalb', 'Pycard', 'Rbbp4', 'Rbfox3', 'S100a10',
            'S100b', 'Satb2', 'Sdc4', 'Sdk2', 'Set', 'Sirt3', 'Slc1a2', 'Slc1a3', 'Slc6a3', 'Slc6a4', 'Snca', 'Sod2', 'Sox2', 'Sox4', 'Sox9', 'Sst',
            'Stat1', 'Stx1a', 'Stx1a/1b/2/3', 'Sun2', 'Syn1', 'Syn2', 'Syp', 'Tardbp', 'Tbr1', 'Th', 'Tmem119', 'Tph1', 'Tph2', 'Tuba', 'Tubb', 'Tubb3',
            'Uchl1', 'Vim']

gene_list_neuromab = ['Adam11', 'Aldh1l1', 'Amigo1', 'Arx', 'Atp7a', 'Bdnf', 'Cacna1h', 'Cadm4', 'Calb1', 'Calb2', 'Clcn4', 'Cntnap1',
                     'Dlg1', 'Dlg2', 'Dlg3', 'Dlg4', 'Drd2', 'Fgf13', 'Gabrb3', 'Gabre', 'Hspa9', 'Kcna1', 'Kcnd2', 'Lrp4',
                     'Lrrk1', 'Lrrtm2', 'Lrrtm4', 'Mff', 'Mog', 'Nos1', 'Npy', 'Nrcam', 'Olfm1', 'Znf746', 'Pex5l', 'Qk',
                     'Rbm17', 'Reep1', 'Reep2', 'Rufy3', 'S100a5', 'Shank1', 'Shank2', 'Shank3', 'Slc38a1', 'Snapin', 'Svop', 'Trpc4',
                     'Vapa', 'Vdac1', 'Tpte']

excluded_genes = ['Col4a3/5/2/1', 'Eea5', 'Gad2/1', 'Gsk3a/b', 'Ins', 
                  'Mapk8/9/10', 'Nefm/h', 'Prkca/b/g', 'Stx1a/1b/2/3', 
                  'Tuba', 'Tubb', 'Znf746', 'Qki']

gene_list = [
    gene for gene in (gene_list_otherVendors + gene_list_neuromab) 
    if gene not in excluded_genes
]

avg_expression_genesAll_div_df = read_file(bucket_name, 'data/avg_expression_div_genesAll_notNormalized_df.csv')
avg_expression_genesAll_div_df = avg_expression_genesAll_div_df[gene_list]
gc.collect()
print_memory_usage("After loading avg_expression_genesAll_div_df")

avg_expression_genesAll_class_df = read_file(bucket_name, 'data/avg_expression_class_genesAll_notNormalized_df.csv') 
avg_expression_genesAll_class_df = avg_expression_genesAll_class_df[gene_list]
gc.collect()
print_memory_usage("After loading avg_expression_genesAll_class_df")

avg_expression_genesAll_subclass_df =  read_file(bucket_name, 'data/avg_expression_subclass_genesAll_notNormalized_df.csv')
avg_expression_genesAll_subclass_df = avg_expression_genesAll_subclass_df[gene_list]
gc.collect()
print_memory_usage("After loading avg_expression_genesAll_subclass_df")

avg_expression_genesAll_supertype_df = read_file(bucket_name, 'data/avg_expression_supertypes_genesAll_notNormalized_df.csv')  
avg_expression_genesAll_supertype_df = avg_expression_genesAll_supertype_df[gene_list]    
gc.collect()
print_memory_usage("After loading avg_expression_genesAll_supertype_df")
    
    
################################################################################################################
# Average expressions of genes in all taxonomy levels
# #expression_df = pd.read_csv('/bil/users/psimko/holis/transcriptomic_analysis/expression_df_types.csv', index_col=0)
# avg_expression_genesAll_div_df = pd.read_csv('/data/avg_expression_div_genesAll_notNormalized_df.csv', index_col=0)
# avg_expression_genesAll_class_df = pd.read_csv('/data/avg_expression_class_genesAll_notNormalized_df.csv', index_col=0)
# avg_expression_genesAll_subclass_df = pd.read_csv('/data/avg_expression_subclass_genesAll_notNormalized_df.csv', index_col=0)
# avg_expression_genesAll_supertype_df = pd.read_csv('/data/avg_expression_supertypes_genesAll_notNormalized_df.csv', index_col=0)

# Load binarized expressions
# with open('/bil/users/psimko/holis/clustering/2025_holis_analysis/expression_cells_bin_df.pickle', 'rb') as file:
#     expression_cells_bin_df_copy = pickle.load(file)
################################################################################################################

# Import taxonomy dictionaries if needed

sample_to_type = read_pickle(bucket_name, '/data/taxonomy_dictionaries/sample_to_type.pkl')

if sample_to_type is None:
    raise RuntimeError("Failed to load pickle — returned None.")
    
print_memory_usage("After loading sample_to_type")
    
sample_to_subclass = {
    sample: [type_to_subclass.get(sample_type, None) for sample_type in sample_types]
    if isinstance(sample_types, list) else type_to_subclass.get(sample_types, None)
    for sample, sample_types in sample_to_type.items()
}

del sample_to_type
gc.collect()

print_memory_usage("After loading sample_to_subclass and del sample_to_type.")

class_to_division = read_pickle(bucket_name, '/data/taxonomy_dictionaries/class_to_division.pkl')
division_to_class = read_pickle(bucket_name, '/data/taxonomy_dictionaries/division_to_class.pkl')
subclass_to_class = read_pickle(bucket_name, '/data/taxonomy_dictionaries/subclass_to_class.pkl')
class_to_subclass = read_pickle(bucket_name, '/data/taxonomy_dictionaries/class_to_subclass.pkl')
subclass_to_division = read_pickle(bucket_name, '/data/taxonomy_dictionaries/subclass_to_division.pkl')
subclass_to_supertype = read_pickle(bucket_name, '/data/taxonomy_dictionaries/subclass_to_supertype.pkl')
supertype_to_subclass = read_pickle(bucket_name, '/data/taxonomy_dictionaries/supertype_to_subclass.pkl')
type_to_subclass = read_pickle(bucket_name, '/data/taxonomy_dictionaries/type_to_subclass.pkl')
    
# # Loading the dictionary from the file
# with open('/data/taxonomy_dictionaries/class_to_division.pkl', 'rb') as file:
#     class_to_division = pickle.load(file)
    
# with open('/data/taxonomy_dictionaries/division_to_class.pkl', 'rb') as file:
#     division_to_class = pickle.load(file)
    
# with open('/data/taxonomy_dictionaries/subclass_to_class.pkl', 'rb') as file:
#     subclass_to_class = pickle.load(file)
    
# with open('/data/taxonomy_dictionaries/class_to_subclass.pkl', 'rb') as file:
#     class_to_subclass = pickle.load(file)
    
# with open('/data/taxonomy_dictionaries/subclass_to_division.pkl', 'rb') as file:
#     subclass_to_division = pickle.load(file)
    
# with open('/data/taxonomy_dictionaries/subclass_to_supertype.pkl', 'rb') as file:
#     subclass_to_supertype = pickle.load(file)
    
# with open('/data/taxonomy_dictionaries/supertype_to_subclass.pkl', 'rb') as file:
#     supertype_to_subclass = pickle.load(file)
    
# with open('/data/taxonomy_dictionaries/sample_to_type.pkl', 'rb') as file:
#     sample_to_type = pickle.load(file) 
    
# with open('/data/taxonomy_dictionaries/type_to_subclass.pkl', 'rb') as file:
#     type_to_subclass = pickle.load(file)   
   

# Flatten lists where applicable
sample_to_subclass = {k: v if isinstance(v, str) else v[0] for k, v in sample_to_subclass.items() if v}

# Repeat for classes and divisions
sample_to_class = {
    sample: subclass_to_class.get(sub, None) for sample, sub in sample_to_subclass.items()
}

sample_to_division = {
    sample: class_to_division.get(cls, None) for sample, cls in sample_to_class.items()
}

class_to_supertype = {}

for cls, subclasses in class_to_subclass.items():
    supertypes = set()  # Use a set to avoid duplicate supertypes
    for sub in subclasses:
        if sub in subclass_to_supertype.keys():  # Check if subclass has an assigned supertype
            supertypes.update(subclass_to_supertype[sub])  # Add all supertypes from subclass
    
    class_to_supertype[cls] = list(supertypes)
    
################################################################################################################
    
# Function to convert row values to a binary string
def row_to_binary(row):
    return ''.join(map(str, row))

def binary_to_decimal(binary_str):
    return int(binary_str, 2)
                                            
################################################################################################################

# # Pick standard genes from expression_df and the average expression files

# # If you want to test only genes from a list use this
# gene_list_otherVendors = ['Acta2', 'Adgre1', 'Abcc8', 'Abcc9', 'Actl6b', 'Aif1', 'Akap5', 'Aldh5a1', 'App', 'Aqp1', 'Aqp4', 'Arg1', 'Bcl11b', 'Calca',
#             'Ccnd1', 'Cd247', 'Cd3e', 'Cd4', 'Cd5', 'Cd68', 'Cd86', 'Cd8a', 'Cdh1', 'Chat', 'Cnp', 'Cntnap1', 'Cntnap2', 'Col4a3/5/2/1',
#             'Creb1', 'Cspg4', 'Ctnnb1', 'Dbh', 'Dcx', 'Ddx5', 'Dlg2', 'Eea1', 'Eea5', 'Egr1', 'Emcn', 'Epm2a', 'Ewsr1', 'Fn1', 'Foxa2', 'Gad1', 'Gad2',
#             'Gad2/1', 'Gap43', 'Gfap', 'Gria2', 'Grin1', 'Grm2', 'Gsk3a', 'Gsk3a/b', 'Gucy1b1', 'Hcls1', 'Hopx', 'Htr2b', 'Htr7', 'Il10', 'Ins',
#             'Itgam', 'Itgax', 'Khdrbs1', 'Lamp1', 'Lyve1', 'Mag', 'Maoa', 'Maob', 'Map2', 'Mapk3', 'Mapk8/9/10', 'Mapt', 'Mbp', 'Mki67', 'Mog',
#             'Mrc1', 'Myb', 'Ncam1', 'Nefh', 'Nefl', 'Nefm', 'Nefm/h', 'Nfasc', 'Nfatc1', 'Nos1', 'Nos3', 'Npy', 'Nr3c2', 'Nrp1', 'Ntrk3', 'Ocrl', 'Oxt',
#             'P2rx4', 'P2ry12', 'Pax6', 'Pax7', 'Pdgfrb', 'Pecam1', 'Plp1', 'Ppp1r1b', 'Prkca/b/g', 'Pvalb', 'Pycard', 'Rbbp4', 'Rbfox3', 'S100a10',
#             'S100b', 'Satb2', 'Sdc4', 'Sdk2', 'Set', 'Sirt3', 'Slc1a2', 'Slc1a3', 'Slc6a3', 'Slc6a4', 'Snca', 'Sod2', 'Sox2', 'Sox4', 'Sox9', 'Sst',
#             'Stat1', 'Stx1a', 'Stx1a/1b/2/3', 'Sun2', 'Syn1', 'Syn2', 'Syp', 'Tardbp', 'Tbr1', 'Th', 'Tmem119', 'Tph1', 'Tph2', 'Tuba', 'Tubb', 'Tubb3',
#             'Uchl1', 'Vim']

# gene_list_neuromab = ['Adam11', 'Aldh1l1', 'Amigo1', 'Arx', 'Atp7a', 'Bdnf', 'Cacna1h', 'Cadm4', 'Calb1', 'Calb2', 'Clcn4', 'Cntnap1',
#                      'Dlg1', 'Dlg2', 'Dlg3', 'Dlg4', 'Drd2', 'Fgf13', 'Gabrb3', 'Gabre', 'Hspa9', 'Kcna1', 'Kcnd2', 'Lrp4',
#                      'Lrrk1', 'Lrrtm2', 'Lrrtm4', 'Mff', 'Mog', 'Nos1', 'Npy', 'Nrcam', 'Olfm1', 'Znf746', 'Pex5l', 'Qk',
#                      'Rbm17', 'Reep1', 'Reep2', 'Rufy3', 'S100a5', 'Shank1', 'Shank2', 'Shank3', 'Slc38a1', 'Snapin', 'Svop', 'Trpc4',
#                      'Vapa', 'Vdac1', 'Tpte']

# excluded_genes = ['Col4a3/5/2/1', 'Eea5', 'Gad2/1', 'Gsk3a/b', 'Ins', 
#                   'Mapk8/9/10', 'Nefm/h', 'Prkca/b/g', 'Stx1a/1b/2/3', 
#                   'Tuba', 'Tubb', 'Znf746', 'Qki']

# gene_list = [
#     gene for gene in (gene_list_otherVendors + gene_list_neuromab) 
#     if gene not in excluded_genes
# ]

# avg_expression_div_df = avg_expression_genesAll_div_df[gene_list]
# avg_expression_class_df = avg_expression_genesAll_class_df[gene_list]
# avg_expression_subclass_df = avg_expression_genesAll_subclass_df[gene_list]
# avg_expression_supertype_df = avg_expression_genesAll_supertype_df[gene_list]

# # gene_names = ['Arx', 'Vdac1', 'Reep1', 'Reep2', 'Actl6b', 'Abcc8', 'Abcc9', 'Clcn4','Aif1',
# #              'Rbm17', 'Epm2a', 'Ocrl', 'Cd3e', 'Sox9', 'Sun2', 'Aldh5a1', 'Sox4',
# #              'Tbr1', 'Tmem119', 'Tardbp', 'Ddx5', 'Rbbp4', 'Khdrbs1', 'Set', 'Dlg4', 'Gsk3a',
# #              'Pecam1', 'Eea1', 'Lamp1', 'Cd68', 'Bdnf', 'Rbfox3', 'Sod2', 'Sun2','Calb1', 'Calb2','Pvalb', 'Qk', 'Gfap', 'Nos1']

################################################################################################################

#genes_to_test = ['Rbfox3','Aif1','Nos1', 'Aldh5a1', 'Rbfox3', 'Arx', 'Calb1']

################################################################################################################

divisions = avg_expression_genesAll_div_df.index.values

neuronal_divs = divisions[0:4]
nonNeuronal_divs = divisions[4:]

# classes = avg_expression_genesAll_class_df.index.values
# subclasses = avg_expression_genesAll_subclass_df.index.values

neuronal_classes = [cls for d in neuronal_divs for cls in division_to_class.get(d, [])]
nonNeuronal_classes = [cls for d in nonNeuronal_divs for cls in division_to_class.get(d, [])]

division_colors = px.colors.qualitative.Set1 

subclass_sample_counts = {subclass: sum(1 for s in sample_to_subclass.values() if s == subclass) for subclass in avg_expression_subclass_df.index}
class_sample_counts = {cls: sum(1 for s in sample_to_class.values() if s == cls) for cls in avg_expression_class_df.index}
division_sample_counts = {div: sum(class_sample_counts.get(cls, 0) for cls in division_to_class.get(div, [])) for div in avg_expression_div_df.index}

neuronal_sample_count = sum(1 for div in sample_to_division.values() if div in neuronal_divs)
nonNeuronal_sample_count = sum(1 for div in sample_to_division.values() if div in nonNeuronal_divs)

################################################################################################################


################################################################################################################

app = dash.Dash(__name__)

# fig = go.Figure()
# fig_bin = go.Figure()
# fig_trin = go.Figure()

# app.layout = html.Div([
#     html.H1("Gene Expression", style={'textAlign': 'center'}),

#     html.Label("Enter genes (comma-separated):"),
#     dcc.Input(id='gene-input', type='text', value='Aif1', debounce=True, style={'width': '100%'}),

#     html.Button('Update Plot', id='update-button', n_clicks=0),
    
#     # Toggle for switching between binary and trinary plots
#         dcc.Checklist(
#             id='toggle-trinary',
#             options=[{'label': 'Show Trinarized Plot', 'value': 'trinary'}],
#             value=[],  # Default is empty (not checked)
#             style={'marginTop': '10px'}
#         ),
    
#     html.Div([
#         html.Div([
#             dcc.Graph(id='sunburst1'),
#             dcc.Graph(id='sunburst2')
#         ], style={'display': 'inline-block', 'width': '48%', 'verticalAlign': 'top'}),

#         html.Div(id='gene-bar-plots', style={'display': 'inline-block', 'width': '50%', 'verticalAlign': 'top', 'paddingLeft': '2%'})       
#     ], style={'display': 'flex', 'justify-content': 'space-between', 'marginTop': '20px'})
# ])

app.layout = html.Div([
    html.H1("Gene Expression", style={'textAlign': 'center'}),
    
    # Stores for precomputed figures (hidden storage)
    dcc.Store(id='store-fig-bin'),
    dcc.Store(id='store-fig-trin'),

    html.Label("Enter genes (comma-separated):"),
    dcc.Input(id='gene-input', type='text', value='Aif1, Gfap', debounce=True, style={'width': '100%'}),

    html.Div([
    html.Button('Update Plot', id='update-button', n_clicks=0)
    ], id="button-container"),
    
    html.Div([  
        # First Sunburst Plot
        html.Div([
            dcc.Graph(id='sunburst1', style={'height': '600px'})  # Fixed height for consistency
        ], style={'display': 'inline-block', 'width': '48%', 'position': 'relative', 'verticalAlign': 'top'}),

        # Second Sunburst Plot with Toggle Overlayed
        html.Div([
            # Toggle Positioned Over the Graph (Absolute Positioning)
            html.Div([
                dcc.Checklist(
                    id='toggle-trinary',
                    options=[{'label': ' Show Trinarized Plot', 'value': 'trinary'}],
                    value=[],
                    style={'textAlign': 'left'}
                )
            ], style={'position': 'absolute', 'top': '0px', 'left': '0px', 'zIndex': '10'}),  

            dcc.Graph(id='sunburst2', style={'height': '600px'})  # Match first sunburst
        ], style={'display': 'inline-block', 'width': '48%', 'position': 'relative', 'verticalAlign': 'top'})

    ], style={'display': 'flex', 'justify-content': 'space-between', 'alignItems': 'center', 'marginTop': '20px'}),

#     html.Div([
#         html.Div([
#             dcc.Graph(id='sunburst1')
#         ], style={'display': 'inline-block', 'width': '48%', 'verticalAlign': 'top'}),

#         html.Div([
#             dcc.Checklist(
#                 id='toggle-trinary',
#                 options=[{'label': ' Show Trinarized Plot', 'value': 'trinary'}],
#                 value=[],  # Default is unchecked
#                 style={'marginBottom': '10px'}
#             ),
#             style={'position': 'absolute', 'top': '0px', 'right': '10px'}
#             dcc.Graph(id='sunburst2')  # The figure affected by the toggle
#         ], style={'display': 'inline-block', 'width': '48%', 'verticalAlign': 'top'}),

#     ], style={'display': 'flex', 'justify-content': 'space-between', 'marginTop': '20px'}),

     html.Div(id='gene-bar-plots', style={'marginTop': '20px'})
])


@app.callback(
    [Output('sunburst1', 'figure'),
     Output('sunburst2', 'figure'),
     Output('gene-bar-plots', 'children'),
     Output('store-fig-bin', 'data'),
     Output('store-fig-trin', 'data')],
    [Input('update-button', 'n_clicks'),
     Input('toggle-trinary', 'value')], 
    [State('gene-input', 'value'),
     State('store-fig-bin', 'data'),
     State('store-fig-trin', 'data')]
)

def update_sunburst(n_clicks, toggle_value, gene_input, stored_fig_bin, stored_fig_trin):
    ctx = dash.callback_context  # Identifies what triggered the callback
    bar_plots = []

    if not ctx.triggered:  # If nothing has triggered the function yet
        return go.Figure(), go.Figure(), [], go.Figure().to_dict(), go.Figure().to_dict()

    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]  # Get the trigger ID
    
    #avg_expression_genesAll_subclass_df = pd.read_csv('/bil/users/psimko/holis/transcriptomic_analysis/avg_expression_subclass_genesAll_notNormalized_df.csv', index_col=0)
    
    # global avg_expression_genesAll_subclass_df
    # if avg_expression_genesAll_subclass_df is None:
    #     return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

    # If button is clicked, update sunburst and store new figures
    if trigger_id == "update-button":
        genes_to_test = [gene.strip() for gene in gene_input.split(',') if gene.strip()]
        print(f"Updating plots for genes: {genes_to_test}")
        print(genes_to_test)

        avg_expression_div_df = avg_expression_genesAll_div_df[genes_to_test]
        avg_expression_class_df = avg_expression_genesAll_class_df[genes_to_test]
        avg_expression_subclass_df = avg_expression_genesAll_subclass_df[genes_to_test]
        avg_expression_supertype_df = avg_expression_genesAll_supertype_df[genes_to_test]

        data_neuronal = []
        data_nonNeuronal = []

        # Add Subclasses (Linked to Classes)
        #for class_ in avg_expression_class_df.index:
        for subclass in avg_expression_subclass_df.index:
            parent_division = subclass_to_division.get(subclass, None) #
            parent_class = subclass_to_class.get(subclass, None)
            total_supertypes = len(subclass_to_supertype.get(subclass, []))
            if parent_class in neuronal_classes:
                nonzero_supertypes = sum(
                    (avg_expression_supertype_df.loc[s] >= 2).any()
                    for s in subclass_to_supertype.get(subclass, []) if s in avg_expression_supertype_df.index
                )
                supertype_proportion = nonzero_supertypes / total_supertypes if total_supertypes > 0 else 0

                adjusted_color = supertype_proportion * 100

                data_neuronal.append({"parent_division": parent_division, "parent_class": parent_class, "id": subclass, "value": subclass_sample_counts.get(subclass, 1), "color": adjusted_color})

        # Convert to DataFrame
        df_neuronal = pd.DataFrame(data_neuronal)

        # Add Subclasses (Linked to Classes)
        #for class_ in avg_expression_class_df.index:
        for subclass in avg_expression_subclass_df.index:
            parent_division = subclass_to_division.get(subclass, None) #
            parent_class = subclass_to_class.get(subclass, None)
            total_supertypes = len(subclass_to_supertype.get(subclass, []))
            if parent_class in nonNeuronal_classes:
                nonzero_supertypes = sum(
                    (avg_expression_supertype_df.loc[s] >= 2).any()
                    for s in subclass_to_supertype.get(subclass, []) if s in avg_expression_supertype_df.index
                )
                supertype_proportion = nonzero_supertypes / total_supertypes if total_supertypes > 0 else 0

                adjusted_color = supertype_proportion * 100

                data_nonNeuronal.append({"parent_division": parent_division, "parent_class": parent_class, "id": subclass, "value": subclass_sample_counts.get(subclass, 1), "color": adjusted_color})


        # Convert to DataFrame
        df_nonNeuronal = pd.DataFrame(data_nonNeuronal)

        # Create subplot layout with 1 row and 2 columns
        fig = make_subplots(rows=1, cols=2, 
                            specs=[[{"type": "domain"}, {"type": "domain"}]],  # "domain" for sunburst
                            subplot_titles=[f"Neuronal Cells {neuronal_sample_count:,}", f"Non-Neuronal Cells {nonNeuronal_sample_count:,}"])  # Add titles

        # Create first sunburst (Neuronal)
        fig1 = px.sunburst(
            df_neuronal,
            #path=["parent", "id"],
            path=["parent_division", "parent_class", "id"],
            values="value",
            color="color",
        )

        # Create second sunburst (Non-Neuronal)
        fig2 = px.sunburst(
            df_nonNeuronal,
            path=["parent_division", "parent_class", "id"],
            values="value",
            color="color",
        )

        # Add both sunburst plots to the figure
        fig.add_trace(fig1.data[0], row=1, col=1)  # Add first chart
        fig.add_trace(fig2.data[0], row=1, col=2)  # Add second chart

        # Adjust layout
        fig.update_layout(
            width=1200,  # Increase width (default is ~700)
            height=800,  # Increase height (default is ~450)
            title_text=f"{genes_to_test}",
            title_x=0.5,
            coloraxis=dict(
                    colorscale='blues',  # Set colorscale
                    cmin=0,  # Ensure minimum value is 0
                    cmax=100  # Ensure maximum value is 1
            ),
            showlegend=True,
        )

        fig.update_traces(marker=dict(line=dict(width=0.1, color='black')))

        ################################################################################################################
        ################################################################################################################

        # fig_to_use = fig_trin if 'trinary' in toggle_value else fig_bin

        expression_cells_bin_df = (avg_expression_subclass_df >= 2).astype(int)

        # Add one column with the cluster number and one with the corresponding binary expression
        if 'binary_signature' in expression_cells_bin_df.columns:
            expression_cells_bin_df.drop('binary_signature', axis=1, inplace=True)
        if 'Cluster' in expression_cells_bin_df.columns:
            expression_cells_bin_df.drop('Cluster', axis=1, inplace=True)

        expression_cells_bin_df['binary_signature'] = expression_cells_bin_df.apply(row_to_binary, axis=1) 
        expression_cells_bin_df['Cluster'] = expression_cells_bin_df['binary_signature'].apply(binary_to_decimal)

        # Use these this line if you want to use binary expression of subclasses, otherwise the one above
        expression_cells_bin_df["subclass"] = expression_cells_bin_df.index

        # Create the dictionary using groupby
        subclass_to_binary_signature = expression_cells_bin_df.set_index('subclass')['binary_signature'].to_dict()

        # Add binary signature column with each entry reflecting the code of its corresponding subclass
        df_neuronal['binary_signature'] = df_neuronal['id'].map(subclass_to_binary_signature)
        df_nonNeuronal["binary_signature"] = df_nonNeuronal['id'].map(subclass_to_binary_signature)


        #Create a colormap

        # Define a fixed set of distinct colors 
        distinct_colors = plotly.colors.qualitative.Safe + plotly.colors.qualitative.Dark24

        # Ensure we have enough colors (if not, repeat or extend)
        num_unique_signatures = expression_cells_bin_df['binary_signature'].nunique()
        if num_unique_signatures > len(distinct_colors):
            raise ValueError(f"Not enough distinct colors for {num_unique_signatures} unique binary signatures.")

        # Get unique binary signatures
        unique_signatures = expression_cells_bin_df['binary_signature'].unique()

        # Create a mapping dictionary for binary_signature → color
        binary_signature_to_color = {sig: distinct_colors[i] for i, sig in enumerate(unique_signatures)}

        color_map = {
        **binary_signature_to_color.copy(),
        **{d: "white" for d in df_neuronal["parent_division"].dropna().unique()},  # Ensure unique & no NaNs
        **{c: "white" for c in df_neuronal["parent_class"].dropna().unique()},     # Ensure unique & no NaNs
        **{'CB Glut': "white"},
        '(?)': 'white'  # Default unknown category color
        }


        #Plot

        # Create subplot layout with 1 row and 2 columns
        fig_bin = make_subplots(
            rows=1, cols=2, 
            specs=[[{"type": "domain"}, {"type": "domain"}]],  # "domain" for sunburst
            subplot_titles=[f"Neuronal Cells {neuronal_sample_count:,}", f"Non-Neuronal Cells {nonNeuronal_sample_count:,}"]
        )

        # Create first sunburst (Neuronal)
        fig1 = px.sunburst(
            df_neuronal,
            path=["parent_division", "parent_class", "id"],  # Keep hierarchy
            values="value",
            color="binary_signature",  # Color by binary signature
            color_discrete_map=color_map  # Use predefined colors
        )

        # Create second sunburst (Non-Neuronal)
        fig2 = px.sunburst(
            df_nonNeuronal,
            path=["parent_division", "parent_class", "id"],  # Keep hierarchy
            values="value",
            color="binary_signature",  # Color by binary signature
            color_discrete_map=color_map  # Use predefined colors
        )

        # Add both sunburst plots to the figure
        fig_bin.add_trace(fig1.data[0], row=1, col=1)  # Add first chart
        fig_bin.add_trace(fig2.data[0], row=1, col=2)  # Add second chart

        # **Manually Add a Custom Legend for Binary Signatures**
        legend_traces = []
        for binary_signature, color in binary_signature_to_color.items():
            legend_traces.append(go.Scatter(
                x=[None], y=[None],  # No actual data points, just legend markers
                mode='markers',
                marker=dict(size=12, color=color),
                name=str(binary_signature),  # Label the legend with binary signature
                showlegend=True
            ))

        # Adjust layout
        fig_bin.update_layout(
            paper_bgcolor='white',  # Removes the grid-like background
            plot_bgcolor='white',   # Ensures no grid in the plot area
            width=1200,  # Increase width
            height=800,  # Increase height
            xaxis=dict(showticklabels=False, showgrid=False, zeroline=False, title=""),
            yaxis=dict(showticklabels=False, showgrid=False, zeroline=False, title=""),
            title_text=f"{genes_to_test}",
            title_x=0.5,
            showlegend=True,  # Ensure legend is displayed
            #margin=dict(t=50, b=10, l=10, r=10),
            legend=dict(title="Binary Clusters", font=dict(size=12))
        )

        # Add legend traces to figure
        for trace in legend_traces:
            fig_bin.add_trace(trace)

        fig_bin.update_traces(marker=dict(line=dict(width=0.1, color='black')))

        ################################################################################################################
        ################################################################################################################


        # avg_expression_genesAll_subclass_df = avg_expression_genesAll_subclass_df[genes_to_test]

        # Apply trinarization
        # expression_cells_trin_df = avg_expression_subclass_df.copy()

        expression_cells_trin_df = avg_expression_subclass_df.applymap(
            lambda x: 0 if x < 2 else (1 if x < 7 else 2)
        )

        # Add one column with the cluster number and one with the corresponding trinary expression

        if 'trinary_signature' in expression_cells_trin_df.columns:
            expression_cells_trin_df.drop('trinary_signature', axis=1, inplace=True)

        expression_cells_trin_df['trinary_signature'] = expression_cells_trin_df.apply(row_to_binary, axis=1) 
        expression_cells_trin_df["subclass"] = expression_cells_trin_df.index

        # Create the dictionary using groupby
        trinary_signature_to_subclass = expression_cells_trin_df.groupby('trinary_signature')['subclass'].unique().to_dict()

        trinary_signature_to_subclass = {k: list(v) for k, v in trinary_signature_to_subclass.items()}

        subclass_to_trinary_signature = expression_cells_trin_df.set_index('subclass')['trinary_signature'].to_dict()

        df_neuronal['trinary_signature'] = df_neuronal['id'].map(subclass_to_trinary_signature)
        df_nonNeuronal["trinary_signature"] = df_nonNeuronal['id'].map(subclass_to_trinary_signature)


        # Ensure we have enough colors (if not, repeat or extend)
        num_unique_signatures = expression_cells_trin_df['trinary_signature'].nunique()
        if num_unique_signatures > len(distinct_colors):
            raise ValueError(f"Not enough distinct colors for {num_unique_signatures} unique trinary signatures.")

        # Get unique trinary signatures
        unique_signatures = expression_cells_trin_df['trinary_signature'].unique()

        # Create a mapping dictionary for trinary_signature → color
        trinary_signature_to_color = {sig: distinct_colors[i] for i, sig in enumerate(unique_signatures)}

        # Ensure that divisions and classes are white, only subclasses are colored
        color_map = trinary_signature_to_color.copy()

        # Modify color map to include transparency (RGBA format)
        color_map = {
            **trinary_signature_to_color.copy(),
            '(?)': 'white'  # Default unknown category color
        }


        # Create subplot layout with 1 row and 2 columns
        fig_trin = make_subplots(
            rows=1, cols=2, 
            specs=[[{"type": "domain"}, {"type": "domain"}]],  # "domain" for sunburst
            subplot_titles=[f"Neuronal Cells {neuronal_sample_count:,}", f"Non-Neuronal Cells {nonNeuronal_sample_count:,}"]
        )

        # Create first sunburst (Neuronal)
        fig1 = px.sunburst(
            df_neuronal,
            path=["parent_division", "parent_class", "id"],  # Keep hierarchy
            values="value",
            color="trinary_signature",  # Color by binary signature
            color_discrete_map=color_map  # Use predefined colors
        )

        # Create second sunburst (Non-Neuronal)
        fig2 = px.sunburst(
            df_nonNeuronal,
            path=["parent_division", "parent_class", "id"],  # Keep hierarchy
            values="value",
            color="trinary_signature",  # Color by binary signature
            color_discrete_map=color_map  # Use predefined colors
        )

        # Add both sunburst plots to the figure
        fig_trin.add_trace(fig1.data[0], row=1, col=1)  # Add first chart
        fig_trin.add_trace(fig2.data[0], row=1, col=2)  # Add second chart

        # **Manually Add a Custom Legend for Binary Signatures**
        legend_traces = []
        for trinary_signature, color in trinary_signature_to_color.items():
            legend_traces.append(go.Scatter(
                x=[None], y=[None],  # No actual data points, just legend markers
                mode='markers',
                marker=dict(size=12, color=color),
                name=str(trinary_signature),  # Label the legend with binary signature
                showlegend=True
            ))

        # Adjust layout
        fig_trin.update_layout(
            paper_bgcolor='white',  # Removes the grid-like background
            plot_bgcolor='white',   # Ensures no grid in the plot area
            width=1200,  # Increase width
            height=800,  # Increase height
            xaxis=dict(showticklabels=False, showgrid=False, zeroline=False, title=""),
            yaxis=dict(showticklabels=False, showgrid=False, zeroline=False, title=""),
            title_text=f"{genes_to_test}",
            title_x=0.5,
            showlegend=True,  # Ensure legend is displayed
            #margin=dict(t=50, b=10, l=10, r=10),
            legend=dict(title="Trinary Clusters", font=dict(size=12))
        )

        # Add legend traces to figure
        for trace in legend_traces:
            fig_trin.add_trace(trace)

        fig_trin.update_traces(marker=dict(line=dict(width=0.1, color='black')))

        ################################################################################################################
        ################################################################################################################

        unique_divisions = avg_expression_genesAll_div_df.index
        unique_classes = avg_expression_genesAll_class_df.index
        print(unique_divisions)

        # Create a color map for divisions
        division_colors = plt.cm.jet(np.linspace(0, 1, len(unique_divisions)))
        division_colors_hex = [mcolors.to_hex(color) for color in division_colors]
        division_colors_dict = {division: division_colors_hex[i] for i, division in enumerate(unique_divisions)}

        bar_plots = []
        for gene in genes_to_test:
            class_expression = avg_expression_class_df[gene]
            sorted_classes = []
            for division in unique_divisions:
                division_classes = [cls for cls in class_expression.index if class_to_division[cls] == division]
                sorted_classes.extend(division_classes)

            class_expression = class_expression.loc[sorted_classes]

            # class_colors = [division_colors_dict[class_to_division[cls]] for cls in sorted_classes]
            # class_colors_hex = [mcolors.to_hex(color) for color in class_colors]
            class_colors_dict = {cls: division_colors_dict[class_to_division[cls]] for cls in class_to_division}


            class_expression['colors'] = class_expression.index.map(class_colors_dict)

            expression_values = class_expression[sorted_classes].values


            # Calculate error bars (asymmetric, prevent negative)
            division_std_dev = avg_expression_class_df[gene].std()
            upper_error = np.full_like(expression_values, division_std_dev)
            lower_error = np.minimum(expression_values, division_std_dev)

            fig_bars = go.Figure()

            # Add bars
            fig_bars.add_trace(go.Bar(
                x=sorted_classes,
                y=expression_values,
                marker_color=class_expression['colors'],
                error_y=dict(
                    type='data',
                    symmetric=False,
                    array=upper_error,
                    arrayminus=lower_error,
                    thickness=1.5,
                    width=3,
                    color='rgba(0,0,0,0.5)'
                )
            ))


            for division in unique_divisions:
                fig_bars.add_trace(go.Scatter(
                    x=[sorted_classes[0]],  # Use a real x-value
                    y=[0],  # Set y to a dummy value like 0
                    mode='markers',
                    marker=dict(color=division_colors_dict[division], size=10),
                    name=division,
                    legendgroup=division,  # Group the bars with the legend
                    showlegend=True
                ))

            fig_bars.update_layout(
                title=f'Gene Expression in Classes: {gene}',
                xaxis_title='Classes',
                yaxis_title='Expression Intensity',
                barmode='group',
                legend=dict(itemsizing='constant',traceorder='normal'),
                legend_title='Divisions',
                yaxis=dict(range=[0, 10]),
                width=900,
                height=500,
                plot_bgcolor='white'
            )

            fig_bars.update_xaxes(tickangle=45, tickfont=dict(size=10))

            bar_plots.append(dcc.Graph(figure=fig_bars))


        return fig, fig_bin, bar_plots, fig_bin.to_dict(), fig_trin.to_dict()
    
    ########################################################################

    # If toggle is switched, retrieve stored figures and update `sunburst2`
    elif trigger_id == "toggle-trinary":
        if stored_fig_bin is None or stored_fig_trin is None:
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        fig_bin = go.Figure(stored_fig_bin)
        fig_trin = go.Figure(stored_fig_trin)
        fig_to_use = fig_trin if 'trinary' in toggle_value else fig_bin

        return dash.no_update, fig_to_use, dash.no_update, dash.no_update, dash.no_update

    ################################################################################################################
    ################################################################################################################

if __name__ == '__main__':
    app.run_server(host="0.0.0.0", port=8050, debug=False)