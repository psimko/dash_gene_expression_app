import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import plotly.graph_objects as go
from dash import dcc
from app_context import * 

################################################################################################################
##### Helper functions 
################################################################################################################

# Function to convert row values to a binary string
def row_to_binary(row):
    return ''.join(map(str, row))

def binary_to_decimal(binary_str):
    return int(binary_str, 2)

def update_bar_plots(toggle_value, 
                     genes_to_test,
                     avg_expression_div_df, 
                     avg_expression_class_df,
                     gene_thresholds,
                     binary_threshold=2, 
                     trinary_threshold_low=2, 
                     trinary_threshold_high=7
                    ):
    
    if not toggle_value:
        toggle_value = ['binary']

    # unique_divisions = avg_expression_genesAll_div_df.index
    # unique_classes = avg_expression_genesAll_class_df.index
    unique_divisions = avg_expression_div_df.index
    unique_classes = avg_expression_class_df.index

    # Create a color map for divisions
    division_colors = plt.cm.jet(np.linspace(0, 1, len(unique_divisions)))
    division_colors_hex = [mcolors.to_hex(color) for color in division_colors]
    division_colors_dict = {division: division_colors_hex[i] for i, division in enumerate(unique_divisions)}
 
    bar_plots = []
    for gene in genes_to_test:
        gene_label = get_gene_label(gene)
        class_expression = avg_expression_class_df[gene_label]
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
        division_std_dev = avg_expression_class_df[gene_label].std()
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
    
        ########################################################################
        # Determine which thresholds to show
        ########################################################################

        # Get gene-specific thresholds, or use global defaults
        thresholds = gene_thresholds.get(gene_label, [binary_threshold, trinary_threshold_high])
    
        if gene_label in gene_thresholds:
            binary_thr = gene_thresholds[gene_label][0]
            trinary_thr_low = gene_thresholds[gene_label][0]
            trinary_thr_high = gene_thresholds[gene_label][1]
        else:
            binary_thr = binary_threshold
            trinary_thr_low = trinary_threshold_low
            trinary_thr_high = trinary_threshold_high
            
        threshold_shapes = []
    
        if 'binary' in toggle_value:
            threshold_shapes.append(dict(
                type='line',
                xref='paper', yref='y',
                x0=0, x1=1,
                y0=binary_thr, y1=binary_thr,
                line=dict(color='red', width=2, dash='dash')
            ))
            
        elif 'trinary' in toggle_value:
            threshold_shapes.extend([
                dict(
                    type='line',
                    xref='paper', yref='y',
                    x0=0, x1=1,
                    y0=trinary_thr_low, y1=trinary_thr_low,
                    line=dict(color='orange', width=2, dash='dash')
                ),
                dict(
                    type='line',
                    xref='paper', yref='y',
                    x0=0, x1=1,
                    y0=trinary_thr_high, y1=trinary_thr_high,
                    line=dict(color='green', width=2, dash='dash')
                )
            ])
    
        #######################################################################
    
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
            plot_bgcolor='white',
            shapes=threshold_shapes
        )
    
        print("Threshold shapes:", threshold_shapes)
    
        fig_bars.update_xaxes(tickangle=45, tickfont=dict(size=10))
    
        bar_plots.append(dcc.Graph(figure=fig_bars))
    
    return bar_plots

# Parse gene input with (possibly) pooled genes
def parse_gene_input(gene_input):
    tokens = re.findall(r'\[[^\]]+\]|[^,\[\]\s]+', gene_input)
    parsed = []
    for token in tokens:
        token = token.strip()
        if token.startswith('[') and token.endswith(']'):
            # Remove brackets and split inside
            pooled_genes = [g.strip() for g in token[1:-1].split(',') if g.strip()]
            parsed.append(pooled_genes)
        else:
            parsed.append(token)
    return parsed

# Pool genes in a dataframe
def pool_genes_in_df(df, genes_to_test):
    
    """Here genes_to_test = ["Rbfox3", ["Gad1", "Acta2"], "Pvalb", ["Gfap", "Aif1"]] for example"""
    
    df_copy = df.copy()
    pooled_columns = []

    for entry in genes_to_test:
        if isinstance(entry, list):
            # Make a column name like "Gad1+Acta2"
            col_name = '+'.join(entry)
            # Sum the selected genes across rows
            df_copy[col_name] = df_copy[entry].sum(axis=1)
            pooled_columns.append(col_name)
        else:
            pooled_columns.append(entry)
    
    return df_copy[pooled_columns]

def get_gene_label(g):
    return '+'.join(g) if isinstance(g, list) else g