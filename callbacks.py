import dash
from dash import callback, dcc, html, Input, Output, State
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from app_context import * 
from plotly.subplots import make_subplots
import plotly.express as px
import plotly.colors
import matplotlib.colors as mcolors
from plotly.colors import sample_colorscale
import matplotlib.pyplot as plt
from numpy import inf
import colorcet as cc
import re
import utils
from utils import row_to_binary, binary_to_decimal, update_bar_plots, parse_gene_input, pool_genes_in_df, get_gene_label


################################################################################################################
##### Callback decorator
################################################################################################################

@callback(
    [Output('sunburst1', 'figure'),
     Output('sunburst2', 'figure'),
     Output('gene-bar-plots', 'children'),
     Output('store-fig-bin', 'data'),
     Output('store-fig-trin', 'data'),
     Output('store-genes', 'data'),
     Output('store-expression-div', 'data'),
     Output('store-expression-class', 'data')],
     # Output('store-expression-subclass', 'data'),
     # Output('store-expression-supertype', 'data')],
    [Input('update-button', 'n_clicks'),
     Input('toggle-trinary', 'value')], 
    [State('gene-input', 'value'),
     State('store-fig-bin', 'data'),
     State('store-fig-trin', 'data'),
     State('store-genes', 'data'),
     State('store-expression-div', 'data'),
     State('store-expression-class', 'data')]
)

################################################################################################################
##### Actual callback
################################################################################################################

                                   
def update_sunburst(n_clicks, toggle_value, gene_input, stored_fig_bin, stored_fig_trin, store_genes, stored_expression_div, stored_expression_class):
    if toggle_value is None:
        toggle_value = ['binary']
    ctx = dash.callback_context  # Identifies what triggered the callback
    bar_plots = []

    if not ctx.triggered:  # If nothing has triggered the function yet
        return go.Figure(), go.Figure(), [], go.Figure().to_dict(), go.Figure().to_dict(), []

    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]  # Get the trigger ID
    
    #avg_expression_genesAll_subclass_df = pd.read_csv('/bil/users/psimko/holis/transcriptomic_analysis/avg_expression_subclass_genesAll_notNormalized_df.csv', index_col=0)
    
    # global avg_expression_genesAll_subclass_df
    # if avg_expression_genesAll_subclass_df is None:
    #     return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

    # If button is clicked, update sunburst and store new figures
    if trigger_id == "update-button":
        # genes_to_test = [gene.strip() for gene in gene_input.split(',') if gene.strip()]
        # num_genes = len(genes_to_test)
        genes_to_test = parse_gene_input(gene_input)
        num_genes = sum(len(g) if isinstance(g, list) else 1 for g in genes_to_test)
        
        print(f"Updating plots for genes: {genes_to_test}")
        print(genes_to_test)

        # avg_expression_div_df = avg_expression_genesAll_div_df[genes_to_test]
        # avg_expression_class_df = avg_expression_genesAll_class_df[genes_to_test]
        # avg_expression_subclass_df = avg_expression_genesAll_subclass_df[genes_to_test]
        # avg_expression_supertype_df = avg_expression_genesAll_supertype_df[genes_to_test]

        avg_expression_div_df = pool_genes_in_df(avg_expression_genesAll_div_df, genes_to_test)
        avg_expression_class_df = pool_genes_in_df(avg_expression_genesAll_class_df, genes_to_test)
        avg_expression_subclass_df = pool_genes_in_df(avg_expression_genesAll_subclass_df, genes_to_test)
        avg_expression_supertype_df = pool_genes_in_df(avg_expression_genesAll_supertype_df, genes_to_test)
        
        ################################################################################################################
        ##### Sunburst 1 - expression (% of expressed supertypes in subclasses for given genes)
        ################################################################################################################

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
        ##### Sunburst 2 - binary
        ################################################################################################################

        # fig_to_use = fig_trin if 'trinary' in toggle_value else fig_bin

        # binary_threshold = 2
        # expression_cells_bin_df = (avg_expression_subclass_df >= binary_threshold).astype(int)

        # Binarize expressions
        expression_cells_bin_df = avg_expression_subclass_df.copy()
        num_columns = expression_cells_bin_df.shape[1]
        print(f'expression_cells_bin_df= {expression_cells_bin_df.columns}')
        for gene in expression_cells_bin_df.columns:
            label = get_gene_label(gene)
            threshold = gene_thresholds.get(label, [2])[0]  # first value, fallback to 2
            expression_cells_bin_df[label] = (expression_cells_bin_df[label] >= threshold).astype(int)

        print(expression_cells_bin_df)

        # Add one column with the cluster number and one with the corresponding binary expression
        if 'binary_signature' in expression_cells_bin_df.columns:
            expression_cells_bin_df.drop('binary_signature', axis=1, inplace=True)
        if 'Cluster' in expression_cells_bin_df.columns:
            expression_cells_bin_df.drop('Cluster', axis=1, inplace=True)

        expression_cells_bin_df['binary_signature'] = expression_cells_bin_df.apply(row_to_binary, axis=1) 
        expression_cells_bin_df['Cluster'] = expression_cells_bin_df['binary_signature'].apply(binary_to_decimal)

        # Use these this line if you want to use binary expression of subclasses, otherwise the one above
        expression_cells_bin_df["subclass"] = expression_cells_bin_df.index

        print(expression_cells_bin_df)

        # Create the dictionary using groupby
        subclass_to_binary_signature = expression_cells_bin_df.set_index('subclass')['binary_signature'].to_dict()

        # Add binary signature column with each entry reflecting the code of its corresponding subclass
        df_neuronal['binary_signature'] = df_neuronal['id'].map(subclass_to_binary_signature)
        df_nonNeuronal["binary_signature"] = df_nonNeuronal['id'].map(subclass_to_binary_signature)


        #Create a colormap

        # Define a fixed set of distinct colors 
        #distinct_colors = plotly.colors.qualitative.Safe + plotly.colors.qualitative.Dark24
        distinct_colors = distinct_colors = cc.glasbey[:256]

        # Ensure we have enough colors (if not, repeat or extend)
        num_unique_signatures = expression_cells_bin_df['binary_signature'].nunique()
        if num_unique_signatures > len(distinct_colors):
            raise ValueError(f"Not enough distinct colors for {num_unique_signatures} unique binary signatures.")

        # Get unique binary signatures
        unique_signatures = expression_cells_bin_df['binary_signature'].unique()
        possible_binary_signatures = 2**num_columns
        
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
            legend=dict(title=f"Binary Clusters<br>{num_unique_signatures}/{possible_binary_signatures}", font=dict(size=12))
        )

        # Add legend traces to figure
        for trace in legend_traces:
            fig_bin.add_trace(trace)

        fig_bin.update_traces(marker=dict(line=dict(width=0.1, color='black')))

        ################################################################################################################
        ##### Sunburst 2 - trinary
        ################################################################################################################

        # Apply trinarization
        expression_cells_trin_df = avg_expression_subclass_df.copy()
   
        # trinary_threshold_low = 2
        # trinary_threshold_high = 7
        # expression_cells_trin_df = avg_expression_subclass_df.applymap(
        #     lambda x: 0 if x < trinary_threshold_low else (1 if x < trinary_threshold_high else 2)
        # )

        def trinarize_column(col, gene_thresholds, default_low=2, default_high=7):
            gene = col.name
            low, high = gene_thresholds.get(gene, [default_low, default_high])
            return col.apply(lambda x: 0 if x < low else (1 if x < high else 2))
        
        expression_cells_trin_df = avg_expression_subclass_df.apply(
            lambda col: trinarize_column(col, gene_thresholds),
            axis=0
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
        possible_trinary_signatures = 3**num_columns

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
            legend=dict(title=f"Trinary Clusters<br>{num_unique_signatures}/{possible_trinary_signatures}", font=dict(size=12))
        )

        # Add legend traces to figure
        for trace in legend_traces:
            fig_trin.add_trace(trace)

        fig_trin.update_traces(marker=dict(line=dict(width=0.1, color='black')))

        ################################################################################################################
        ##### Bar plots - expression level in classes
        ################################################################################################################

        unique_divisions = avg_expression_genesAll_div_df.index
        unique_classes = avg_expression_genesAll_class_df.index

        # Create a color map for divisions
        division_colors = plt.cm.jet(np.linspace(0, 1, len(unique_divisions)))
        division_colors_hex = [mcolors.to_hex(color) for color in division_colors]
        division_colors_dict = {division: division_colors_hex[i] for i, division in enumerate(unique_divisions)}

        return (
            fig,
            fig_bin,
            update_bar_plots(
                toggle_value,
                genes_to_test,
                avg_expression_div_df,
                avg_expression_class_df,
                gene_thresholds,
                binary_threshold=2,
                trinary_threshold_low=2,
                trinary_threshold_high=7
            ), 
            fig_bin.to_dict(),
            fig_trin.to_dict(),
            genes_to_test,
            avg_expression_div_df.to_dict(),     
            avg_expression_class_df.to_dict(),   
        )
    
    ########################################################################

    # If toggle is switched, retrieve stored figures and update `sunburst2`
    elif trigger_id == "toggle-trinary":
        if stored_fig_bin is None or stored_fig_trin is None:
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        fig_bin = go.Figure(stored_fig_bin)
        fig_trin = go.Figure(stored_fig_trin)
        fig_to_use = fig_trin if 'trinary' in toggle_value else fig_bin
        genes_to_test = store_genes

        avg_expression_div_df = pd.DataFrame(stored_expression_div)
        avg_expression_class_df = pd.DataFrame(stored_expression_class)

        return (dash.no_update, 
                fig_to_use, 
                update_bar_plots(
                    toggle_value, 
                    genes_to_test, 
                    avg_expression_div_df, 
                    avg_expression_class_df, 
                    gene_thresholds, 
                    binary_threshold=2, 
                    trinary_threshold_low=2, 
                    trinary_threshold_high=7
                ), 
                dash.no_update,
                dash.no_update,
                dash.no_update,
                dash.no_update,
                dash.no_update
        )

    ################################################################################################################
    ################################################################################################################
