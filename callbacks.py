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
import matplotlib.pyplot as plt
from numpy import inf

################################################################################################################
##### Helper functions (can be moved to utils later)
################################################################################################################

# Function to convert row values to a binary string
def row_to_binary(row):
    return ''.join(map(str, row))

def binary_to_decimal(binary_str):
    return int(binary_str, 2)

################################################################################################################
##### Callback decorator
################################################################################################################

@callback(
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

################################################################################################################
##### Actual callback
################################################################################################################

                                   
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
        ##### Sunburst 2 - trinary
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
        ##### Bar plots - expression level in classes
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
