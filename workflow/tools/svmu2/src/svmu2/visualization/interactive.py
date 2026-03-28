'''
Docstring for visualization.interactive
Plotly interactive dotplot
'''

import plotly.graph_objects as go
import plotly.io as pio

def plot_interactive_dotplot(best_delta, new_SVs=None, x_marker=None):
    plotly_data = []

    # Add SVs as scatter lines
    if new_SVs:
        for sv in new_SVs:
            color = 'yellow' if sv.sv_type == 'BND' else 'orange'
            trace = go.Scatter(
                x=[sv.reference_start, sv.reference_end],
                y=[sv.query_start, sv.query_end],
                mode='lines',
                line=dict(color=color),
                hovertext=f"SV type: {sv.sv_type}",
                hoverinfo='text'
            )
            plotly_data.append(trace)

    for block in best_delta.alignment_blocks:
        # Always plot the original block
        block_color = (
            'grey' if (not block.is_repeat and not block.part_of_primary_synteny) else
            'blue' if block.part_of_primary_synteny else
            'black'
        )
        hover_text = (
            f"Index: {block.index}<br>"
            f"Ref: {block.reference_start}-{block.reference_end}<br>"
            f"Query: {block.query_start}-{block.query_end}<br>"
            f"Repeat: {block.is_repeat}<br>"
            f"Primary Synteny: {block.part_of_primary_synteny}<br>"
            f"Reflected: {getattr(block, 'reflected', False)}"
        )
        orig_trace = go.Scatter(
            x=[block.reference_start, block.reference_end],
            y=[block.query_start, block.query_end],
            mode='lines',
            line=dict(color=block_color),
            hovertext=hover_text,
            hoverinfo='text',
        )
        plotly_data.append(orig_trace)

        # If reflected, also plot the reflected coordinates
        if block.reflected:
            reflect_trace = go.Scatter(
                x=[block.reference_start, block.reference_end],
                y=[block.y1_reflection, block.y2_reflection],
                mode='lines',
                line=dict(color='green'),
                hovertext = (
                    f"Index: {block.index}<br>"
                    f"Ref: {block.reference_start}-{block.reference_end}<br>"
                    f"Query: {block.y1_reflection}-{block.y2_reflection}<br>"
                    f"Repeat: {block.is_repeat}<br>"
                    f"Primary Synteny: {block.part_of_primary_synteny}<br>"
                    f"Reflected: {getattr(block, 'reflected', False)}"
                ),
                hoverinfo='text',
            )
            print('working')
            plotly_data.append(reflect_trace)

    # Add vertical line at x_marker (e.g., x_target for regression line)
    if x_marker is not None:
        miny = min(b.query_start for b in best_delta.alignment_blocks)
        maxy = max(b.query_end for b in best_delta.alignment_blocks)
        line_trace = go.Scatter(
            x=[x_marker, x_marker],
            y=[miny, maxy],
            mode='lines',
            line=dict(color='red', dash='dash'),
            name='X Marker',
            hoverinfo='skip'
        )
        plotly_data.append(line_trace)

    layout = go.Layout(
        title=f"Interactive Dotplot: {best_delta.reference}",
        xaxis=dict(title='Reference'),
        yaxis=dict(title='Query'),
        showlegend=False
    )
    
    plotly_fig = go.Figure(data=plotly_data, layout=layout)
    pio.write_html(plotly_fig, file=f'{best_delta.reference}_interactive.html', auto_open=True)
