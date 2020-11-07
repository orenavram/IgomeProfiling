"""
Union distinctive and identified motifs to single table
"""
from os import path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

colors_map = {
    'perfect': 'green',
    'mixed': 'blue',
    'incorrect': 'red',
    'unexpected': 'orange',
    'artifact': 'red',
    'negative': 'orange',
    'irrelevant': 'grey'
}

def render_mpl_table(data, col_width=3.0, row_height=0.625, font_size=14,
                     header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                     bbox=[0, 0, 1, 1], header_columns=0, color_from_index=2,
                     output_path=None, ax=None, **kwargs):
    if ax is None:
        size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
        fig, ax = plt.subplots(figsize=size)
        ax.axis('off')
    mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, cellLoc="left", **kwargs)
    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(font_size)

    for k, cell in mpl_table._cells.items():
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        elif k[1] >= color_from_index:
            cell.set_facecolor(colors_map[data.iat[k[0] - 1, k[1]]])
        else:
            cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
    if output_path:
        fig.savefig(output_path)
    return ax.get_figure(), ax

# TODO refactor
def process(results, top_motifs, output_base_path: str):
    motifs = {}
    # Index
    for result in results:
        source = result[0]
        df = pd.read_csv(result[1])
        df = df[['motif', 'label']]
        for _, row in df.iterrows():
            motif = row['motif']
            label = row['label']          
            motif_data = None
            try:
                motif_data = motifs[motif]
            except:
                motif_data = { 'is_top': motif in top_motifs }
                motifs[motif] = motif_data
            motif_data[source] = label
    # TODO extract details?
    # Flatten
    motif_keys = list(motifs.keys())
    first_key = motif_keys[0]
    columns = list(motifs[first_key].keys())
    data = []
    for motif in motif_keys:
        motif_data = motifs[motif]
        motif_flat_data = [motif]
        for column in columns: motif_flat_data.append(motif_data[column])
        data.append(motif_flat_data)
    # Save to csv
    print('Saving results...')
    output_path = f'{output_base_path}.csv'
    df = pd.DataFrame(data, columns=['motif'] + columns)
    df.to_csv(output_path, index=False)
    # Save table (svg)
    # TODO group by is_top
    output_path = f'{output_base_path}_all.svg'
    render_mpl_table(df, output_path=output_path)
    df_groupped = df.groupby('is_top')
    for is_top, group in df_groupped:
        name = 'top' if is_top else 'others'
        output_path = f'{output_base_path}_{name}.svg'
        group = group.drop('is_top', axis=1)
        render_mpl_table(group, output_path=output_path, color_from_index=1)


if __name__ == '__main__':
    base_path = '/mnt/d/workspace/data/webiks/igome/results/exp4+9_all_half_gaps/model_fitting'
    top_motifs = ['CSTTTTTRTC', 'CDTDGWIRMDPC', 'CPSHISLRNPC', 'CTRPTRPSLWMSPAC', 'CTNLCGTLAC', 'CCHHIFQRPPC', 'CPACSLFTC', 'CDAHSLFFPC', 'CGGTAGSFSRC', 'CLSSLFTRC']
    results = [
        ['base', path.join(base_path, 'ferret_texas/distinctive_motifs_all.csv')],
        ['reduced', path.join(base_path,'ferret_texas/exp4_texas_motifs_on_exp8_reduced.csv')],
        ['complete', path.join(base_path,'ferret_texas/exp4_texas_motifs_on_exp8_complete.csv')],
    ]
    # TODO add union
    output_base_path = path.join(base_path, 'ferret_texas/cross_summary')
    process(results, top_motifs, output_base_path)