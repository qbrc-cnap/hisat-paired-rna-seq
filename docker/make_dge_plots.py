import argparse
import sys
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from bokeh.plotting import figure, output_file, save
from bokeh.models import ColumnDataSource, HoverTool

sns.set_style('darkgrid')

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', 
        dest='dg_table', 
        required=True, 
        help="Path to a DESeq2-format output table (tsv)"
    )

    parser.add_argument('-c',
        dest='counts_table',
        required=True,
        help="Path to a normalized count matrix"
    )

    parser.add_argument('-s',
        dest='sample_annotations',
        required=True,
        help="Path to a sample annotation file.  No header line!"
    )

    parser.add_argument('-n', 
        dest='num_genes',
        default=20,
        type=int,
        help="The maximum number of differentially-expressed genes to plot (default: %(default)s)"
    )

    parser.add_argument('-p', 
        dest='padj_threshold',
        default=0.05,
        type=float,
        help="The \"significance threshold\" (adjusted p-value) for differential expression. (default: %(default)s)"
    )

    parser.add_argument('-x',
        dest='contrast_id',
        help= 'The name of the contrast, which makes the prefix for the output figures.'
    )

    parser.add_argument('-o',
        dest='output_dir',
        help= 'The name of the output directory where the figures are placed.'
    )
    
    parser.add_argument('-a',
        dest='base_condition',
        help= 'The base/control condition'
    )

    parser.add_argument('-b',
        dest='experimental_condition',
        help= 'The experimental condition'
    )

    args = parser.parse_args()
    return vars(args)


def jitter(n, mid, delta=0.2):
    return delta*np.random.random(n)-0.5*delta + mid


def make_scatter_plot_matrix(dge_df, nc, annotations, nmax, padj_threshold, contrast_id, output_dir):
    '''
    This makes a scatter plot for the top genes
    '''

    unique_conditions = annotations['condition'].unique()
    if len(unique_conditions) != 2:
        print('Should only be two groups.  Simple contrasts only!')
        sys.exit(1)
    
    group1_name = unique_conditions[0]
    group2_name = unique_conditions[1]
    group1_samples = annotations.loc[annotations.condition == group1_name]['sample'].values
    group2_samples = annotations.loc[annotations.condition == group2_name]['sample'].values

    # set some plot params:
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.size'] = 14
    ncols=3

    # get the number of significantly diff. exp. genes.
    N = np.min([nmax, np.sum(dge_df['padj']<=padj_threshold)])
    if N == 0:
        print('No differentially expressed genes at the padj <= %s threshold' % padj_threshold)
        fig, ax = plt.subplots(figsize=(9,9))
        left, width = .25, .5
        bottom, height = .25, .5
        right = left + width
        top = bottom + height
        ax.text(0.5*(left+right), 0.5*(bottom+top), 'Zero genes passed\nsignificance at\npadj<%.2f threshold' % padj_threshold,
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=20, color='red',
                transform=ax.transAxes)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        plt.tight_layout()
        fig.savefig(os.path.join(output_dir, '%s.scatter_plot.pdf' % contrast_id), bbox_inches='tight')
        fig.savefig(os.path.join(output_dir, '%s.scatter_plot.png' % contrast_id), bbox_inches='tight')
        plt.close()
        return
    nrows = int(np.ceil(N/ncols))
    fig, axarray = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20,5*nrows))
    for i in range(N):
        r = i // ncols
        c = i % ncols
        ax = axarray[r,c]
        dge_row = dge_df.iloc[i]
        gene_name = dge_row['Gene']
        counts = nc.loc[gene_name]
        group1_counts = counts[group1_samples]
        group2_counts = counts[group2_samples]
        ax.scatter(
            jitter(group1_counts.shape[0], 0.0),
            group1_counts,
            alpha=0.5,
            s=100
        )
        ax.scatter(
            jitter(group2_counts.shape[0], 1.0),
            group2_counts,
            alpha=0.5,
            s=100
        )
        ax.set_xticks([0,1])
        ax.set_xlim([-0.25, 1.25])
        ax.set_xticklabels([group1_name, group2_name])
        ax.set_title('%s (p=%.2e, p-adj=%.2e)' % (gene_name, dge_row['pvalue'], dge_row['padj']))

    # now remove any excess empty axes:
    total_axes = nrows*ncols
    if total_axes > N:
        index = N
        while index < total_axes:
            r = index // ncols
            c = index % ncols
            fig.delaxes(axarray[r,c])
            index += 1

    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, '%s.scatter_plot.pdf' % contrast_id), bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, '%s.scatter_plot.png' % contrast_id), bbox_inches='tight')


def make_volcano_plot(dge_df, padj_threshold, contrast_id, output_dir):
    '''
    Makes a basic, non-interactive volcano plot
    '''
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.size'] = 14
    fig, ax = plt.subplots(figsize=(10,10))
    min_fc = dge_df['log2FoldChange'].min()
    max_fc = dge_df['log2FoldChange'].max()
    sig = dge_df['padj'] <= padj_threshold
    scatter1 = ax.scatter(dge_df.loc[sig, 'log2FoldChange'], -np.log10(dge_df.loc[sig,'padj']), alpha=0.5)
    scatter2 = ax.scatter(dge_df.loc[~sig, 'log2FoldChange'], -np.log10(dge_df.loc[~sig,'padj']), alpha=0.5)
    ax.hlines(-np.log10(padj_threshold), min_fc, max_fc)
    ax.set_xlim([min_fc, max_fc])
    ax.set_xlabel('$\log_2$(Fold-change)', fontsize=20)
    ax.set_ylabel('$-\log_{10}$(Adjusted p-value)', fontsize=20)
    #sig_gene_names = dge_df.loc[sig, 'Gene'].tolist()
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, '%s.volcano_plot.pdf' % contrast_id), bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, '%s.volcano_plot.png' % contrast_id), bbox_inches='tight')


def interactive_volcano(dge_df, padj_threshold, contrast_id, output_dir):
    '''
    Makes a dynamic plot using Bokeh
    '''

    # split the data to sig and insig:
    sig_row_idx = dge_df['padj'] <= padj_threshold
    sig_rows = dge_df.loc[sig_row_idx]
    unsig_rows = dge_df.loc[~sig_row_idx]

    # subsample the insignificant ones so we don't send too much to the front-end
    fraction = 0.1
    n_sample = int(unsig_rows.shape[0]*fraction)
    keep_idx = np.random.permutation(~sig_row_idx)[:n_sample]
    unsig_rows = unsig_rows.loc[keep_idx]

    # extract the data and make a ColumnDataSource for Bokeh
    sig_x_vals = sig_rows['log2FoldChange'].values
    unsig_x_vals = unsig_rows['log2FoldChange'].values
    padj_vals = sig_rows['padj'].values
    sig_y_vals = -np.log10(sig_rows['padj'].values)
    unsig_y_vals = -np.log10(unsig_rows['padj'].values)
    sig_gene_names = sig_rows['Gene'].values
    sig_source = ColumnDataSource(data=dict(x=sig_x_vals, y=sig_y_vals, gene=sig_gene_names, padj=padj_vals))
    unsig_source = ColumnDataSource(data=dict(x=unsig_x_vals, y=unsig_y_vals))

    # behavior of the hover:
    hover = HoverTool(
        tooltips=[
            ('Gene', "@gene"),
            ("log2 fold-change", "@x{%.2f}"),
            ('pval (adj)','@padj{%.2e}')
        ],
        formatters = {
            'padj': 'printf',
            'x': 'printf',
        },
        names=['sig']
    )
    p = figure(plot_width=1500, 
        plot_height=800, 
        tools=[hover, 'pan', 'wheel_zoom','box_zoom','reset', 'save'], 
        title='Differentially-expressed genes',
        x_axis_label='log2(fold-change)', y_axis_label='-log10(Adjusted p-value)')

    p.circle('x', 'y', size=10, name='sig', source=sig_source, fill_alpha=0.5)
    p.circle('x', 'y', size=10, source=unsig_source, fill_alpha=0.5)

    output_file(os.path.join(output_dir, '%s.volcano.html' % contrast_id))
    save(p)


if __name__ == '__main__':
    args = get_arguments()
    dg_table = pd.read_table(args['dg_table'])
    nc_table = pd.read_table(args['counts_table'], index_col=0)
    annotations = pd.read_table(args['sample_annotations'], header=None, names=['sample', 'condition'])
    
    # subset the annotations for this contrast:
    conditions = np.array([args['base_condition'], args['experimental_condition']])
    annotations = annotations.loc[annotations.condition.isin(conditions)]
    
    make_scatter_plot_matrix(
        dg_table, 
        nc_table, 
        annotations, 
        args['num_genes'], 
        args['padj_threshold'],
        args['contrast_id'],
        args['output_dir']
    )

    #make_volcano_plot(dg_table, args['padj_threshold'], contrast_id, output_dir)
    interactive_volcano(
        dg_table, 
        args['padj_threshold'],
        args['contrast_id'],
        args['output_dir']
    )
