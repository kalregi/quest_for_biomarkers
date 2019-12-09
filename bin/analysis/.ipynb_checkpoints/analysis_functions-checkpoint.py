#!/usr/bin/env python
# -*-coding: utf8 -*-
""".
Module for helper functions used in the analysis process.

Copyright 2017 by Attila Jády, Balázs Ligeti (obalasz@gmail.com). All rights reserved.
ds dgf

"""
import pandas
import numpy as np
import math
from math import floor


def get_filtered_taxon_sample_table(species_cum, tax_data_pure, min_pct_hit_count=0.0001):
    """ Filtering the taxa based on the raw hit counts. Taxa recieveing lower percentage of a hits will be excluded"""
    print("Filtering the taxa using pct treshold")
    from math import floor

    sc = species_cum.stack().to_frame()
    sc.reset_index(inplace=True)
    new_columnnames = ['taxon_id', 'sample_id', 'hitcount']
    sc.columns = new_columnnames
    sc.columns = new_columnnames
    taxon2hitcount = sc.merge(tax_data_pure, how='left', left_on='taxon_id', right_index=True)
    root2maxhit_count = taxon2hitcount[taxon2hitcount["taxon_id"] == 1][["sample_id", "hitcount"]]
    root2maxhit_count["min_count"] = root2maxhit_count.apply(lambda row: floor(row[1] * min_pct_hit_count), axis=1)
    taxon2hitcount_min_count = taxon2hitcount.merge(root2maxhit_count[["sample_id", "min_count"]], left_on='sample_id', right_on='sample_id')
    taxon2hitcount_min_count = taxon2hitcount_min_count[taxon2hitcount_min_count["hitcount"] >= taxon2hitcount_min_count["min_count"]]
    return [taxon2hitcount, taxon2hitcount_min_count]


def normalize_taxahits_taxonomy_level(taxon2hitcount, taxon2hitcount_min_count):
    """ Normilizing the the hit counts at each taxonomy level i.e. species, genus etc"""
    sample_rank2count_sum = taxon2hitcount.groupby(['sample_id', 'rank']).agg('sum')['hitcount'].reset_index()
    sample_taxon2hitcounts_normalized_count = taxon2hitcount_min_count.merge(sample_rank2count_sum, how='left', left_on=['sample_id', 'rank'], right_on=['sample_id', 'rank'])
    sample_taxon2hitcounts_normalized_count['normalized_count'] = sample_taxon2hitcounts_normalized_count['hitcount_x'].divide(sample_taxon2hitcounts_normalized_count['hitcount_y'])
    sample_taxon2hitcounts_normalized_count = sample_taxon2hitcounts_normalized_count.rename(columns={'hitcount_y': 'all_hitcount', 'hitcount_x': 'hitcount'})
    # Calculating how many times a taxon appear in the samples
    taxon2_avg_count = sample_taxon2hitcounts_normalized_count.groupby(['taxon_id', 'taxon_name', 'rank']).agg(['mean', 'count'])['normalized_count'].reset_index()
    return [sample_taxon2hitcounts_normalized_count, taxon2_avg_count]


def get_taxon_sample_dataframe(sample_taxon2hitcounts_normalized_count, taxon2_avg_count,
                               Npresent=20, avg_abd=0, rank='species',
                               dominant_abd=0.04, dominant_Nabd=20, dominant_phages_only=False, mic_factor=10e-45):
    """    Constructs a pivoted data table based on the given filtering paramaters
    
    """
    print("AAAHeldsdsadlo4dffd")
    # import pandas
    # import numpy as np
    filter_conditions = (taxon2_avg_count['count'] > Npresent) & \
                        (taxon2_avg_count['mean'] > avg_abd) & \
                        (taxon2_avg_count['rank'] == rank)
    taxon2_avg_count_filtered = taxon2_avg_count[filter_conditions]
    # filtering to dominant taxa: i.e. they present at least in dominant_Nabd taxa with abundance dominant_abd
    st2nrc = sample_taxon2hitcounts_normalized_count[sample_taxon2hitcounts_normalized_count["normalized_count"] > dominant_abd]
    st2nrc_tsh_counted = st2nrc[["taxon_id", "sample_id"]].groupby(["taxon_id"]).agg('count').reset_index()
    dominant_s = st2nrc_tsh_counted[st2nrc_tsh_counted['sample_id'] > dominant_Nabd]  # domináns fágok
    if dominant_phages_only:
        taxon2_avg_count_filtered = taxon2_avg_count_filtered[taxon2_avg_count_filtered.taxon_id.isin(dominant_s['taxon_id'])]

    N_distinct_taxa = len(taxon2_avg_count_filtered['taxon_id'].unique())
    print('number of taxa: ', N_distinct_taxa)

    taxon_id_sample_normalized_count = sample_taxon2hitcounts_normalized_count.merge(taxon2_avg_count_filtered, how='inner', left_on='taxon_id', right_on='taxon_id')[
        ['sample_id', 'taxon_name_x', 'normalized_count']].drop_duplicates()
    taxon_id_sample_normalized_count_pivot = taxon_id_sample_normalized_count.pivot(index='taxon_name_x', columns='sample_id', values='normalized_count').fillna(0)
    taxon_id_sample_normalized_count_pivot = taxon_id_sample_normalized_count_pivot[taxon_id_sample_normalized_count_pivot.notnull()]

    random_numbers = pandas.DataFrame(mic_factor * np.random.rand(taxon_id_sample_normalized_count_pivot.shape[0], taxon_id_sample_normalized_count_pivot.shape[1]),
                                      columns=taxon_id_sample_normalized_count_pivot.columns, index=taxon_id_sample_normalized_count_pivot.index)
    taxon_id_sample_normalized_count_pivot = taxon_id_sample_normalized_count_pivot.add(random_numbers)
    # print(taxon_id_sample_normalized_count_pivot.head())
    return taxon_id_sample_normalized_count_pivot


def plot_cluster_map_with_classlabel(classlabel, data_table, data_labels, all_colors, pdffile=None, fig_size_y_per_entitity=0.3):
    """ Plot SNS cluster map with labels"""
    import seaborn as sns;
    sns.set()
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    res = 300
    ## adding legend information for the better readibility

    # figsize=[fig_size_x, len(tt.index)*fig_size_y_per_entitity]
    target_name2color = {}
    used_colors = []
    values = []
    indeces = []
    y = data_labels[classlabel].copy()
    # print(y)
    target_names = list(y.unique())
    print("target_names,", target_names)
    for target_name, color in zip(target_names, all_colors[:len(target_names)]):
        y[y == target_name] = color
        target_name2color[target_name] = color
    print(target_name2color)
    for k, v in target_name2color.items():
        used_colors.append(v)
        values.append(1)
        indeces.append(k)
    legend_data = pandas.DataFrame(values, index=indeces).transpose()
    f = plt.figure()
    plt.title('Colors', color='black')
    legend_data.plot.barh(colors=used_colors, ax=f.gca())
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    # plt.tight_layout()
    if pdffile is not None:
        pdffile.savefig(f, bbox_inches='tight', dpi=res)

    g = sns.clustermap(data_table, metric='braycurtis', figsize=[16, len(data_table.index) * fig_size_y_per_entitity], col_colors=y, xticklabels=False)
    plt.setp(g.ax_heatmap.set_xlabel(''))
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    # plt.show()
    # plt.tight_layout()
    if pdffile is not None:
        pdffile.savefig(g.fig, bbox_inches='tight', dpi=res)

    return target_name2color, g


def plot_pca(data, data_labels, data_label, colors, pca_components):
    from sklearn.decomposition import PCA, KernelPCA
    import matplotlib.pyplot as plt
    lw = 2
    y = data_labels[data_label]
    target_names = list(y.unique())
    pca = PCA(n_components=pca_components)
    pca.fit(data)
    # print(pca.explained_variance_ratio_)
    plt.figure(1)
    # print(pca.components_)
    if pca_components == 2:
        for color, target_name in zip(colors, target_names):
            print(target_name)
            X_sub = data[y == target_name]
            X_sub_proj = pca.transform(X_sub)
            plt.scatter(X_sub_proj[:, 0], X_sub_proj[:, 1], color=color, alpha=.4, lw=lw, label=target_name)
        plt.legend(loc='best', shadow=False, scatterpoints=1)
        plt.title('PCA of Phage')
        plt.show()

    elif pca_components == 3:
        for color, target_name in zip(colors, target_names):
            print(target_name)
            X_sub = data[y == target_name]
            X_sub_proj = pca.transform(X_sub)

            plt.subplot(221)
            plt.scatter(X_sub_proj[:, 0], X_sub_proj[:, 1], color=color, alpha=.3, lw=lw, label=target_name)
            plt.legend(loc='best', shadow=False, scatterpoints=1)
            plt.title('PCA coord 0 vs PCA coord 1')

            plt.subplot(222)
            plt.scatter(X_sub_proj[:, 0], X_sub_proj[:, 2], color=color, alpha=.3, lw=lw, label=target_name)
            plt.legend(loc='best', shadow=False, scatterpoints=1)
            plt.title('PCA coord 0 vs PCA coord 2')

            plt.subplot(223)
            plt.scatter(X_sub_proj[:, 1], X_sub_proj[:, 2], color=color, alpha=.3, lw=lw, label=target_name)
            plt.legend(loc='best', shadow=False, scatterpoints=1)
            plt.title('PCA coord 1 vs PCA coord 2')
        plt.show()
    else:
        print("to many coordinates, discard plotting! ")

    return pca


def plot_kernel_pca(data, data_labels, data_label, colors, pca_components):
    '''
    So many bugs in that particular function, should be restored or checked!
    :param data:
    :param data_labels:
    :param data_label:
    :param colors:
    :param pca_components:
    :return:
    '''
    from sklearn.decomposition import PCA, KernelPCA
    lw = 2
    y = data_labels[data_label]
    target_names = list(y.unique())
    # pca = PCA(n_components=pca_components)
    # pca.fit(data)

    kpca = KernelPCA(kernel="rbf", fit_inverse_transform=True, gamma=10)
    X_kpca = kpca.fit_transform(data)

    # print(kpca.explained_variance_ratio_)
    plt.figure(1)
    if pca_components == 2:
        for color, target_name in zip(colors, target_names):
            print(target_name)
            X_sub = data[y == target_name]
            X_sub_proj = kpca.transform(X_sub)
            plt.scatter(X_sub_proj[:, 0], X_sub_proj[:, 1], color=color, alpha=.3, lw=lw, label=target_name)
        plt.legend(loc='best', shadow=False, scatterpoints=1)
        plt.title('Kernel PCA of Phage')
        plt.show()

    elif pca_components == 3:
        for color, target_name in zip(colors, target_names):
            print(target_name)
            X_sub = data[y == target_name]
            X_sub_proj = kpca.transform(X_sub)

            plt.subplot(221)
            plt.scatter(X_sub_proj[:, 0], X_sub_proj[:, 1], color=color, alpha=.3, lw=lw, label=target_name)
            plt.legend(loc='best', shadow=False, scatterpoints=1)
            plt.title('Kernel PCA coord 0 vs PCA coord 1')

            plt.subplot(222)
            plt.scatter(X_sub_proj[:, 0], X_sub_proj[:, 2], color=color, alpha=.3, lw=lw, label=target_name)
            plt.legend(loc='best', shadow=False, scatterpoints=1)
            plt.title('Kernel PCA coord 0 vs PCA coord 2')

            plt.subplot(223)
            plt.scatter(X_sub_proj[:, 1], X_sub_proj[:, 2], color=color, alpha=.3, lw=lw, label=target_name)
            plt.legend(loc='best', shadow=False, scatterpoints=1)
            plt.title('Kernel coord 1 vs PCA coord 2')
        plt.show()
    else:
        print("to many coordinates!")


def generate_string_from_params(**kwargs):
    """ Generates a string from paramters"""
    print(kwargs)
    description_str = ""
    for k, v in kwargs.items():
        description_str += "." + str(k) + "-" + str(v)

    return description_str


def pca_3d_plotly(data, data_labels, class_name, colors, pca_components):
    from sklearn.decomposition import PCA, KernelPCA
    from mpl_toolkits.mplot3d import Axes3D
    import plotly.plotly as py
    import plotly.graph_objs as go
    from plotly.graph_objs import Scatter, Figure, Layout
    from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
    import plotly
    # from jupyterlab_plotly import Plotly

    plotly.offline.init_notebook_mode()

    lw = 2
    y = data_labels[class_name]
    target_names = list(y.unique())
    pca = PCA(n_components=pca_components)
    pca.fit(data)
    plot_data = []
    for color, target_name in zip(colors, target_names):
        X_sub = data[y == target_name]
        X_sub_proj = pca.transform(X_sub)

        trace_act = go.Scatter3d(
            x=X_sub_proj[:, 0],
            y=X_sub_proj[:, 1],
            z=X_sub_proj[:, 2],
            mode='markers',
            name=target_name,
            marker=dict(
                size=8,
                line=dict(
                    color=color,
                    width=0.5
                ),
                opacity=0.8
            )
        )
        plot_data.append(trace_act)
    layout = go.Layout(
        showlegend=True,
        margin=dict(
            l=0,
            r=0,
            b=0,
            t=0
        )
    )
    fig = go.Figure(data=plot_data, layout=layout)
    iplot(fig, filename='simple-3d-scatter')
    return plot_data, layout
    # Plotly(fig)


##  Function for transforming the ids
def bac_name_transform(act_name):
    act_name = act_name.lower()
    act_name = act_name.replace('-', ' ', 10)
    act_name = act_name.replace('_', ' ', 10)
    act_name = act_name.replace('/', ' ', 10)
    act_name = act_name.replace(r'\\', ' ', 10)
    act_name = act_name.replace(r':', ' ', 10)
    act_name = act_name.replace(r"'", '', 10)
    act_name = act_name.replace(r",", '', 10)
    act_name = act_name.replace(r"\.", '', 10)
    return act_name


# get the ncbi id from the marker id
# It could be GeneID or refseq, genbank or WGS accession

def get_ncbi_id(metaphlan_id):
    # print("Processing ....")
    # print(metaphlan_id)
    prefix = metaphlan_id.split('|')
    ncbi_id = None
    if len(prefix) == 5 and ('|ref|' in metaphlan_id or '|gb|'):
        ncbi_id = prefix[3]
    elif 'NC' in metaphlan_id:
        ncbi_id = metaphlan_id[0]
    elif 'GeneID:' in metaphlan_id:
        ncbi_id = metaphlan_id[7:]
    else:
        ncbi_id = None
    # print(ncbi_id)
    return ncbi_id


def get_ncbi_taxon_id_entrez(gene_id, database='Gene'):
    from Bio import Entrez
    import re

    regexp = re.compile(r"(taxon\:)(\d+)")

    Entrez.email = "ligeti.balazs@itk.ppke.hu"
    handle = Entrez.efetch(db=database, id=str(gene_id), retmode="xml")
    if '.' in gene_id:
        gene_id = gene_id.split('.')[0]

    # records = Entrez.parse(handle)
    # records
    taxon_id = None
    for record in Entrez.parse(handle):
        # print(record)
        if gene_id not in str(record):
            return None

        if database == 'Gene':
            refdbrecord = record['Entrezgene_source']['BioSource']['BioSource_org']['Org-ref']['Org-ref_db']
            dbdatai = 0
            for dbdata in refdbrecord:
                if dbdata['Dbtag_db'] == 'taxon':
                    taxon_id = dbdata['Dbtag_tag']['Object-id']['Object-id_id']
                    break
        else:
            if 'taxon:' in str(record):
                # print("There is taxon information:")
                mo = regexp.search(str(record))
                if mo is not None:
                    #   print("MATCH")
                    print(mo.group(2))
                else:
                    pass
        #           print("No match")

    return taxon_id
