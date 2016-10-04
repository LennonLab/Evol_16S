from __future__ import division
from skbio import TabularMSA, RNA
import os
import skbio as sk

mydir = os.path.expanduser("~/github/Evol_16S/")



import matplotlib.pyplot as plt

def msa_to_heatmap(msa, value_map, legend_labels=('Low', 'Medium', 'High'), fig_size=(15,10), cmap='YlGn', sequence_order=None):
    """Plot a multiple sequence alignment as a heatmap.

    Parameters
    ----------
    msa : skbio.TabularMSA
        The multiple sequence alignment to be plotted
    value_map : dict, collections.defaultdict
        Dictionary mapping characters in the alignment to values. KeyErrors are not
        caught, so all possible values should be in this dict, or it should be a
        collections.defaultdict which can, for example, default to ``nan``.
    legend_labels : iterable, optional
        Labels for the min, median, and max values in the legend.
    fig_size : tuple, optional
        Size of figure in inches.
    cmap : matplotlib colormap, optional
        See here for choices: http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
    sequence_order : iterable, optional
        The order, from top-to-bottom, that the sequences should be plotted in.

    Raises
    ------
    KeyError
        If a character in ``msa`` is not in ``value_map``, and ``value_map`` is not a
        ``collections.defaultdict``.

    """
    if sequence_order is None:
        sequence_order = msa.index

    # fill a data matrix by iterating over the alignment and mapping
    # characters to values
    mtx = []
    for label in sequence_order:
        seq = str(msa.loc[label])
        mtx.append([value_map[aa] for aa in seq])

    # build the heatmap, this code derived from the Matplotlib Gallery
    # http://matplotlib.org/examples/pylab_examples/colorbar_tick_labelling_demo.html
    fig, ax = plt.subplots()
    fig.set_size_inches(fig_size)

    cax = ax.imshow(mtx, interpolation='nearest', cmap=cmap)

    # Add colorbar and define tick labels
    values = list(value_map.values())
    cbar = fig.colorbar(cax,
                        ticks=[min(values),
                               np.nanmedian(values),
                               max(values)],
                        orientation='horizontal')
    ax.set_yticks([0] + list(range(3, msa.shape.sequence - 3, 3)) + [msa.shape.sequence - 1])
    ax.set_yticklabels(sequence_order)
    ax.set_xticks(range(msa.shape.position))
    ax.set_xticklabels(msa.consensus(), size=7)
    cbar.ax.set_xticklabels(legend_labels) # horizontal colorbar
    fig.savefig(mydir + 'figs/align.png', bbox_inches='tight',  dpi = 600)


msa = TabularMSA.read(mydir + 'data/nmicrobiol201648-s7.txt', constructor=RNA)
msa.reassign_index(minter='id')


import numpy as np
from collections import defaultdict
hydrophobicity_idx = defaultdict(lambda: np.nan)
hydrophobicity_idx.update({'A': 0.61, 'C': 0.61, 'G': 0.61, 'U': 0.61,})
hydrophobicity_labels=['Hydrophilic', 'Medium', 'Hydrophobic']

msa_to_heatmap(msa, hydrophobicity_idx, legend_labels=hydrophobicity_labels)

#msa = sk.alignment.TabularMSA.read(mydir + 'data/nmicrobiol201648-s7.txt', format = 'fasta')
#msa.reassign_index(minter='id')


#msa = TabularMSA.read('data/nmicrobiol201648-s7.txt', constructor=Protein)
#msa.reassign_index(minter='id')
