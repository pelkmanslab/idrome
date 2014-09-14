import pandas as _pd, numpy as _np, seaborn as _sns, matplotlib.pyplot as _plt
from numpy.lib.stride_tricks import as_strided as _ast
import mischelperfuncs as _mhf
from matplotlib.patches import Rectangle as _rect

"""
Miscellaneous functions for running and analyzing Hidden Markov Models.
"""

_aas = _mhf.get_aas()

def run_HMM_viterbi_map(model, proteome, window1, window2, label):
    """
    Run a given HMM on a proteome, filtering first for regions where the viterbi 
    output matches a given state label for at least a specified window length. 
    Then the posterior probabilities are calculated for these viterbi regions 
    and their sum over a second window length is computed. These values are 
    returned as a array that is the same length as the size of the proteome. A 
    perfect match would have a score equal to the value of `window2`.
    """
    similarity = []

    indices = {state.name : i for i,state in enumerate(model.states)}

    for index, row in proteome.iterrows():
        seq = [char for char in row['SEQ'] if char in _aas]

        viterbi_path = model.viterbi(seq)[1]

        viterbi_pred = _np.array([state[1].name 
                                 for state in viterbi_path[1:]]) == label
        starts = _np.nonzero(viterbi_pred & ~_np.roll(viterbi_pred, 1))[0]
        ends = _np.nonzero(viterbi_pred & ~_np.roll(viterbi_pred, -1))[0]
        sgtr_regions = [(i, j) for i, j in zip(starts, ends) if j+1 >= i + window1]

        ems = _np.exp(model.forward_backward(list(seq))[1])
        probs = ems/_np.sum(ems, axis=1)[:, _np.newaxis]
        out = probs[:, indices[label]]
        padding = _np.zeros(window2 - 1)
        padded_out = _np.concatenate([padding,out])
        strided = _ast(padded_out,shape = (len(padded_out) + 1 - window2, window2),
                      strides = padded_out.strides * 2)
        scores = [0]
        for i,j in sgtr_regions:
            scores.append(strided.sum(1)[i:j+1].max())
        similarity.append(max(scores))

    return similarity

def plot_HMM(model, proteome, window1, window2, label, background_label, uniprot_id, regions=[], name="", mutations=[], seq=""):
    """
    
    """
    if len(seq) == 0:
        seq = [char for char in proteome.loc[uniprot_id, 'SEQ'] if char in _aas]
    else:
        if type(seq) == str:
            seq = list(seq)
        seq = [char for char in seq if char in _aas]
    
    for pos,aa in mutations:
        print("%s --- %s ---> %s"%(seq[pos-1], pos, aa))
        seq[pos-1] = aa
        
    indices = {state.name : i for i,state in enumerate(model.states)}
    
    # TODO: Duplicate functionality, encapsulate in separate functions
    # Calculate viterbi regions
    viterbi_pred = _np.array([state[1].name for state in model.viterbi(seq)[1][1:]]) == label
    starts = _np.nonzero(viterbi_pred & ~_np.roll(viterbi_pred, 1))[0]
    ends = _np.nonzero(viterbi_pred & ~_np.roll(viterbi_pred, -1))[0]
    vitb_regions = [(i, j) for i, j in zip(starts, ends) if j+1 >= i + window1]
    
    # Calculate MAP regions
    map_pred = _np.array([state[1].name for state in 
                         model.maximum_a_posteriori(seq)[1][1:-1]]) == label
    
    # Calculate posterior probabilities
    ems = _np.exp(model.forward_backward(list(seq))[1])
    probs = ems/_np.sum(ems, axis=1)[:, _np.newaxis]
    out = probs[:, indices[label]]
    padding = _np.zeros(window2 - 1)
    padded_out = _np.concatenate([padding,out])
    strided = _ast(padded_out,shape = (len(padded_out) + 1 - window2, window2),
                  strides = padded_out.strides * 2)
    
    map_pred = out.T >= 0.5
    starts = _np.nonzero(map_pred & ~_np.roll(map_pred, 1))[0]
    ends = _np.nonzero(map_pred & ~_np.roll(map_pred, -1))[0]
    map_regions = [(i, j) for i, j in zip(starts, ends) if j+1 >= i + 5]
    
    scores = [0]
    for i,j in vitb_regions:
        scores.append(strided.sum(1)[i:j+1].max())

    colors = _sns.color_palette('deep', n_colors=6, desat=0.5)
    fig = _plt.figure(figsize=(14,6));
    _plt.axhline(y=1.1, c=colors[0], alpha=0.7, figure=fig)
    _plt.axhline(y=0.5, c=colors[1], alpha=0.7, figure=fig)
    _plt.axhline(y=1.16, c=colors[1], alpha=0.7, figure=fig)
    ax = fig.get_axes()[0]
    ax.set_xlim([1, len(seq)+1])
    ax.set_ylim([0,1.2])
    ax.set_ylabel(r'posterior probability $\gamma_k$')
    ax.set_xlabel(r'amino acid position $k$')
    ax.set_title(name if name != '' else uniprot_id)
    
    # Plot viterbi predicted regions
    for start, end in vitb_regions:
        ax.add_patch(_rect((start+1, 1.075), end-start+1, 0.05, 
                                 facecolor=colors[0], alpha=0.7))
    ax.text(len(seq) - 20,1.13, 'viterbi')
        
    # Plot map predicted regions
    for start, end in map_regions:
        ax.add_patch(_rect((start+1, 1.135), end-start+1, 0.05, 
                                 facecolor=colors[1], alpha=0.7))
#         axis.text(start, .81, "%s"%start, rotation='vertical', fontsize=7)
#         axis.text(end, .81, "%s"%end, rotation='vertical', fontsize=7)
#         print(''.join(seq[start-5:end+5]))
    ax.text(len(seq) - 20,1.19, 'map')

    for start, end, color, descrip in regions:
        ax.add_patch(_rect((start, .475), end-start+1, 0.05, 
                                 facecolor=colors[color], alpha=0.7))
        ax.text(start, .53, descrip)

    unique_states = 0
    for idx, state in enumerate(indices.keys()):
        if 'end' not in state and 'start' not in state and state not in background_label:
            _plt.plot(range(1, len(seq)+1), probs[:, indices[state]], 
                     c=colors[2+unique_states], alpha=0.7, figure=fig)
            unique_states += 1
    
    return fig
