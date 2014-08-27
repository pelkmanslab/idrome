from scipy.cluster.hierarchy import *
import warnings as _warnings
from collections import Counter as _counter
import numpy as _np, pandas as _pd

"""
Miscellaneous helper functions, mainly related to sequence processing
"""

_aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
        'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def retrieve_disordered_regions(x):
    """
    Retrieves the sequences of disordered regions in a given protein sequence
    (as a Pandas Dataframe via `apply`) and returns them as a single string 
    with distinct regions separated by an underscore. NaN is returned if a 
    protein has no disordered regions longer than the threshold.
    """
    if _pd.isnull(x['Disordered Threshold: 30']):
        return _np.nan
    
    regions = [tuple(region.split('_')) for region 
               in x['Disordered Threshold: 30'].split(';')]
    return '_'.join([x['Sequence'][int(i)-1:int(j)] for i,j in regions])
    
def calc_disordered_regions(limits, seq):
    """
    Returns the sequence of disordered regions given a string of
    starts and ends of the regions and the sequence.
    
    Example
    -------
    limits = 1_5;8_10
    seq = AHSEDQNAAANTH...
    
    This will return `AHSED_AAA`
    """
    seq = seq.replace(' ', '')
    
    regions = [tuple(region.split('_')) for region 
               in limits.split(';')]
               
    return '_'.join([seq[int(i)-1:int(j)] for i,j in regions])

def build_tree(links, threshold=20):
    """
    Takes a NumPy-style linkage matrix and a threshold for the minimum number of 
    children for a node to be considered a "leaf". It returns the root of the 
    tree and list of the `ClusterNode` "leaves". 
    
    A link to the parent of each `ClusterNode` is added, as well as a list of the
    "true leaf" children nodes and the "true leaf" children ids. "true leaf" is 
    different from the aforementioned "leaf" because the true leaves are the 
    original IDRs or sequences. "Leaves", on the other hand, are clusters containing
    at least the threshold amount of "true leaves".  
    """
    
    (root,nodes) = to_tree(links, rd=True)
    _warnings.warn("This whole thing is a super fragile hack. This is could break \
hard with future SciPy releases. There needs to be a better way to map the \
original data to the ClusterNodes. For now I'm relying on this line to get \
this mapping: \
https://github.com/scipy/scipy/blob/v0.14.0/scipy/cluster/hierarchy.py#L880")
    orig_data_to_cluster_ids = [node.get_id() for node in nodes[:len(links)+1]]
    assert all([_id in orig_data_to_cluster_ids for _id in range(0, len(links)+1)]), "Something \
funky may be going on. Check the first nodes in ``d`` are the \
original data https://github.com/scipy/scipy/blob/master/scipy/cluster/hierarchy.py#L880"
    assert all([node.get_count() == 1 for node in nodes[:len(links)+1]]), "Something \
funky may be going on. Check the counts of the first batch of nodes in ``d``. The counts should be 1 because \
they are the original data https://github.com/scipy/scipy/blob/master/scipy/cluster/hierarchy.py#L880"
    
    t_leaves = []

    # Yes, recursion is horribly inefficient in python. This should be feasible
    # for smaller trees. 
    def build_parents(parent, t, t_leaves):
        """
        Adds parent information to the hierarchy tree. ``t`` is the threshold
        for the minimum number of elements needed to be called a leaf. Leaf
        nodes will be appended to ``t_leaves`` 
        """
        l,r = parent.get_left(), parent.get_right()
        if l is not None:
            l.parent = parent
            build_parents(l, t, t_leaves)
        if r is not None:
            r.parent = parent
            build_parents(r, t, t_leaves)
        # ugh
        if getattr(l, 'count', None) < t and getattr(r, 'count', None) < t and parent.get_count() >= t:
            t_leaves.append(parent)

    def build_child_info(parent):
        """
        For each node adds information of the leaves of original data contained in the
        subtree with that node as the root.
        """
        l,r = parent.get_left(), parent.get_right()
        if l is None and r is None:
            parent.children = []
            return [parent], [parent.id]
        left_children, left_ids = build_child_info(l)
        right_children, right_ids = build_child_info(r)
        parent.children = left_children + right_children
        parent.child_ids = left_ids + right_ids
        return parent.children, parent.child_ids

    root.parent = root
    build_parents(root, threshold, t_leaves)
    build_child_info(root)
    all_nodes = nodes
    
    return root, t_leaves, all_nodes
    
def get_go_counts(entries, proteome):
    """
    
    """
    go_terms = proteome.loc[entries.str.extract('([\w]+)_'), 'Gene ontology (GO)']
    go_terms = go_terms[go_terms.notnull()].str.split('; ').tolist()
    return _counter([term for prot in go_terms for term in _np.unique(prot)])

def get_freqs(seq, pseudocount=0, verbose=True):
    """
    Determines the frequencies of amino acids in a string representation
    of a sequence.

    Parameters
    ----------
    seq : str
        A string of representation of a protein sequence
    pseudocount : float, default 0
        Add a pseudocount to all emission values. When calculating emission
        values for HMMs trained on limited data 0.01, i.e.
        1 count per 100 observations, is good default.
    verbose : bool, default True
        Print removed values and final counts (for pseudocounts)

    Returns
    -------
    dict
        A dictionary of the amino acids mapped to their frequencies

    See Also
    --------
    get_normed_freqs : calculate normalized frequencies from a list of seqs

    """
    counts = _counter(seq)
    for aa in list(counts.keys()):
        if aa not in _aas:
            removed = counts.pop(aa, -1)
            if verbose: print('Extraneous char %s occurred %s times; removed.'%(aa,
                  removed))
    aa_count = float(sum(counts.values()))
    freqs = {}
    for aa in _aas:
        freqs[aa] = counts[aa] + pseudocount * aa_count
    final_counts = float(sum(freqs.values()))
    for aa in _aas:
        freqs[aa] /= final_counts
    if verbose: print final_counts
    return freqs

def get_normed_freqs(seqs, pseudocount=0.01):
    """
    Determine the frequencies of amino acids in a given list of sequences. The
    frequencies will be calculated per-sequence and then averaged together so
    all sequences have an equal contribution to the final frequencies.

    Parameters
    ----------
    seqs : list of strings
        A list of sequences in string form
    pseudocount : float, default 0.01
        Add a pseudocount for missing emission values. Defaults to 0.01, i.e.
        1 count per 100 observations. Necessary if training a HMM on a limited
        number of sequences.

    Returns
    -------
    dict
        A dictionary of the amino acids mapped to their frequencies

    """
    seq_counts = [get_freqs(seq, pseudocount) for seq in seqs]
    num_seqs = float(len(seq_counts))

    normed_freqs = _counter()
    for aa in _aas:
        for seq_count in seq_counts:
            normed_freqs[aa] += seq_count[aa]
        normed_freqs[aa] /= num_seqs
    return normed_freqs

def get_aas():
    """
    Returns a list of the one-letter amino acid codes

    """
    return _aas
    
def get_percent_enrichments(new, old):
    """
    Calculate the change in frequencies of amino acids between two groups as
    their percent enrichment in the `new` group as compared to the `old` group.
    
    Parameters
    ----------
    new : dict
        A dictionary of amino acid frequencies for a region of interest
    old : dict
        A dictionary of amino acid frequencies to compare against.
    
    Returns
    -------
    dict
        A dictionary of the amino acids mapped to their percent enrichment in
        `new` versus `old`
        
    See Also
    --------
    get_fold_enrichments : calculate fold enrichment between two groups

    """
    return {aa:new[aa]/old[aa]*100 for aa in _aas}

def get_fold_enrichments(new, old):
    """
    Calculate amino acid fold enrichment between two groups. Enriched and
    depleted values are normalized around 1 which makes comparisons between them
    more straightforward.
    
    Parameters
    ----------
    new : dict
        A dictionary of amino acid frequencies for a region of interest
    old : dict
        A dictionary of amino acid frequencies to compare against.
    
    Returns
    -------
    dict
        A dictionary of the amino acids mapped to their fold enrichment in
        `new` versus `old`. Values are centered around 1.
    
    """
    return {aa:new[aa]/old[aa] if new[aa] > old[aa] else -1*old[aa]/new[aa] for aa in _aas}