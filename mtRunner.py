import os, re, shutil, subprocess
import multiprocessing as mp
import pandas as pd, numpy as np
from multiprocessing.dummy import Pool
from functools import partial
from subprocess import CalledProcessError

"""
Runs IUPred and CAST in parallel on a given CSV file. This CSV file
must have the sequence identities (like Uniprot IDs) in the first
column and contain a 'Sequence' column. 
"""

prog_commands = {
    'iupred':'./iupred ../{0} "long"',
    'cast':'./cast/cast_MacOSX {0} -verbose'
}

def process_cast(cast_output, seq):
    """
    Extracts all relevant information from cast output

    Takes the cast output and the original sequence and extracts the all the 
    information for each region determined by cast. It gets the enriched AA, 
    the start, the end, the score, and the positions of the AAs belonging to 
    the region. 

    Parameters:
    ----------

    cast_output: str
        The raw output from the cast algorithm
    seq: str
        The original sequence

    Returns:
    -------
    out: list of ints, str
        For each region, the enriched AA, the start, the end, the CAST score,
        followed by a string representation of a list of positions of each
        residue belonging to the region. 
        
    """
    # Compiling the regexes gives a trivial speedup in test cases, but makes the code more readable 
    get_masked_seq = re.compile(r'([A-Z]{2,})')
    get_region_info = re.compile(r'([A-Z])-rich.*?([\d]+)\s.*?([\d]+)\s.*?([\d]+)')
    
    # exit quickly if cast didn't find anything
    if 'region' not in cast_output: return [np.nan] 
    
    # Clean up strings
    seq_stripped = seq.replace(' ', '')
    output = cast_output.replace('\n', '')
    masked_seq = get_masked_seq.findall(output)[0]
    
    # Get the CAST region info: the enriched AA, start, end, CAST score, and add an empty list for positions
    region_info = [element for region in get_region_info.findall(cast_output) 
                   for element in [region[0], int(region[1]), int(region[2]), int(region[3]), []]]
    
    # Bin the positions of masked amino acids into the enriched AA groups
    # This corrects a bug resulting from the presence of multiple regions of the same AA
    for i in range(len(masked_seq)):
        aa = seq_stripped[i]
        if masked_seq[i] is not aa: 
            # Let's find out what region this AA belongs to and add its position to that region's list
            for j in range(0, len(region_info), 5): 
                if region_info[j] is aa and region_info[j+1] <= i+1 and region_info[j+2] >= i+1:
                    region_info[j+4].append(i+1) 
    # convert the position lists to strings for processing simplicity
    region_info[4::5] = [str(positions) for positions in region_info[4::5]]
    return region_info

def process_iupred(iupred_output, seq):
    """
    Extracts relevant information from IUPred output. 
    """

def run(params, prog='iupred', tmpdir='disordered_tmps'):
    seq, seq_name = params

    tmpfile = "%s/tmp_%s.fasta"%(tmpdir, seq_name)
    f = open(tmpfile, 'w')
    f.write('>%s\n%s\n'%(seq_name, seq))
    f.close()

    if prog is 'iupred':
        output = subprocess.check_output(prog_commands[prog].format(tmpfile), 
            shell=True, cwd=prog, stderr=subprocess.STDOUT)
        #print(output)
    elif prog is 'cast':
        try:
            subprocess.check_output(prog_commands[prog].format(tmpfile).split(' '))
        except CalledProcessError as e:
            print(process_cast(str(e.output, 'utf-8'), seq))
    
    # results.append(pd.Series(probs_list if len(probs_list) > 0 else [np.nan]))

def runAnalysis(input_file, prog):

    polyprots = pd.read_csv(input_file, index_col=0)


    tmp_dir = 'disordered_tmps'
    os.environ['IUPred_PATH'] = os.getcwd() + '/iupred'

    seq_pairs = zip(polyprots['Sequence'].values.tolist(), 
        polyprots.index.values.tolist())
    
    try: 
        os.makedirs(tmp_dir)
    except OSError:
        if not os.path.isdir(tmp_dir):
            raise

    pool = Pool(max(mp.cpu_count() - 1, 1))

    for i, returncode in enumerate(pool.imap_unordered(partial(run, prog=prog, 
      tmpdir=tmp_dir), seq_pairs)):
        pass
    
    # shutil.rmtree(tmp_dir)

