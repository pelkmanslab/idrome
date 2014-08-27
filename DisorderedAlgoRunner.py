# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 12:20:02 2014

A Python wrapper around IUPred and CAST. The IUPred and CAST
executables should be located in the iupred/ and cast/
folders, respectively

@author: Tamas Nagy <tamas at tamasnagy dot com>
"""
import os, sys, subprocess, warnings, re
import pandas as pd, numpy as np

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
    for i in xrange(len(masked_seq)):
        aa = seq_stripped[i]
        if masked_seq[i] is not aa:
            # Let's find out what region this AA belongs to and add its position to that region's list
            for j in xrange(0, len(region_info), 5):
                if region_info[j] is aa and region_info[j+1] <= i+1 and region_info[j+2] >= i+1:
                    region_info[j+4].append(i+1)
    # convert the position lists to strings for processing simplicity
    region_info[4::5] = [str(positions) for positions in region_info[4::5]]
    return region_info


def IUPredRunner(seqs, seq_names, filename):
    """
    Runs IUPred on a list of sequences and saves the output
    to the given filename.
    """
    assert(len(seqs) == len(seq_names))
    if len(seqs) > 5000:
        warnings.warn("Running IUPred sequentially on >5000 sequences takes awhile. Memory pressure" +\
            " will also be quite high because the entire input and output datasets are in memory simultaneously")
        sys.stderr.flush()
    
    print("Running IUPRED on %s sequences..."%len(seqs))
    sys.stdout.flush()
    
    probs = pd.DataFrame()
    
    results = []
    
    try:
        os.chdir('iupred/')
        tmpfile = 'tmpfile.fasta'
        
        # Set the path to IUPred
        os.environ['IUPred_PATH'] = os.getcwd()
        
        # Hackish way of extracting information from IUPred output
        f = lambda x, seq: [float(a.split('     ')[1]) for a in str(x).rsplit('#', 1)[1].split('\n')[1:-1]]
        
        command = lambda: subprocess.check_output('./iupred %s "long"'%tmpfile, shell=True)
        
        results = getRawOutput(seqs, tmpfile, command, f)
        
        # Glue all the columns together
        probs = pd.concat(results, axis=1).T
        probs.index = seq_names
            
    except subprocess.CalledProcessError as e:
        print(e.output)
    except KeyboardInterrupt:
        print("Exiting...")
    finally:
        os.chdir('../')
        probs.to_csv(filename)
        
def CastRunner(seqs, seq_names, filename):
    """
    Runs CAST on the given list of sequences and saves the output to the
    given filename.
    """
    assert(len(seqs) == len(seq_names))
    if len(seqs) > 10000:
        warnings.warn("Running cast sequentially on >10000 sequences takes awhile. Memory pressure" +\
            " will also be quite high because the entire input and output datasets are in memory simultaneously")
        sys.stderr.flush()
    
    print("Running CAST on %s sequences..."%len(seqs))
    sys.stdout.flush()
    
    probs = pd.DataFrame()
    
    results = []
    
    try:
        os.chdir('cast/')
        tmpfile = 'tmp.fasta'
        
        # f = lambda x : [b.split(' ')[c] for b in x.splitlines() if 'region' in b for c in [0, 3, 5]]
        
        command = lambda: subprocess.Popen(
        ('./cast_MacOSX %s -verbose'%tmpfile).split(' '), stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE).communicate()[0]
        
        # TODO: Test this on a Linux machine, hopefully this won't break.
        results = getRawOutput(seqs, tmpfile, command, process_cast)
        
        # Glue all the columns together
        probs = pd.concat(results, axis=1).T
        probs.index = seq_names
    
    except subprocess.CalledProcessError as e:
        print(e.output)
    except KeyboardInterrupt:
        print("Exiting...")
    finally:
        os.chdir('../')
        probs.to_csv(filename)
        
def runGBA(seqs, seq_names, filename):
    """
    Experimental
    """
    assert(len(seqs) == len(seq_names))

    warnings.warn("This is experimental and broken. Don't use", FutureWarning)
    
    probs = pd.DataFrame()
    
    results = []
    
    with open('%s.fasta'%filename, 'w') as f:
        for seq, name in zip(seqs, seq_names):
            f.write('>%s\n%s\n\n'%(name, seq))


def runDisorderedAnalysis(input_file, runCAST=True, runIUPred=True, forceIUPred=False, forceCAST=False):
    """
    Runs a suite of disorder prediction algorithms (e.g. IUPred, CAST) on the
    given input_file. This input_file should be a csv and have a column named
    `Sequences` and a column named `Entry`. The force options, if set to true,
    will re-run the analysis and override previous LCR and IDR information
    present in the csv, respectively.
    """

    polyprots = pd.read_csv(input_file, index_col=0)

    if all(['LCRs' in polyprots.columns, 'IDRs' in polyprots.columns, not forceIUPred, not forceCAST]):
        print('Nothing to do. Use `force=True` to force generation.')
        return

    if ('LCRs' not in polyprots.columns and runCAST) or forceCAST:

        cast_output_file = "%s_cast.csv"%input_file.rsplit('.csv', 1)[0]

        if forceCAST:
            if 'LCRs' in polyprots.columns: polyprots.drop('LCRs', inplace=True, axis=1)
            try:
                os.remove(cast_output_file)
                print('Old cast output found and deleted.')
            except:
                # File doesn't exist, that's okay
                pass

        try:
            cast_results = pd.read_csv(cast_output_file, index_col=0)
            print("Found cast ouput, using it.")
        except IOError as e:
            print("No cast output found. Generating...\nHold tight this can take awhile if the dataset is big...")
            sys.stdout.flush()
            try:
                CastRunner(polyprots['Sequence'].values.tolist(), polyprots.index.values.tolist(), cast_output_file)
                cast_results = pd.read_csv(cast_output_file, index_col=0)
                print('CAST run complete. Loading file and processing...')
                sys.stdout.flush()
            except Exception as e:
                raise

        # Clean up the results; compress into a single line; remove extraneous characters
        tmp = cast_results.applymap(lambda x: str(x) if pd.notnull(x) else np.nan)
        lcrs = (tmp + ':' + tmp.shift(-1, axis=1) + '_' + tmp.shift(-2, axis=1) + '$'+ tmp.shift(-3, axis=1)
          + '@' + tmp.shift(-4, axis=1) + ';').iloc[:, ::5].fillna('').sum(1)
        # Set blanks to null
        lcrs[lcrs.str.len() == 0] = np.nan
        polyprots.insert(polyprots.columns.get_loc('Protein names')+1, 'LCRs', lcrs)

    else:
        print('CAST output already in spreadsheet. Use forceCAST=True to force regeneration.')

    if ('IDRs' not in polyprots.columns and runIUPred) or forceIUPred:
	    
        iupred_output_file = "%s_iupred.csv"%input_file.rsplit('.csv', 1)[0]
        
        if forceIUPred:
            if 'IDRs' in polyprots.columns: polyprots.drop('IDRs', inplace=True, axis=1)
            try:
                os.remove(iupred_output_file)
                print('Old iupred output found and deleted.')
            except:
                # File doesn't exist, that's okay
                pass

        try:
            iupred_results = pd.read_csv(iupred_output_file, index_col=0)
            print('Found iupred output, using it.')
        except:
            print('No iupred output found. Generating...\nHold tight this can take awhile if the dataset is big...')
            sys.stdout.flush()
            try:
                IUPredRunner(polyprots['Sequence'].values.tolist(), polyprots.index.values.tolist(), iupred_output_file)
                print('IUPred run complete. Loading file and running thresholding...')
                sys.stdout.flush()
                iupred_results = pd.read_csv(iupred_output_file, index_col=0)
            except Exception as e:
                raise

        # Minimum length of disordered region to analyze
        thresholds = [1, 5, 10, 30, 50, 100]
        results = {k:[] for k in thresholds}

        above_threshold = iupred_results >= 0.5
        starts = ~above_threshold.shift(1, axis=1).fillna(False)
        ends = ~above_threshold.shift(-1, axis=1).fillna(False)

        for i in xrange(len(iupred_results)):
            # The first AA in a region is preceded by an AA not above the threshold
            fst = iupred_results.columns[above_threshold.iloc[i, :] & starts.iloc[i, :]].astype(int)
            # The last AA in a region is followed by an AA not above the threshold
            lst = iupred_results.columns[above_threshold.iloc[i, :] & ends.iloc[i, :]].astype(int)

            for threshold in thresholds:
                # Save regions longer than the minimum length as given by the threshold
                pr = [(i+1, j+1) for i, j in zip(fst, lst) if j >= i + threshold-1]
                results[threshold].append(';'.join(['%s_%s'%reg for reg in pr]) if len(pr) > 0 else np.nan)

        disordered = pd.DataFrame(results, columns=thresholds, index=iupred_results.index)
        disordered.columns=['Disordered Threshold: %d'%threshold for threshold in thresholds]
        polyprots = polyprots.join(disordered)

    print('\nWriting results to file.')
    polyprots.to_csv(input_file)


def getRawOutput(seqs, tmpfile, command, func):
    """
    Returns output from a given subprocess command that is run iteratively
    on a given sequence list. It writes to a temporary file on the disk that
    should be specified. Also, a text processing function can be passed that
    takes the command line program's output + the original sequence and
    translates theme into a list of values. This is than cast into a Pandas
    Series so care should be taken when using this function.
    """
    results = []
    raw = []
    
    # TODO: parallelize this. Sequentially is too slow for >10000 runs
    f = open(tmpfile, 'w')
    for i in xrange(len(seqs)):
        
        seq = seqs[i]
        # Write the sequence to a temporary file
        f.seek(0)
        f.write('>seq\n%s'%(seq))
        f.truncate()
        
        # If we're missing the sequence data then just add a blank column
        if(pd.isnull(seq)):
            results.append(pd.Series([np.nan]))
            continue
    
        # Run command on the temporary file; Process text
        probs_list = func(command(), seq)
        results.append(pd.Series(probs_list if len(probs_list) > 0 else [np.nan]))
        
    f.close()
    return results
    