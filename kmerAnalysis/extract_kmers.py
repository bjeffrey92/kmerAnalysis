import pandas as pd
import itertools
import math
import random
import os
from shutil import rmtree

def make_kmer_df(fasta, k, binary = False):
    '''Takes a fasta as input and return dataframe of all kmers in the sequence

    Parameters:
    fasta (str): Valid fasta file
    k (int): Length of kmers to be extracted
    binary (bool): If true output will only record whether each kmer is present rather than the number of times is appears

    Returns:
    kmer_df: dataframe with the kmer counts for the fasta file
    '''

    #make temp directory 
    counter = 0
    out_dir = str(random.getrandbits(64)) + str(counter)
    while os.path.isdir(out_dir):
        counter += 1
    out_dir = out_dir[:-1] + str(counter)
    os.makedirs(out_dir)

    #uses subprocesses to run kmc tool for kmer counting
    name = name = fasta.split('/')[-1]
    shell_command = 'kmc -k{k} -fa -ci1  {fasta} {out_dir}/{name} {out_dir}'.format(k = str(k),
                                                                                    fasta = fasta,
                                                                                    name = name,
                                                                                    out_dir = out_dir)
    os.system(shell_command)
    shell_command = 'kmc_dump {out_dir}/{name} {out_dir}/{name}.{k}.kmrs'.format(out_dir = out_dir,
                                                                            name = name,
                                                                            k = k)
    os.system(shell_command)

    kmrs_file = '{out_dir}/{name}.{k}.kmrs'.format(out_dir = out_dir,
                                                    name = name,
                                                    k = k)
    kmer_data = pd.read_csv(kmrs_file, sep = '\t', header = None)
    rmtree(out_dir) #delete temp directory

    #save results to pandas dataframe
    bases = ['A','T','C','G']
    all_kmers = []
    counts = []
    for i in itertools.product(bases, repeat=k):
        kmer = list(i)
        kmer = ''.join(kmer)
        all_kmers.append(kmer)
        if kmer in list(kmer_data[0]):
            if binary:
                count = 1
            else:
                count = list(kmer_data[kmer_data[0] == kmer][1])[0]
        else:
            count = 0
        counts.append(count)
    d = {'Kmer':all_kmers,'Count':counts}
    kmer_df = pd.DataFrame(d)
    
    return(kmer_df)


def calculate_optimal_k(seq_length, q):
    '''Uses equation from https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x
    to calculate the optimal length of kmer to use for a nucleotide sequence of given length 

    Parameters:
    seq_length (int): Number of mucleotide bases in a sequence
    q (int): Desired probability of observing a random kmer in sequence of length seq_length

    Returns:
    k: length of a kmer that will appear in sequence of length seq_length with probability q
    '''
    k = math.log(seq_length * (1 - q)/q, 4)

    return(k)