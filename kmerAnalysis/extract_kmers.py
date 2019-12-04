import pandas as pd
import itertools
import math
import random
import os
from ctypes import cdll
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

    #load c functions
    basedir = os.path.abspath(os.path.dirname(__file__))
    libpath = os.path.join(basedir, 'combinekmerdata.so')
    combinekmers = cdll.LoadLibrary(libpath)
    
    #strings have to be bytes to be passed to c 
    kmrs_file_bytes = kmrs_file.encode('utf-8')
    tmp_file = kmrs_file + '.out'
    tmp_file_bytes = tmp_file.encode('utf-8')
    k_bytes = str(k).encode('utf-8')
    combinekmers.ck(kmrs_file_bytes, tmp_file_bytes, k_bytes)

    kmer_df = pd.read_csv(tmp_file, sep = '\t', names = ['Kmer', 'Count'])
    rmtree(out_dir) #delete temp directory

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