import pandas as pd
import itertools
import math


def make_kmer_df(fasta, k, binary = False):
    '''Takes a fasta as input and return dataframe of all kmers in the sequence

    Parameters:
    fasta (str): Valid fasta file
    k (int): Length of kmers to be extracted
    binary (bool): If true output will only record whether each kmer is present rather than the number of times is appears

    Returns:
    kmer_df: dataframe with the kmer counts for the fasta file
    '''
    #create dataframe for kmer counts to be stored
    bases = ['A','T','C','G']
    all_kmers = []
    counts = []
    for i in itertools.product(bases, repeat=k):
        kmer = list(i)
        kmer = ''.join(kmer)
        all_kmers.append(kmer)
        counts.append(0)
    d = {'Kmer':all_kmers,'Count':counts}
    kmer_df = pd.DataFrame(d)

    #extract the coding sequence as a string
    lines = []
    with open(fasta, 'r') as a:
        for line in a:
            if line[0] == '>':
                if 'contig' in locals() and len(contig) != 0: #if one contig has already been completed add it to list
                    contig = ''.join(contig)     
                    lines.append(contig)
                contig = []  
            elif line[-1] == '\n': #remove line breaks
                line = line[:-1] 
                contig.append(line)
            else:
                contig.append(line)
        if len(lines) == 0:
            contig = ''.join(contig)
            lines.append(contig)

    #fill in kmer dataframe 
    for line in lines:
        num_kmers = len(line) - k + 1
        for i in range(num_kmers):
            kmer = line[i:i+k] #Slice the string to get the kmer
            if binary:
                kmer_df.loc[kmer_df['Kmer'] == kmer, 'Count' ] = 1 #add to df
            else:
                kmer_df.loc[kmer_df['Kmer'] == kmer, 'Count' ] += 1 #add to df

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