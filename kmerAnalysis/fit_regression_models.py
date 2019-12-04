import pandas as pd 
import numpy as np
from bloom_filter import BloomFilter
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Lasso
import matplotlib.pyplot as plt

def summarise_kmer_data(kmer_dfs, y_dict):
    '''Summarises multiple kmer dataframes in format outputted from make_kmer_df() 
    so regression models can be fitted

    Parameters:
    kmer_dfs (dict): dictionary of kmer dataframes in format defined by make_kmer_df(), 
    key is name of dataframe with corresponding key in y_dict
    y_dict (dict): dictionary of y data as values and keys corresponding to those in kmer_dfs

    Returns:
    counts_array, y_array, represnted_kmers: three arrays of the counts of each kmer and 
    y values in same order and the kmers represented in the data
    '''

    max_kmers = len(kmer_dfs[list(kmer_dfs.keys())[0]]) #below code checks all are the same length

    represented_kmers = []
    bloom = BloomFilter(max_elements=max_kmers, error_rate=0.1)
    for name in kmer_dfs:
        kmer_df = kmer_dfs[name]
        if len(kmer_df) != max_kmers:
            raise ValueError('kmer_dfs are not of equal length')
        for kmer in list(kmer_df.loc[kmer_df['Count'] != 0,]['Kmer']): 
            if kmer not in bloom:
                bloom.add(kmer)
                represented_kmers.append(kmer)

    counts_array = []
    y_array = []
    for name in kmer_dfs:
        if name in y_dict:
            kmer_df = kmer_dfs[name]
            counts = list(kmer_df.loc[kmer_df['Kmer'].isin(represented_kmers),]['Count'])
            counts_array.append(counts)
            y_array.append(y_dict[name])

    return(counts_array, y_array, represented_kmers)


def fit_lasso_model(x, y, alpha = 0.1, plot = True, test_size = 0.3):
    '''Fits lasso regresson model 

    Parameters: 
    x (array): numpy array of values of x in same order as y
    y (array): numpy array of value of y in same order as x
    alpha (float): regularization parameter for lasso regression
    plot (bool): true or false plot model predictions against real data
    test_size (float): proportion of sample to be in test split

    Returns: 
    train_score, test_score, coeff_used, lasso, plot (if plot = True):
    respectively r squared of training data with model predictions, r squared of testing 
    data with model predictions, number of items in x whose coefficient is not 0, the fitted
    model, and a plot of the model predictions against real data
    '''

    X_train,X_test,y_train,y_test = train_test_split(x, y, test_size=test_size, 
                                                    random_state=31)

    lasso = Lasso(alpha = alpha)
    lasso.fit(X_train,y_train)
    train_score = lasso.score(X_train,y_train) #R squared data
    test_score = lasso.score(X_test,y_test) 
    coeff_used = np.sum(lasso.coef_!=0) #number of points in x used by model

    if plot:
        plot = plt.scatter(y,lasso.predict(x))
        return(train_score, test_score, coeff_used, lasso, plot)
    else:
        return(train_score, test_score, coeff_used, lasso)