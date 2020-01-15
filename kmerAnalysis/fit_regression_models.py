import pandas as pd 
import numpy as np
import pickle
import warnings
from bloom_filter import BloomFilter
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Lasso
import sklearn.svm
import matplotlib.pyplot as plt

def summarise_kmer_data(kmer_dfs, y_dict, names):
    '''Summarises multiple kmer dataframes in format outputted from make_kmer_df() 
    so regression models can be fitted

    Parameters:
    kmer_dfs (list): list of paths to pickled kmer dataframes in format defined by make_kmer_df()
    y_dict (dict): dictionary of y data as values and keys corresponding to those in kmer_dfs
    names (list): list of names of kmer_dfs in same order as kmer_dfs list

    Returns:
    counts_array, y_array, represnted_kmers: three arrays of the counts of each kmer and 
    y values in same order and the kmers represented in the data
    '''
    with open(kmer_dfs[0], 'rb') as a:
        kmer_df = pickle.load(a)
        max_kmers = len(kmer_df) #below code checks all are the same length
    bloom = BloomFilter(max_elements=max_kmers, error_rate=0.1)
    represented_kmers = []
    counts_array = [None]*len(names)
    y_array = [None]*len(names)
  
    for kmer_df_file in kmer_dfs: 
        with open(kmer_df_file, 'rb') as a:
            kmer_df = pickle.load(a)
        if len(kmer_df) != max_kmers:
            raise ValueError('kmer_dfs are not of equal length')
        for kmer in list(kmer_df.loc[kmer_df['Count'] != 0,]['Kmer']): 
            if kmer not in bloom:
                bloom.add(kmer)
                represented_kmers.append(kmer)
        
    pos = 0
    for kmer_df_file in kmer_dfs: 
        with open(kmer_df_file, 'rb') as a:
            kmer_df = pickle.load(a)
        name = names[pos]
        if name in y_dict:
            y_array[pos] = y_dict[name]
            counts = list(kmer_df.loc[kmer_df['Kmer'].isin(represented_kmers),]['Count'])
            counts_array[pos] = counts
        else:
            warnings.warn(name + ' missing from y_dict')
        pos += 1

    counts_array = [i for i in counts_array if i != None]
    y_array = [i for i in y_array if i != None]

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

def fit_svr_model(x, y, C = 0.1, epsilon = 0.1, plot = True, test_size = 0.3):
    '''Fits support vector regression model 

    Parameters: 
    x (array): numpy array of values of x in same order as y
    y (array): numpy array of value of y in same order as x
    c (float): regularization parameter for svr model
    epsilon (float): width of hyperplane in svr model
    plot (bool): true or false plot model predictions against real data
    test_size (float): proportion of sample to be in test split

    Returns: 
    train_score, test_score, lasso, plot (if plot = True):
    respectively r squared of training data with model predictions, r squared of testing 
    data with model predictions, number of items in x whose coefficient is not 0, the fitted
    model, and a plot of the model predictions against real data
    '''

    X_train,X_test,y_train,y_test = train_test_split(x, y, test_size=test_size, 
                                                    random_state=31)

    svr = sklearn.svm.SVR(gamma = 'scale', C = C, epsilon = epsilon)  
    svr.fit(X_train, y_train)
    train_score = svr.score(X_train, y_train)
    test_score = svr.score(X_test, y_test)

    if plot:
        plot = plt.scatter(y,svr.predict(x))
        return(train_score, test_score, svr, plot)
    else:
        return(train_score, test_score, svr)
