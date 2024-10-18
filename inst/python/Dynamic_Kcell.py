import numpy as np
from sklearn.metrics import pairwise_distances
import pandas as pd

def cell_similarity(data_similar,parameter):
    similar_method = parameter['similar_method']
    kmax = int(parameter['kmax'])
    kmin = int(parameter['kmin'])
    decay = parameter['decay']
    alpha = parameter['alpha']
    rows, cols = data_similar.shape
    d = pairwise_distances(data_similar, metric=similar_method)
    d2 = np.sort(d, axis=1)[...,::1]
    top = d2[:, 1:] - d2[:, :-1]
    top = top[:,1:]
    triup = np.triu(np.ones(top.shape[1]))
    dif_ave = top @ triup / np.arange(1, top.shape[1]+1)
    top = top[:,1:]
    dif_ave = dif_ave[:,:-1]
    #idx = pd.DataFrame(np.argwhere(top-dif_ave>dif_ave[:,-1].reshape(dif_ave.shape[0],1)),columns=['row', 'col'])
    diff_cutoff = top-dif_ave
    idx = pd.DataFrame(np.argwhere(diff_cutoff>0),columns=['row', 'col'])
    values = pd.DataFrame(diff_cutoff[diff_cutoff > 0],columns=['value'])
    idx = pd.concat([idx, values], axis=1)
    idx['mean_value'] = idx.groupby('row')['value'].transform('mean')
    idx = idx[idx['value'] > idx['mean_value']]
    idx_select = idx.groupby('row')['col'].min().reset_index().rename(columns={'index': 'row'}).to_numpy()

    idx_select[:,1]=idx_select[:,1]+2
    idx_select[:,1][np.where(idx_select[:,1]>kmax)]=kmax
    idx_select[:,1][np.where(idx_select[:,1]<kmin)]=kmin
    #cells do not satisfy top>dif_ave
    idx_select_all=np.array([np.arange(0,rows),np.repeat(kmax,rows)]).T
    idx1=idx_select[:,0]
    idx2=np.repeat(1,idx_select.shape[0])
    idx_select_all[idx1,idx2]=idx_select[:,1]
    k=idx_select_all[:,1]
    thred=d2[idx_select_all[:,0],idx_select_all[:,1]]
    sigma = 1/(np.sort(d,axis=1)[idx_select_all[:,0],np.ceil(idx_select_all[:,1]/3).astype(int)])
    a = np.exp(-(np.square(d * sigma.reshape(rows,1))))
    a = a + a.T
    a_kth = np.sort(a,axis=1)[np.arange(0,rows),rows-k-1]
    a[np.where(a<a_kth.reshape(rows,1))] = 0
    m = a/((a.sum(axis=1)).reshape(rows,1))
    m1 = alpha * (np.dot(m,m)) + (1-alpha) * m
    sse = np.sum(np.square(m - m1))
    ssr = np.sum(np.square(m1 - np.mean(m,axis=None)))
    delta = sse/(sse+ssr)
    delta_all = np.array([delta])
    while delta > decay:
        m2 = alpha * np.dot(m1,m) + (1-alpha) * m
        sse = np.sum(np.square(m2 - m1))
        ssr = np.sum(np.square(m2 - np.mean(m1,axis=None)))
        delta = sse/(sse+ssr)
        m1 = m2
        delta_all = np.append(delta_all,delta)
    m1_kth = np.sort(m1,axis=1)[np.arange(0,rows),rows-k-1]
    m1[np.where(m1<m1_kth.reshape(rows,1))] = 0
    m2 = m1/((m1.sum(axis=1)).reshape(rows,1))
    print("Delta:",delta_all)
    return(m2,k)
