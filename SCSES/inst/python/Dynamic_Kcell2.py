import numpy as np
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import KNeighborsTransformer
from scipy import sparse
from scipy import io
import pandas as pd

def cell_similarity(data_similar,parameter):
    similar_method = parameter['similar_method']
    kmax = int(parameter['kmax'])
    kmin = int(parameter['kmin'])
    decay = int(parameter['decay'])
    alpha = int(parameter['alpha'])
    # d = pairwise_distances(data_similar, metric=similar_method)
    print(kmax)
    transformer = KNeighborsTransformer(
            metric=similar_method,
            n_jobs=-1,
            n_neighbors=kmax,
            mode='distance'
        )
    d = transformer.fit_transform(data_similar)
    block_size = 1000  # 每个块的大小
    rows, cols = data_similar.shape
    num_blocks = (rows // block_size) + 1  # 计算需要划分的块数

    a_result = np.zeros((rows, rows))  # 存储最终结果
    k_result = np.zeros(rows)  # 存储最终结果

    for i in range(num_blocks):
        start_row = i * block_size
        end_row = min((i + 1) * block_size, rows)
        # 划分矩阵为当前块
        current_block = d[start_row:end_row, :]
        current_block = current_block.toarray()
        rows_block, cols_block = current_block.shape
        d2 = np.sort(current_block, axis=1)
        d2 = d2[:,cols_block-kmax-1:cols_block]
        top = d2[:, 1:] - d2[:, :-1]
        top = top[:,1:]
        triup = np.triu(np.ones(top.shape[1]))
        dif_ave = top @ triup / np.arange(1, top.shape[1]+1)
        top = top[:,1:]
        dif_ave = dif_ave[:,:-1]
        diff_cutoff = top-dif_ave
        idx = pd.DataFrame(np.argwhere(diff_cutoff>0),columns=['row', 'col'])
        values = pd.DataFrame(diff_cutoff[diff_cutoff > 0],columns=['value'])
        idx = pd.concat([idx, values], axis=1)
        idx['mean_value'] = idx.groupby('row')['value'].transform('mean')
        idx = idx[idx['value'] > idx['mean_value']]
        idx_select = idx.groupby('row')['col'].min().reset_index().rename(columns={'index': 'row'}).to_numpy()
        idx_select[:,1]=idx_select[:,1]+2
        k=np.full(rows_block,kmax)
        k[idx_select[:,0]]=idx_select[:,1]
        k[np.where(k<kmin)]=kmin
        k[np.where(k>kmax)]=kmax
        sigma = 1/(d2[np.arange(rows_block),np.ceil(k/3).astype(int)])
        a = np.exp(-(np.square(current_block * sigma.reshape(rows_block,1))))
        a[np.where(current_block==0)]=0
        a[np.arange(rows_block),np.arange(start_row,end_row)]=1
        a_result[start_row:end_row, :] = a
        k_result[start_row:end_row] = k

    a_result = a_result + a_result.T
    k_result=k_result.astype(int)
    a_kth = np.sort(a_result,axis=1)[np.arange(0,rows),rows-k_result-1]
    a_result[np.where(a_result<a_kth.reshape(rows,1))] = 0
    m = a_result/((a_result.sum(axis=1)).reshape(rows,1))
    del a_result
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
    m1_kth = np.sort(m1,axis=1)[np.arange(0,rows),rows-k_result-1]
    m1[np.where(m1<m1_kth.reshape(rows,1))] = 0
    m2 = m1/((m1.sum(axis=1)).reshape(rows,1))
    ##save cs.mtx.gz
    m2_sparse_mtx = sparse.coo_matrix(m2)
    # io.mmwrite(work_path,m2_sparse_mtx)
    return(m2_sparse_mtx,k_result)
