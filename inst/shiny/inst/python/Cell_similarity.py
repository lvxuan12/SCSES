import numpy as np
from sklearn.metrics import pairwise_distances
def cell_similarity(data_similar,parameter):
    similar_method = parameter['similar_method']
    rows, cols = data_similar.shape
    d = pairwise_distances(data_similar, metric=similar_method)
    a = np.exp(-(np.square(d)))
    a = a + a.T
    m = a/((a.sum(axis=1)).reshape(rows,1))
    k = np.full(rows, rows)
    return(m)
