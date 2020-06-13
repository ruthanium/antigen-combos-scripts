# **********************************
# Build sketches per cancer type.
# **********************************
import numpy as np
import pandas as pd
from fbpca import pca
from geosketch import gs

# sketch size
sk_sz = 0.2

# read in input matrix
matAll = pd.read_csv("../data/test-normalized-matrix.txt")
cans = set(matAll.loc[matAll["type"] == "cancer", "tissue.cancer"].tolist())

# iterate through the cancer types to make sketches
for can in cans:
    print(can)

    #subset and label
    mat = matAll.loc[(matAll["tissue.cancer"] == can),]
    mat = mat.append(matAll.loc[(matAll["type"] == "tissue"),])
    mat["cluster"] = 2
    mat.loc[mat['tissue.cancer'] == can, "cluster"] = 1
    N_samples = len(mat["tissue.cancer"])

    # remove unnecessary fields
    matt = mat.copy()
    matt.drop("batch", axis=1, inplace=True)
    matt.drop("type", axis=1, inplace=True)
    matt.drop("cluster", axis=1, inplace=True)
    matt.drop("tissue.cancer", axis=1, inplace=True)
    matt.set_index('dataset', inplace=True)
    mattm = matt.values

    # compute the PCs - necessary input for the sketching
    U, s, Vt = pca(mattm, k=100)
    X_dimred = U[:, :100] * s[:100]

    # sketch
    N = int(N_samples * sk_sz) # Number of samples to obtain from the dataset
    sketch_index = gs(X_dimred, N, replace=False)
    X_sketch = X_dimred[sketch_index]

    # get the samples selected in the sketch and output
    reduced = pd.DataFrame(X_sketch)
    pca_out = pd.DataFrame(X_dimred)
    pca_out["dataset"] = list(matt.index)
    red_with_labs = pd.merge(pca_out, reduced, how="inner", on=list(reduced.columns.values))
    selected = list(red_with_labs["dataset"])

    out = open("../results/sketches/" + can.lower().replace(" ", "-") + "-sketch.txt", "w")
    for elm in selected:
        out.write(elm + "\n")
    out.close()
