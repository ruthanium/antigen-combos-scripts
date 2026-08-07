# **********************************
# Build sketches per cancer type.
# **********************************
import numpy as np
import pandas as pd
from fbpca import pca
from geosketch import gs
import os

# sketch size
sk_sz = 0.2

# more robust file paths
currdir = os.path.dirname(__file__)
if "src" in currdir:
    currdir = currdir + "/"

# read in input matrix
matAll = pd.read_csv(currdir + "../data/test-normalized-matrix.txt")
cans = sorted(set(matAll.loc[matAll["type"] == "cancer", "tissue.cancer"].tolist()))

# make relevant directory structure
if not os.path.exists(currdir + "../results"):
    os.makedirs(currdir + "../results")

if not os.path.exists(currdir + "../results/sketches"):
    os.makedirs(currdir + "../results/sketches")

# iterate through the cancer types to make sketches
for can in cans:
    print(can)

    #subset and label
    mat = matAll.loc[(matAll["tissue.cancer"] == can),]
    mat = pd.concat([mat, matAll.loc[matAll["type"] == "tissue"]], ignore_index=True)
    N_samples = len(mat["tissue.cancer"])

    # remove unnecessary fields
    matt = mat.copy()
    matt.drop("batch", axis=1, inplace=True)
    matt.drop("type", axis=1, inplace=True)
    matt.drop("tissue.cancer", axis=1, inplace=True)
    matt.set_index('dataset', inplace=True)
    mattm = matt.values

    # compute the PCs - necessary input for the sketching
    n_components = min(100, min(mattm.shape))
    U, s, Vt = pca(mattm, k=n_components)
    X_dimred = U[:, :n_components] * s[:n_components]

    # sketch
    N = int(N_samples * sk_sz) # Number of samples to obtain from the dataset
    sketch_index = gs(X_dimred, N, replace=False, seed=0)
    # Get the sample identifiers directly from the selected row indices.
    # This avoids ambiguous floating-point joins when PCA coordinates repeat.
    selected = list(matt.index[sketch_index])

    out = open(currdir + "../results/sketches/" + can.lower().replace(" ", "-") + "-sketch.txt", "w")
    for elm in selected:
        out.write(elm + "\n")
    out.close()
