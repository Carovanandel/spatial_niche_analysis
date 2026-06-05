import numpy as np
from spatialdata import read_zarr
import anndata as ad
import scanpy as sc
import squidpy as sq
import scvi
from pathlib import Path
import matplotlib.pyplot as plt
import cellcharter as cc
import os
import pandas as pd
import itertools
from sklearn.metrics import adjusted_rand_score, fowlkes_mallows_score
from datetime import datetime
import cellcharter as cc

import logging
logger = logging.getLogger('pytorch_lightning.utilities.rank_zero')
logger.setLevel(logging.ERROR)

from lightning.pytorch import seed_everything
seed_everything(12345)


def clustering_stability_GMM(n_runs, adata, use_rep, n_clusters):

    base_random_state = 12345

    def pairwise_stability(clusterings):
        """Compute pairwise ARI and FMI across clusterings."""
        ari_scores = []
        fmi_scores = []

        for c1, c2 in itertools.combinations(clusterings, 2):
            ari_scores.append(adjusted_rand_score(c1, c2))
            fmi_scores.append(fowlkes_mallows_score(c1, c2))

        return np.array(ari_scores), np.array(fmi_scores)


    results = []

    for n in n_clusters:
        clusterings = []
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f'[{now}] Starting {n} clusters GMM runs',flush=True)

        for run in range(n_runs):
            print(f'Computing GMM clustering run {run+1}/{n_runs} ',flush=True)

            gmm = cc.tl.Cluster(
                n_clusters=n,
                random_state=base_random_state + run
            )
            gmm.fit(adata, use_rep=use_rep)
            labels = gmm.predict(adata, use_rep=use_rep)
            clusterings.append(labels.astype(int))

        ari, fmi = pairwise_stability(clusterings)

        ari_mean = ari.mean()
        ari_sd   = ari.std()
        fmi_mean = fmi.mean()
        fmi_sd   = fmi.std()

        results.append({
            "n_clusters": n,
            "ari_mean": ari_mean,
            "ari_sd": ari_sd,
            "ari_ci_min": ari_mean - 1.96 * ari_sd,
            "ari_ci_max": ari_mean + 1.96 * ari_sd,
            "fmi_mean": fmi_mean,
            "fmi_sd": fmi_sd,
            "fmi_ci_min": fmi_mean - 1.96 * fmi_sd,
            "fmi_ci_max": fmi_mean + 1.96 * fmi_sd,
        })

    df_stability = pd.DataFrame(results)

    return df_stability

def clustering_stability_leiden(n_runs, adata, use_rep, resolutions):

    def pairwise_stability(clusterings):
        """Compute pairwise ARI and FMI across clusterings."""
        ari_scores = []
        fmi_scores = []

        for c1, c2 in itertools.combinations(clusterings, 2):
            ari_scores.append(adjusted_rand_score(c1, c2))
            fmi_scores.append(fowlkes_mallows_score(c1, c2))

        return np.array(ari_scores), np.array(fmi_scores)


    results = []

    # Dictionary to store clusterings and number of clusters per resolution
    clusterings_dict = {res: [] for res in resolutions}
    n_clusters_dict  = {res: [] for res in resolutions}

    for run in range(n_runs):
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f'[{now}] Starting Leiden run {run+1}/{n_runs}',flush=True)
        # Compute neighbors once per run
        sc.pp.neighbors(
            adata,
            use_rep=use_rep,
            n_neighbors=15,
            random_state= 1234 + run
        )

        # Compute Leiden clustering for all resolutions
        for res in resolutions:
            print(f'Computing Leiden clustering at resolution {res}',flush=True)
            key = f"spatial_clusters_{res}"
            sc.tl.leiden(
                adata,
                resolution=res,
                key_added=key,
                random_state= 1234 + run
            )
            labels = adata.obs[key].to_numpy()
            clusterings_dict[res].append(labels.astype(int))
            n_clusters_dict[res].append(len(np.unique(labels)))

    # Compute pairwise stability metrics for each resolution
    results = []
    for res in resolutions:
        ari, fmi = pairwise_stability(clusterings_dict[res])
        n_clusters = np.array(n_clusters_dict[res])

        results.append({
            "resolution": res,
            "ari_mean": ari.mean(),
            "ari_sd": ari.std(),
            "ari_ci_min": ari.mean() - 1.96 * ari.std(),
            "ari_ci_max": ari.mean() + 1.96 * ari.std(),
            "fmi_mean": fmi.mean(),
            "fmi_sd": fmi.std(),
            "fmi_ci_min": fmi.mean() - 1.96 * fmi.std(),
            "fmi_ci_max": fmi.mean() + 1.96 * fmi.std(),
            "n_clusters_mean": n_clusters.mean(),
            "n_clusters_sd": n_clusters.std(),
            "n_clusters_ci_min": n_clusters.mean() - 1.96 * n_clusters.std(),
            "n_clusters_ci_max": n_clusters.mean() + 1.96 * n_clusters.std(),
        })

    df_stability = pd.DataFrame(results)
    return df_stability


n_clusters_GMM = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
resolutions_Leiden = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
n_runs = 30


##### 480 GENES DATASET #####

adata_fibnuc_480 = sc.read_h5ad('cellcharter_clusters_anndata/480_dataset/Dimensionality_reduced_not_clustered/adata_480_fibnuc_trvae.h5ad')

# Compute neighborhood graph and aggregate neighbors
sq.gr.spatial_neighbors(adata_fibnuc_480, library_key='sample_id_new', coord_type='generic', radius=[0,400], delaunay=True)
cc.gr.aggregate_neighbors(adata_fibnuc_480, n_layers=3, use_rep='X_trVAE', out_key='X_cc_trVAE', sample_key='sample_id_new')

# Run GMM clustering on trVAE CellCharter embedding
print('starting 480 fibnuc stability run - GMM',flush=True)
df_stability = clustering_stability_GMM(n_runs=n_runs, adata=adata_fibnuc_480, use_rep='X_cc_trVAE', n_clusters=n_clusters_GMM)
df_stability.to_csv('clustering_stability_480_fibnuc_trVAE_GMM.csv')
print('saved clustering_stability_480_fibnuc_trVAE_GMM.csv',flush=True)

# Run Leiden clustering on trVAE Cellcharter embedding
print('starting 480 fibnuc stability run - Leiden clustering',flush=True)
df_stability = clustering_stability_leiden(n_runs=n_runs, adata=adata_fibnuc_480, use_rep='X_cc_trVAE', resolutions=resolutions_Leiden)
df_stability.to_csv('clustering_stability_480_fibnuc_trVAE_Leiden.csv')
print('saved clustering_stability_480_fibnuc_trVAE_Leiden.csv',flush=True)

##### 5000 GENES DATASET #####

adata_fibnuc_5000 = sc.read_h5ad('cellcharter_clusters_anndata/5000_dataset/Dimensionality_reduced_not_clustered/adata_5000_fibnuc_trvae.h5ad')

# Compute neighborhood graph and aggregate neighbors
sq.gr.spatial_neighbors(adata_fibnuc_5000, library_key='sample_id_new', coord_type='generic', radius=[0,400], delaunay=True)
cc.gr.aggregate_neighbors(adata_fibnuc_5000, n_layers=3, use_rep='X_trVAE', out_key='X_cc_trVAE', sample_key='sample_id_new')

# Run GMM clustering on trVAE CellCharter embedding
print('starting 5000 fibnuc stability run - GMM',flush=True)
df_stability = clustering_stability_GMM(n_runs=n_runs, adata=adata_fibnuc_5000, use_rep='X_cc_trVAE', n_clusters=n_clusters_GMM)
df_stability.to_csv('clustering_stability_5000_fibnuc_trVAE_GMM.csv')
print('saved clustering_stability_5000_fibnuc_trVAE_GMM.csv',flush=True)

# Run Leiden clustering on trVAE Cellcharter embedding
print('starting 5000 fibnuc stability run - Leiden clustering',flush=True)
df_stability = clustering_stability_leiden(n_runs=n_runs, adata=adata_fibnuc_5000, use_rep='X_cc_trVAE', resolutions=resolutions_Leiden)
df_stability.to_csv('clustering_stability_5000_fibnuc_trVAE_Leiden.csv')
print('saved clustering_stability_5000_fibnuc_trVAE_Leiden.csv',flush=True)