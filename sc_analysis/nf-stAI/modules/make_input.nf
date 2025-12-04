#!/usr/bin/env nextflow

process make_input {
    clusterOptions { "--job-name=nf-stAI-input-[${task.hash}]" }

    input:
        path(spdata_path)
        path(scdata_path)

    output:
        path("batches/*"), emit: batches
    
    script:
	"""
#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import sys
import anndata as ad
import scanpy as sc
from scipy import sparse
import yaml

def parse_null(x):
    return None if x == 'null' else x

batch_key = parse_null('${params.batch_key}')
sc_layer = parse_null('${params.sc_layer}')
sp_layer = parse_null('${params.sp_layer}')
sc_label_key = parse_null('${params.sc_label_key}')
sp_label_key = parse_null('${params.sp_label_key}')

params = {'data_paths': {'rna_data': None,
                         'spatial_data': None,
                         'save_path': None},
          'model_parameters': {'d_hidden': ${params.d_hidden},
                               'd_latent': ${params.d_latent},
                               'lam_clf': ${params.lam_clf},
                               'lam_cos': ${params.lam_cos},
                               'lam_genegraph':${params.lam_genegraph},
                               'lam_impute': ${params.lam_impute},
                               'lam_mmd': ${params.lam_mmd},
                               'lam_recon': ${params.lam_recon}},
          'training_parameters': {'device_id': ${params.device_id},
                                  'eval_step': ${params.eval_step},
                                  'k_folds': ${params.k_folds},
                                  'lr': ${params.lr},
                                  'num_epoch': ${params.num_epoch},
                                  'rna_batchsize': ${params.rna_batchsize},
                                  'spatial_batchsize': ${params.spatial_batchsize},
                                  'spatial_knn': ${params.spatial_knn}}}

spdata = sc.read_h5ad("${spdata_path}")
scdata = sc.read_h5ad("${scdata_path}")

if sp_layer is not None:
    spdata.X = spdata.layers[sp_layer].copy()

if sc_layer is not None:
    scdata.X = scdata.layers[sc_layer].copy()


if batch_key is None:
    batch_key = 'single_batch_key'
    scdata.obs[batch_key] = '1'
    spdata.obs[batch_key] = '1'

if sc_label_key is not None:
    scdata.obs['celltype'] = scdata.obs[sc_label_key].copy()
if sp_label_key is not None:
    spdata.obs['celltype'] = spdata.obs[sp_label_key].copy()
else:
    del spdata.obs['celltype']


# keep only training features
test_features = "${params.test_features}"
test_features = test_features.split(',')
keep = ~spdata.var_names.isin(test_features)
if sum(~keep) > 0:
    spdata = spdata[:,keep].copy()


# write data
datadir = ''
batches = list(set(scdata.obs[batch_key]).intersection(spdata.obs[batch_key]))
for b in batches:
    params_batch = params.copy()
    batchdir = os.path.join(datadir, 'batches', b)
    os.makedirs(batchdir, exist_ok=True)
    os.makedirs(os.path.join(batchdir, 'model'), exist_ok=True)
    scdata_batch = scdata[scdata.obs[batch_key] == b].copy()
    spdata_batch = spdata[spdata.obs[batch_key] == b].copy()
    scdata_batch.write_h5ad(os.path.join(batchdir, 'scdata.h5ad'))
    spdata_batch.write_h5ad(os.path.join(batchdir, 'spdata.h5ad'))
    params_batch['data_paths']['rna_data'] = os.path.join(b, 'scdata.h5ad')
    params_batch['data_paths']['spatial_data'] = os.path.join(b, 'spdata.h5ad')
    params_batch['data_paths']['save_path'] = os.path.join(b, 'model')
    with open(os.path.join(batchdir, 'model', 'train_config.yaml'), 'w') as ymlfile:
        yaml.dump(params_batch, ymlfile, default_flow_style=False)
	"""

}









