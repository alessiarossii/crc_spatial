#!/usr/bin/env nextflow

process run_stai {
    tag "${batchdir.baseName}"
    clusterOptions { "--job-name=nf-stAI-impute-[${task.tag}]" }

    input:
        path(batchdir)

    output:
        path("batches/*")
    
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
import torch
import sys

sys.path.insert(0, "$baseDir")
import bin.stAI.stAI as stAI

torch.cuda.empty_cache()

batchdir_str = "${batchdir}"

params = yaml.safe_load(open(os.path.join(batchdir_str, "model", "train_config.yaml"), "r"))
scdata_batch = sc.read_h5ad(os.path.join(batchdir_str, "scdata.h5ad"))
spdata_batch = sc.read_h5ad(os.path.join(batchdir_str, "spdata.h5ad"))

stAI.main.run_impute(model_parameters=params['model_parameters'],
                     training_parameters=params['training_parameters'],
                     adata_spatial=spdata_batch,
                     adata_rna=scdata_batch,
                     save_dir=params['data_paths']['save_path'])

genes_to_impute = list(scdata_batch.var.index)
imputed = stAI.main.impute_unmeasured_genes(genes_to_impute=genes_to_impute, model_dir=params['data_paths']['save_path'], agg=True, device='cpu')
spdata_batch_imp = ad.AnnData(
    X=None,
    obs=spdata_batch.obs.copy(),
    var=pd.DataFrame(data={'spatial': [gene in spdata_batch.var_names for gene in genes_to_impute]}, index=genes_to_impute)
)
spdata_batch_imp.layers['imputed'] = sparse.csr_matrix(imputed)
outdir = os.path.join("batches", os.path.basename(batchdir_str))
os.makedirs(outdir, exist_ok=True)
spdata_batch_imp.write_h5ad(os.path.join(outdir, "spdata_imputed.h5ad"))

	"""

}









