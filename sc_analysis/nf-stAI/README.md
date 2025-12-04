# nf-stAI

Nextflow pipeline to run stAI for spatial transcriptomics data imputation



### stAI
Imputation of spatial expression data from scRNA-seq data with [stAI](https://doi.org/10.1093/nar/gkaf158)

Environment:
Clone https://github.com/gszou99/stAI into 'bin'
```bash
conda activate 'nf-stAI'
python -m ipykernel install --user --name "nf-stAI" --display-name "nf-stAI"
conda env export | grep -v "^prefix: " > 'envs/nf-stAI.yml'
```

Run pipeline:
```bash
nextflow    -config conf/params.config \
            -log 'logs/nf-stai.log' \
            run main.nf  \
            -output-dir results_test
```

```bash
nextflow    -config conf/params.config \
            -log '../crc-cubes-spatial/data/30_imputed/stai/nf-stai.log' \
            run main.nf  \
            --adata_sp '../crc-cubes-spatial/data/30_imputed/sp_adata.h5ad' \
            --adata_sc '../crc-cubes-spatial/data/30_imputed/sc_refdata.h5ad' \
            --sc_layer 'counts' \
            --sp_layer 'counts' \
            --batch_key 'pt' \
            --spatial_knn 5 \
            --rna_batchsize 2048 \
            --spatial_batchsize 2048 \
            --d_hidden 64 \
            --d_latent 32 \
            --lam_clf 0.01 \
            --lam_cos 1 \
            --lam_genegraph 0 \
            --lam_impute 3 \
            --lam_mmd 0.1 \
            --lam_recon 1.5 \
            --eval_step 10 \
            --k_folds 5 \
            --lr 0.0005 \
            --num_epoch 5000 \
            --spatial_knn 8 \
            --test_features 'EPCAM' \
            --sc_label_key 'celltype' \
            --sp_label_key 'celltype' \
            -output-dir ../crc-cubes-spatial/data/30_imputed/stai
```



