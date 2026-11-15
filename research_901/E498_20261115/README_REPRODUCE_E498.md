# E498 original-source reproducibility and quantum/optical status

## Authentic original source

Official published study Gervais et al. 2026 *Molecular Systems Biology*, DOI `10.1038/s44320-026-00226-5`, Zenodo version `18059846` with full ZIP of 842127932 bytes. Original full-ZIP MD5 from provider is `207baa767ed43b767f2c97d322866e2d`, **NOT independently verified because full archive was not downloaded**. Four exact authentic member Feather files were directly acquired using strict HTTP 206 Range, central-directory local-header location checks, CRC32, byte length and SHA256. Cloud workflow in `dairongzhen3-creator/math`, branch `cloud-901-e498-starvation-degradation-20261115`, `.github/workflows/901-e498-original-degradation.yml`; receipt `E498_original_source_receipts.json`.

Downloaded exact original source SHA256: 

- `mydegrates.feather`: `2a36a01facd1bf79ab0568e88cf07ff6cfebf099ecf9419ba85a0ead19cf7676`.
- `myframes_degradation.feather`: `e005351f269c85804c57f5f28d893f1db3f353ed139da4069a87820d491a4131`.
- `myframes_bleaching.feather`: `6e1733d6b641f5a47dca2851520fb9a0edfed2aad6faecbaf1ec5eece99917f7`.
- `myproteomicsdata_multistrain.feather`: `d0e83674a86d248f29cf7404a86251faca98b76cb4102e528ab3c4d48cc37386`.

Author original complete `Gervais_et_al_analysis.Rmd` SHA `7f915a34db6c3326a0be833bc1c04e57127d90684ce281c9fe4b0cd65a48de70`. README SHA `20ed23f5cdf995ed22536ad487ed8df6f1a46ab51d954098e400efdb27c318af`. Preserve both. For a feather reader in source conversion, Github Actions installed pyarrow and pandas and exported **ALL rows** to equivalent gzip CSV in original column order (`%.17g` float), producing `E498_original_to_csv_gzip_receipt.json` with SHA for 4 outputs; no source filters, no sort and no inferred values during conversion. Both source artifacts have separately archived original GitHub workflow artifact ZIPs and immutable split/rehydration manifests in the Google Drive E498 warehouse.

## Main and independent reproducibility

Place `analyze_e498.py`, `independent_qa_e498.py`, the four exact `*.feather.csv.gz` exported files, `E498_original_to_csv_gzip_receipt.json`, and four exact original `*.feather` bytes in the same directory; execute:

```
python analyze_e498.py
python independent_qa_e498.py
```

The main script needs numpy/pandas. Independent QA needs only Python standard library (CSV/gzip/hashlib/statistics/Decimal) and DOES NOT import numpy/pandas or main functions. Both validate original four Feather byte hashes and all four gzip exported byte SHA; 1,338,282 source rows; source original date, cell, protein, color channels; 403 paired source cells/368 early loss greater; 1659 bleaching WLS cell-phases; source proteomics 16543 missing; and all author-four-pigment model survival curves independently.

Only source author-expressed fluor channels for degradation curves; nonexpressed spectral channels are not biological replicates. Source actual 2 dates degradation and 3 dates bleaching. No technical independent camera calibration and no absolute molecule copies. Neither physical transmission MI, information rate nor bit/J is identified. Do not add previous E488/E489/E494 source optical predictive scores.

Research paper `PUBLICATION_NO_GO` until the required whole-article citation review and PaperQA plus physical identifiability. Source 2018 bacterial fluorescence-to-protein conversion from another instrument must not be applied here.
