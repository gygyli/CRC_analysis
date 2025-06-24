Scripts for running the analysis include:
1) Training-test set splitting
2) Alpha diversity analysis
3) Beta diversity analysis
4) MaAsLin2-based differential feature analysis
5) Random forest — within-cohort validation
6) Random forest — cross-cohort validation
7) GTDB_processing_scripts: Run the data_preprocessing script. The final output will be the `db/` folder. The path to the custom reference database is `db/library/genomes_info.txt`, which corresponds to GTDB:206273.  

```bash
bash /scripts/data_preprocessing -r gtdb_genomes_all/GTDB_complete/files --gm bac120_metadata.tsv --cluster

```
The genome files under the path gtdb_genomes_all/GTDB_complete/files correspond to all the genome entries listed in genomes_info_origin.txt, which is available via Zenodo at https://doi.org/10.5281/zenodo.15704298.
