mkdir data

# download cDNA of Human (Homo sapiens) and Mouse (Mus musculus)
wget ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
mv *gz data

# extract cDNA of protein coding genes
python script/extract_protein_coding_cdna.py data/Homo_sapiens.GRCh38.cdna.all.fa.gz data/human_cdna_protein_coding.fa.gz
python script/extract_protein_coding_cdna.py data/Mus_musculus.GRCm38.cdna.all.fa.gz data/mouse_cdna_protein_coding.fa.gz

# download Biomart data
Rscript script/biomart.R human
Rscript script/biomart.R mouse

# transform the Biomart data
cat data/human_biomart.tsv|sort -k2,2 -k4,4n|groupBy -inheader -g 1,2,3 -c 4,5,6 -o collapse,collapse,collapse > data/human_biomart_groupby.tsv
cat data/mouse_biomart.tsv|sort -k2,2 -k4,4n|groupBy -inheader -g 1,2,3 -c 4,5,6 -o collapse,collapse,collapse > data/mouse_biomart_groupby.tsv

# extract start and stop codon region (5' and 3' flank 400 bp region)
python script/extract_start_stop_coden_region.py data/human_cdna_protein_coding.fa.gz  data/human_biomart_groupby.tsv  data/human_start_stop_coden_region.tsv
python script/extract_start_stop_coden_region.py data/mouse_cdna_protein_coding.fa.gz  data/mouse_biomart_groupby.tsv  data/mouse_start_stop_coden_region.tsv

# search motif
python script/scan_motif.py GAC[AT] 4 data/human_start_stop_coden_region.tsv data/human_start_coden_region_GAC[AT]_matrix.tsv data/human_stop_coden_region_GAC[AT]_matrix.tsv
python script/scan_motif.py GAC[AT] 4 data/mouse_start_stop_coden_region.tsv data/mouse_start_coden_region_GAC[AT]_matrix.tsv data/mouse_stop_coden_region_GAC[AT]_matrix.tsv

# search motif
python script/scan_motif.py [AGT][AG]AC[ACT] 5 data/human_start_stop_coden_region.tsv data/human_start_coden_region_DRACH_matrix.tsv data/human_stop_coden_region_DRACH_matrix.tsv
python script/scan_motif.py [AGT][AG]AC[ACT] 5 data/mouse_start_stop_coden_region.tsv data/mouse_start_coden_region_DRACH_matrix.tsv data/mouse_stop_coden_region_DRACH_matrix.tsv

# plot motif profile
Rscript script/m6A_motif_profile.R data/mouse_start_coden_region_GAC[AT]_matrix.tsv data/mouse_stop_coden_region_GAC[AT]_matrix.tsv table/MeCP2_KO_mice.tsv figure/MeCP2_KO_mice_GAC[AT].pdf
Rscript script/m6A_motif_profile.R data/mouse_start_coden_region_DRACH_matrix.tsv data/mouse_stop_coden_region_DRACH_matrix.tsv table/MeCP2_KO_mice.tsv figure/MeCP2_KO_mice_DRACH.pdf

Rscript script/m6A_motif_profile.R data/human_start_coden_region_GAC[AT]_matrix.tsv data/human_stop_coden_region_GAC[AT]_matrix.tsv table/MeCP2_293T.tsv figure/MeCP2_293T_GAC[AT].pdf
Rscript script/m6A_motif_profile.R data/human_start_coden_region_DRACH_matrix.tsv data/human_stop_coden_region_DRACH_matrix.tsv table/MeCP2_293T.tsv figure/MeCP2_293T_DRACH.pdf
