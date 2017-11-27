library(biomaRt)

argv <- commandArgs(T)
species <- argv[1]
# species <- 'human'

if (species == 'human') {
  ensembl <- useMart('ensembl', dataset = "hsapiens_gene_ensembl")
}else if (species == 'mouse') {
  ensembl <- useMart('ensembl', dataset = "mmusculus_gene_ensembl")
}

# filters <- listFilters(ensembl)
# attributes <- listAttributes(ensembl)

bm <- getBM(attributes=c('ensembl_gene_id', 
                         'ensembl_transcript_id',
                         'transcript_length',
                         'rank', 
                         'cdna_coding_start', 
                         'cdna_coding_end'), 
            filters = 'transcript_biotype', 
            values = 'protein_coding', 
            mart = ensembl)

write.table(bm, paste0('data/', species, '_biomart.tsv'), sep='\t', row.names = F, quote = F)
