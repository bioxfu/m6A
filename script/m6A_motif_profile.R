library(RColorBrewer)
cols <- brewer.pal(3, 'Set1')
options(stringsAsFactors = F)

argv <- commandArgs(T)
start_file <- argv[1]
stop_file <- argv[2]
gene_list <- argv[3]
output <- argv[4]

# start_mat <- read.table('data/mouse_start_coden_region_GAC[AT]_matrix.tsv', sep='\t')
# stop_mat <- read.table('data/mouse_stop_coden_region_GAC[AT]_matrix.tsv', sep='\t')

start_mat <- read.table(start_file, sep='\t')
stop_mat <- read.table(stop_file, sep='\t')

set.seed(1111)
start_random <- start_mat[sample(nrow(start_mat), 10000), -c(1,2)]

set.seed(2222)
stop_random <- stop_mat[sample(nrow(stop_mat), 10000), -c(1,2)]

dfm <- read.table(gene_list)
dfm$V1 <- sub('\\..+', '', dfm$V1)
up <- dfm$V1[dfm$V2 > 0]
dn <- dfm$V1[dfm$V2 < 0]

up_stop <- stop_mat[stop_mat$V1 %in% up, -c(1,2)]
dn_stop <- stop_mat[stop_mat$V1 %in% dn, -c(1,2)]

up_start <- start_mat[start_mat$V1 %in% up, -c(1,2)]
dn_start <- start_mat[start_mat$V1 %in% dn, -c(1,2)]

pdf(output, hei=4, wid=8)
par(mfrow=c(1,2))
plot(smooth.spline(colMeans(up_start)*1000), col=cols[1], xlab='', ylab='Average Occurence of Motif', xaxt='n')
points(smooth.spline(colMeans(dn_start)*1000), col=cols[2])
points(smooth.spline(colMeans(start_random)*1000), col=cols[3])
legend('bottomright', c('up', 'down', 'random'), fill=cols, bty='n', border=NA)
axis(1, at=seq(0, 800, by=200), lab=c(-400,-200,'start',200,400))
abline(v=400, lty=2, col='gray')

plot(smooth.spline(colMeans(up_stop)*1000), col=cols[1], xlab='', ylab='Average Occurence of Motif', xaxt='n')
points(smooth.spline(colMeans(dn_stop)*1000), col=cols[2])
points(smooth.spline(colMeans(stop_random)*1000), col=cols[3])
axis(1, at=seq(0, 800, by=200), lab=c(-400,-200,'stop',200,400))
abline(v=400, lty=2, col='gray')
dev.off()


