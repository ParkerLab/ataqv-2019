# Processing of public ATAC-seq data
# TODO: finish this
1. First, the data needs to be downloaded (`make fastq`)
2. Next, launch atacseq primary processing (`make atacseq`)
3. Launch ataqv (`make atacseq`)
4. Create ataqv sessions for each project separately, and for the bulk and single-cell datasets separately (`make ataqv_sessions`)
5. Heatmap data can be prepared using `make heatmaps`
