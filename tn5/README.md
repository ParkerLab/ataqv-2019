# Running
1. First, carry out the primary ATAC-seq processing (adapter trimming, fastqc, mapping, duplicate marking, filtering, peak calling, ataqv): `make atacseq`
2. You can have a look at the generic ataqv session: `make ataqv_session`
3. Next, make the master peaks: `make master_peaks`
4. Then get the read counts in the master peaks: `make master_peak_counts`
5. In order to make certain statistics comparable across libraries (e.g., % of reads in peaks), we want to use ataqv metrics created with the master peaks. Therefore, re-run ataqv using master peaks: `make ataqv_master_peaks`
6. ...and then create the corresponding session: `make ataqv_session_master_peaks`
7. Run PCA: `make pca`
8. Normalize the bigwig tracks: `make normalization`
9. Run the NB-GLM: `make tn5_sensitive_peaks`
10. Do the chromHMM overlap analysis: `make read_chromhmm_overlaps`
11. Do the TF ChIP-seq overlap analysis: `make read_tf_overlaps`
12. Make heatmap data: `make heatmaps`
13. Calculate TSS enrichment under different methods (`make tss_enrichment.ataqv; make tss_enrichment.ENCODE; make tss_enrichment.cutsite`)
