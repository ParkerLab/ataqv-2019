# Running
1. First, trim read lengths to 36 bps (makes comparable with the other experiment): `make trim-read-length`
2. Carry out the primary ATAC-seq processing (adapter trimming, fastqc, mapping, duplicate marking, filtering, peak calling, ataqv): `make atacseq`. This will submit jobs to the queue, and you need to wait for these to finish before proceeding.
3. You can have a look at the generic ataqv session: `make ataqv_session`
4. Next, make the master peaks: `make master_peaks`
5. Then get the read counts in the master peaks: `make master_peak_counts`
6. In order to make certain statistics comparable across libraries (e.g., % of reads in peaks), we want to use ataqv metrics created with the master peaks. Therefore, re-run ataqv using master peaks: `make ataqv_master_peaks`
7. ...and then create the corresponding session: `make ataqv_session_master_peaks`
8. Run PCA: `make PCA`
9. Run the NB-GLM: `make tn5_sensitive_peaks`
10. Do the chromHMM overlap analysis: `make read_chromhmm_overlaps`
11. Do the TF ChIP-seq overlap analysis: `make read_tf_overlaps`
