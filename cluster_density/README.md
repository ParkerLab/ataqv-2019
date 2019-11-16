# Running
1. First, trim read lengths to 36 bps (makes comparable with the other experiment): `make trim-read-length`
2. Carry out the primary ATAC-seq processing (adapter trimming, fastqc, mapping, duplicate marking, filtering, peak calling, ataqv): `make atacseq`. This will submit jobs to the queue, and you need to wait for these to finish before proceeding.
3. Do all the downstream processing: First, set the ROOT variable in downstream.config to point to your ATAQV_HOME; then, run `make downstream`
