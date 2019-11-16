# Running
1. First, carry out the primary ATAC-seq processing (adapter trimming, fastqc, mapping, duplicate marking, filtering, peak calling, ataqv): `make atacseq`
2. Do all the downstream processing: First, set the ROOT variable in downstream.config to point to your ATAQV_HOME; then, run `make downstream`
