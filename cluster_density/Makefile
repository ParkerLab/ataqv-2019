include $(ATAQV_HOME)/Makefile
LOCAL_ROOT=$(ROOT)/cluster_density
SAMPLE_INFO=$(ROOT)/cluster_density/sample_info/sample_info.txt
FIGURE = $(LOCAL_ROOT)/figures

ALL_LIBRARIES=$(shell cut -f6 $(SAMPLE_INFO) | awk 'NR>1' | sort | uniq)
LIBRARIES=$(shell cut -f6 $(SAMPLE_INFO) | awk 'NR>1' | sort | uniq | grep -v 81997) # the good libraries, since one was excluded -- 0.5X Tn5, rep 3
LIBRARIES_1X=$(shell grep '1X' $(SAMPLE_INFO) | cut -f6 | sort | uniq | grep -v 81997)
HIGH_CD_LIBRARIES=$(shell grep 1833 $(SAMPLE_INFO) | cut -f6 | sort | uniq | grep -v 81997) # the good libraries, since one was excluded -- 0.5X Tn5, rep 3
NAME=$(CD_NAME)

READGROUPS = $(shell ls $(DATA)/fastq/run_1830/* $(DATA)/fastq/run_1833/* | perl -pe 's@.*run_\d+/(.*)\.[12].fastq.gz@$$1@g' | sort | uniq)
HIGH_CD_READGROUPS = $(shell ls $(DATA)/fastq/run_1833/* | perl -pe 's@.*run_\d+/(.*)\.[12].fastq.gz@$$1@g' | sort | uniq)
LOW_CD_READGROUPS = $(shell ls $(DATA)/fastq/run_1830/* | perl -pe 's@.*run_\d+/(.*)\.[12].fastq.gz@$$1@g' | sort | uniq)
LIBRARY2NAME=$(shell grep $(1) $(SAMPLE_INFO) | cut -f2) Tn5, $(shell grep $(1) $(SAMPLE_INFO) | cut -f3), $(shell grep $(1) $(SAMPLE_INFO) | cut -f7) cluster density

trim-read-length: TARGET_LENGTH = 36
trim-read-length:
	mkdir -p $(ANALYSIS)
	echo "# drmr:job working_directory=$(ANALYSIS)" > $(PIPELINE)
	$(foreach r,$(HIGH_CD_READGROUPS),echo "zcat $(DATA)/fastq/run_1833/$(r).1.fastq.gz | fastx_trimmer -l $(TARGET_LENGTH) -z -o $(ANALYSIS)/$(r).1.fastq.gz" >> $(PIPELINE)$(NL))
	$(foreach r,$(HIGH_CD_READGROUPS),echo "zcat $(DATA)/fastq/run_1833/$(r).2.fastq.gz | fastx_trimmer -l $(TARGET_LENGTH) -z -o $(ANALYSIS)/$(r).2.fastq.gz" >> $(PIPELINE)$(NL))
	$(foreach r,$(LOW_CD_READGROUPS),echo "zcat $(DATA)/fastq/run_1830/$(r).1.fastq.gz | fastx_trimmer -l $(TARGET_LENGTH) -z -o $(ANALYSIS)/$(r).1.fastq.gz" >> $(PIPELINE)$(NL))
	$(foreach r,$(LOW_CD_READGROUPS),echo "zcat $(DATA)/fastq/run_1830/$(r).2.fastq.gz | fastx_trimmer -l $(TARGET_LENGTH) -z -o $(ANALYSIS)/$(r).2.fastq.gz" >> $(PIPELINE)$(NL))
	cd $(ANALYSIS) && drmrarray -j $@ -s 15 $(PIPELINE)

atacseq:
	mkdir -p $(ATACSEQ_DIR)
	python $(CONTROL)/work/atacseq/make_config.py $(ROOT) $(ATACSEQ_DIR)/results $(WORK)/trim-read-length $(WORK)/trim-read-length/*.fastq.gz > $(ATACSEQ_DIR)/config.json
	cd $(ATACSEQ_DIR) && nohup snakemake -p -j 20 --configfile $(ATACSEQ_DIR)/config.json --snakefile $(ROOT)/src/ATACseq-Snakemake/src/Snakefile --cluster-config $(ROOT)/src/ATACseq-Snakemake/examples/cluster.config --cluster "sbatch -t {cluster.time} -n {cluster.N} --mem-per-cpu={cluster.mem}" &

pca:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && Rscript $(ROOT)/src/pca.R --counts-dir $(WORK)/master_peak_counts/counts --out pca.Rda
	cd $(ANALYSIS) && Rscript $(ROOT)/src/plot_pca.R --pca pca.Rda --sample-info $(SAMPLE_INFO) --color-by tn5 --max-pc 8 --out pca.color_by_tn5.pdf
	cd $(ANALYSIS) && Rscript $(ROOT)/src/plot_pca.R --pca pca.Rda --sample-info $(SAMPLE_INFO) --color-by replicate --max-pc 8 --out pca.color_by_replicate.pdf
	cd $(ANALYSIS) && Rscript $(ROOT)/src/plot_pca.R --pca pca.Rda --sample-info $(SAMPLE_INFO) --color-by pcr_cycles --max-pc 8 --out pca.color_by_pcr_cycles.pdf
	cd $(ANALYSIS) && Rscript $(ROOT)/src/plot_pca.R --pca pca.Rda --sample-info $(SAMPLE_INFO) --color-by sequencing_run --max-pc 8 --out pca.color_by_sequencing_run.pdf

sample_info:
	@cd $(ROOT)/cluster_density/sample_info && Rscript make_sample_info.R

peak_chromhmm_overlaps:
	mkdir -p $(ANALYSIS)
	$(foreach l,$(LIBRARIES),bedtools intersect -a $(ATACSEQ_DIR)/results/macs2/$(l)_peaks.broadPeak.noblacklist -b $(DATA)/chromhmm/2013/*.bed -filenames -wo | cut -f4,10,14 | perl -pe 's@\t.*GM12878.(.*).bed@\t$$1@' | sort -k1,1 -k3n,3 | bedtools groupby -g 1 -c 2,3 -o last | cut -f2 | sort | uniq -c > $(ANALYSIS)/$(l).counts.txt$(NL))

read_tf_overlaps: LIBRARIES = $(HIGH_CD_LIBRARIES)

tn5_sensitive_peaks: LIBRARIES = $(HIGH_CD_LIBRARIES)

heatmaps: LIBRARIES = $(HIGH_CD_LIBRARIES)
