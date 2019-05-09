include $(ATAQV_HOME)/Makefile
LOCAL_ROOT=$(ROOT)/tn5
FIGURES = $(LOCAL_ROOT)/figures

LIBRARIES=$(shell cut -f5 $(SAMPLE_INFO) | awk 'NR>1' | sort | uniq)
LIBRARIES_1X=$(shell grep '1X' $(SAMPLE_INFO) | cut -f5 | sort | uniq)
NAME=$(TN5_NAME)
READGROUPS = $(shell ls $(DATA)/fastq/run_1964/* $(DATA)/fastq/run_1979/* | perl -pe 's@.*run_\d+/(.*)\.[12].fastq.gz@$$1@g' | sort | uniq)
LIBRARY2NAME=$(shell grep $(1) $(SAMPLE_INFO) | cut -f2) Tn5, $(shell grep $(1) $(SAMPLE_INFO) | cut -f4)


sample_info:
	cd $(LOCAL_ROOT)/sample_info && Rscript make_sample_info.R

# Primary, generic processing
atacseq:
	mkdir -p $(ATACSEQ_DIR)
	python $(CONTROL)/work/atacseq/make_config.py $(ROOT) $(ATACSEQ_DIR)/results $(DATA)/fastq/run_19*/*.fastq.gz > $(ATACSEQ_DIR)/config.json
	cd $(ATACSEQ_DIR) && nohup snakemake -p -j 20 --configfile $(ATACSEQ_DIR)/config.json --snakefile $(ROOT)/src/ATACseq-Snakemake/src/Snakefile --cluster-config $(ROOT)/src/ATACseq-Snakemake/examples/cluster.config --cluster "sbatch -t {cluster.time} -n {cluster.N} --mem-per-cpu={cluster.mem}" &

# PCA
pca:
	@mkdir -p $(WORK)/pca
	@cd $(WORK)/pca && Rscript $(ROOT)/src/pca.R --scale --counts-dir $(WORK)/master_peak_counts/counts --out pca.Rda
	@cd $(WORK)/pca && Rscript $(ROOT)/src/plot_pca.R --pca pca.Rda --sample-info $(SAMPLE_INFO) --color-by tn5 --max-pc 8 --out pca.color_by_tn5.pdf
	@cd $(WORK)/pca && Rscript $(ROOT)/src/plot_pca.R --pca pca.Rda --sample-info $(SAMPLE_INFO) --color-by replicate --max-pc 8 --out pca.color_by_replicate.pdf
	@cd $(WORK)/pca && Rscript $(ROOT)/src/plot_pca.R --pca pca.Rda --sample-info $(SAMPLE_INFO) --color-by experiment --max-pc 8 --out pca.color_by_experiment.pdf
