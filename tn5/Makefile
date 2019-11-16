include $(ATAQV_HOME)/Makefile
LOCAL_ROOT=$(ROOT)/tn5
FIGURES = $(LOCAL_ROOT)/figures

NAME=$(TN5_NAME)
LIBRARY2NAME=$(shell grep $(1) $(SAMPLE_INFO) | cut -f2) Tn5, $(shell grep $(1) $(SAMPLE_INFO) | cut -f4)

# Primary, generic processing
atacseq:
	mkdir -p $(ATACSEQ_DIR)
	python $(CONTROL)/work/atacseq/make_config.py $(ROOT) $(ATACSEQ_DIR)/results $(DATA)/fastq/run_19*/*.fastq.gz > $(ATACSEQ_DIR)/config.json
	cd $(ATACSEQ_DIR) && nohup snakemake -p -j 20 --configfile $(ATACSEQ_DIR)/config.json --snakefile $(ROOT)/src/ATACseq-Snakemake/src/Snakefile --cluster-config $(ROOT)/src/ATACseq-Snakemake/examples/cluster.config --cluster "sbatch -t {cluster.time} -n {cluster.N} --mem-per-cpu={cluster.mem}" &

# Everything else
downstream:
	mkdir -p $(ANALYSIS)/results
	cp $(LOCAL_ROOT)/downstream.config $(ANALYSIS)/nextflow.config
	ln -s $(LOCAL_ROOT)/bin $(ANALYSIS)/bin
	python $(CONTROL)/work/downstream/make_config.py $(SAMPLE_INFO) > $(ANALYSIS)/config.json
	cd $(ANALYSIS) && nohup nextflow run -resume -params-file config.json --name 'ataqv-tn5_series' --results $(ANALYSIS)/results --sample_info $(SAMPLE_INFO) -with-trace --refseq_housekeeping_tss $(ROOT)/src/ataqv/data/tss/hg19.tss.refseq.housekeeping.all.bed.gz --refseq_tss $(ROOT)/src/ataqv/ata/tss/hg19.tss.refseq.bed.gz  -with-report -with-dag downstream.png -with-conda $(ROOT)/downstream.yml -qs 50 $(LOCAL_ROOT)/downstream.nf &
