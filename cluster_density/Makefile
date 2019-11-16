include $(ATAQV_HOME)/Makefile
LOCAL_ROOT=$(ROOT)/cluster_density
SAMPLE_INFO=$(ROOT)/cluster_density/sample_info/sample_info.txt
FIGURES=$(LOCAL_ROOT)/figures
NAME=$(CD_NAME)

trim-read-length: TARGET_LENGTH = 36
trim-read-length:
	mkdir -p $(ANALYSIS)
	$(foreach r,$(HIGH_CD_READGROUPS),echo "zcat $(DATA)/fastq/run_1833/$(r).1.fastq.gz | fastx_trimmer -l $(TARGET_LENGTH) -z -o $(ANALYSIS)/$(r).1.fastq.gz" >> $(PIPELINE)$(NL))
	$(foreach r,$(HIGH_CD_READGROUPS),echo "zcat $(DATA)/fastq/run_1833/$(r).2.fastq.gz | fastx_trimmer -l $(TARGET_LENGTH) -z -o $(ANALYSIS)/$(r).2.fastq.gz" >> $(PIPELINE)$(NL))
	$(foreach r,$(LOW_CD_READGROUPS),echo "zcat $(DATA)/fastq/run_1830/$(r).1.fastq.gz | fastx_trimmer -l $(TARGET_LENGTH) -z -o $(ANALYSIS)/$(r).1.fastq.gz" >> $(PIPELINE)$(NL))
	$(foreach r,$(LOW_CD_READGROUPS),echo "zcat $(DATA)/fastq/run_1830/$(r).2.fastq.gz | fastx_trimmer -l $(TARGET_LENGTH) -z -o $(ANALYSIS)/$(r).2.fastq.gz" >> $(PIPELINE)$(NL))
	cd $(ANALYSIS) && bash $(PIPELINE)

atacseq:
	mkdir -p $(ATACSEQ_DIR)
	python $(CONTROL)/work/atacseq/make_config.py $(ROOT) $(ATACSEQ_DIR)/results $(WORK)/trim-read-length $(WORK)/trim-read-length/*.fastq.gz > $(ATACSEQ_DIR)/config.json
	cd $(ATACSEQ_DIR) && nohup snakemake -p -j 20 --configfile $(ATACSEQ_DIR)/config.json --snakefile $(ROOT)/src/ATACseq-Snakemake/src/Snakefile --cluster-config $(ROOT)/src/ATACseq-Snakemake/examples/cluster.config --cluster "sbatch -t {cluster.time} -n {cluster.N} --mem-per-cpu={cluster.mem}" &

# Everything else
downstream:
	mkdir -p $(ANALYSIS)/results
	cp $(LOCAL_ROOT)/downstream.config $(ANALYSIS)/nextflow.config
	ln -s $(LOCAL_ROOT)/bin $(ANALYSIS)/bin
	python $(CONTROL)/work/downstream/make_config.py $(LOCAL_ROOT)/sample_info/sample-info-no-exclude.txt > $(ANALYSIS)/config.json
	cd $(ANALYSIS) && nohup nextflow run -resume -params-file config.json --sample_info_high_and_low $(LOCAL_ROOT)/sample_info/sample-info-no-exclude.txt --name 'ataqv-cluster_density' --results $(ANALYSIS)/results --sample_info $(LOCAL_ROOT)/sample_info/sample-info-high-cd-no-exclude.txt -with-trace --refseq_housekeeping_tss $(ROOT)/src/ataqv/data/tss/hg19.tss.refseq.housekeeping.all.bed.gz --refseq_tss $(ROOT)/src/ataqv/data/tss/hg19.tss.refseq.bed.gz  -with-report -with-dag downstream.png -with-conda $(ROOT)/downstream.yml -qs 50 $(LOCAL_ROOT)/downstream.nf &
