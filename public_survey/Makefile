include $(ATAQV_HOME)/Makefile
LOCAL_ROOT=$(ROOT)/public_survey
SRC=$(ROOT)/public_survey/src

atacseq: CONFIG = $(ANALYSIS)/config.json
atacseq:
	mkdir -p $(ANALYSIS)
	python $(SRC)/make_atacseq_config.py $(SAMPLE_INFO) $(ANALYSIS)/results $(DATA)/fastq/public > $(CONFIG)
	cd $(ANALYSIS) && nohup snakemake -p -j 20 --configfile $(CONFIG) --snakefile $(ROOT)/src/ATACseq-Snakemake/src/Snakefile --cluster-config $(ROOT)/src/ATACseq-Snakemake/examples/cluster.config --cluster "sbatch -t {cluster.time} -n {cluster.N} --mem-per-cpu={cluster.mem}" &

downstream:
	mkdir -p $(ANALYSIS)/results
	#python $(CONTROL)/downstream/make_config.py $(ROOT) $(WORK)/atacseq/results $(SAMPLE_INFO) $(ANALYSIS)/results > $(ANALYSIS)/config.json
	cd $(ANALYSIS) && nohup nextflow run -params-file config.json --project_colors $(LOCAL_ROOT)/project-colors.txt -with-conda $(LOCAL_ROOT)/downstream.yml -with-trace -with-report -with-dag downstream.png -resume -qs 300 $(LOCAL_ROOT)/downstream.nf & # TODO: delete key from the pipeline
