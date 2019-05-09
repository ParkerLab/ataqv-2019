include $(ATAQV_HOME)/Makefile
LOCAL_ROOT=$(ROOT)/public_survey
SRC=$(ROOT)/public_survey/src

PROJECTS = $(shell cat $(SAMPLE_INFO) | cut -f6 | awk 'NR>1' | sort | uniq)
RUNS = $(shell cut -f1 $(SAMPLE_INFO) | awk 'NR>1' | sort | uniq)

BULK_LIBRARIES = $(shell grep -v 'TRUE' $(SAMPLE_INFO) | cut -f3 | grep -v experiment | sort | uniq)
SC_LIBRARIES = $(shell grep 'TRUE' $(SAMPLE_INFO) | cut -f3 | grep -v experiment | sort | uniq)
HUMAN_LIBRARIES = $(shell grep 'Homo sapiens' $(SAMPLE_INFO) | cut -f3 | sort | uniq)
MOUSE_LIBRARIES = $(shell grep 'Mus musculus' $(SAMPLE_INFO) | cut -f3 | sort | uniq)
HUMAN_BULK_LIBRARIES = $(filter $(HUMAN_LIBRARIES),$(BULK_LIBRARIES))

atacseq: CONFIG = $(ANALYSIS)/config.json
atacseq:
	mkdir -p $(ANALYSIS)
	python $(SRC)/make_atacseq_config.py $(SAMPLE_INFO) $(ANALYSIS)/results $(DATA)/fastq/public > $(CONFIG)
	cd $(ANALYSIS) && nohup snakemake -p -j 20 --configfile $(CONFIG) --snakefile $(ROOT)/src/ATACseq-Snakemake/src/Snakefile --cluster-config $(ROOT)/src/ATACseq-Snakemake/examples/cluster.config --cluster "sbatch -t {cluster.time} -n {cluster.N} --mem-per-cpu={cluster.mem}" &

ataqv:
	mkdir -p $(ANALYSIS)
	echo "# drmr:job working_directory=$(ANALYSIS)" > $(PIPELINE)
	$(foreach l,$(HUMAN_LIBRARIES),echo "ataqv --peak-file $(WORK)/atacseq/results/macs2/$(l)_peaks.broadPeak --name $(l) --metrics-file $(ANALYSIS)/$(l).ataqv.json.gz --excluded-region-file $(DATA)/mappability/hg19.blacklist.1.bed.gz --excluded-region-file $(DATA)/mappability/hg19.blacklist.2.bed.gz --tss-file $(TSS) --ignore-read-groups human $(WORK)/atacseq/results/mark_duplicates/$(l).md.bam > $(ANALYSIS)/$(l).ataqv.out" >> $(PIPELINE)$(NL))
	$(foreach l,$(MOUSE_LIBRARIES),echo "ataqv --peak-file $(WORK)/atacseq/results/macs2/$(l)_peaks.broadPeak --name $(l) --metrics-file $(ANALYSIS)/$(l).ataqv.json.gz --excluded-region-file $(DATA)/mappability/mm9.blacklist.1.bed.gz --tss-file $(MM9_TSS) --ignore-read-groups mouse $(WORK)/atacseq/results/mark_duplicates/$(l).md.bam > $(ANALYSIS)/$(l).ataqv.out" >> $(PIPELINE)$(NL))
	cd $(ANALYSIS) && drmrarray -j public-$@ -s 25 $(PIPELINE)

biosample_info: SRS = $(shell cut -f4 $(SAMPLE_INFO) | awk 'NR>1' | sort | uniq)
biosample_info:
	mkdir -p $(ANALYSIS)
	$(foreach s,$(SRS),esearch -db biosample -query $(s) | esummary > $(ANALYSIS)/$(s).xml$(NL))
	$(foreach s,$(SRS),python $(SRC)/parse_biosample_info.py $(ANALYSIS)/$(s).xml >> $(ANALYSIS)/parsed.txt$(NL))

heatmaps: LIBRARIES = $(HUMAN_BULK_LIBRARIES)

ataqv_session.public_bulk: LIBRARIES = $(BULK_LIBRARIES)
ataqv_session.public_bulk: SESSION = $(ATAQV_SESSIONS)/$(ROOT_NAME)-public-bulk
ataqv_session.public_bulk:
	mkarv -r SRR891268 $(SESSION) $(foreach l,$(LIBRARIES),$(WORK)/atacseq/results/ataqv/$(l).ataqv.json.gz)

ataqv_session.public_single_cell: LIBRARIES = $(SC_LIBRARIES)
ataqv_session.public_single_cell: SESSION = $(ATAQV_SESSIONS)/$(ROOT_NAME)-public-single-cell
ataqv_session.public_single_cell:
	mkarv -r SRR891268 $(SESSION) $(foreach l,$(LIBRARIES),$(WORK)/atacseq/results/ataqv/$(l).ataqv.json.gz)

ataqv_session.PRJ%: PROJECT = $(subst ataqv_session.,,$@)
ataqv_session.PRJ%: SESSION = $(ATAQV_SESSIONS)/$(PROJECT)
ataqv_session.PRJ%: LIBRARIES = $(shell cat $(SAMPLE_INFO) | grep $(PROJECT) | cut -f3 | sort | uniq)
ataqv_session.PRJ%:
	mkarv -r SRR891268 $(SESSION) $(foreach l,$(LIBRARIES),$(WORK)/atacseq/results/ataqv/$(l).ataqv.json.gz)

ataqv_sessions: ataqv_session.public_bulk ataqv_session.public_single_cell $(foreach p,$(PROJECTS),ataqv_session.$(p))
