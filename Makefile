ROOT=$(ATAQV_HOME)
SRC=$(ROOT)/src
TEX=$(ROOT)/tex
DATA=$(ROOT)/data
FIGURES=$(ROOT)/figures
ATAQV_SESSIONS=$(ROOT)/ataqv_sessions
ROOT_NAME = ataqv

GENOMES = mm9 hg19

TN5_NAME=$(ROOT_NAME)-tn5_series
CD_NAME=$(ROOT_NAME)-cluster_density

# Variables to be used in downstream Makefiles
ANALYSIS = $(WORK)/$@
PIPELINE = $(ANALYSIS)/pipeline
TSS=$(DATA)/tss/hg19.tss.bed.gz
CHROM_SIZES = $(DATA)/chrom_sizes/hg19.chrom.sizes
MM9_CHROM_SIZES = $(ROOT)/data/chrom_sizes/mm9.chrom.sizes
MM9_TSS=$(DATA)/tss/mm9.tss.bed.gz
WORK = $(LOCAL_ROOT)/work
CONTROL = $(LOCAL_ROOT)/control
ATACSEQ_DIR = $(WORK)/atacseq
MASTER_PEAKS = $(WORK)/master_peaks/master_peaks.bed
SAMPLE_INFO = $(LOCAL_ROOT)/sample_info/sample_info.txt
JOBNAME = $(NAME)-$@
MKARV_REFERENCE=SRR891268
MKARV_TEMPLATE=$(ROOT)/src/ataqv/src/web

.PHONY: data GEO sample_info

##### https://stackoverflow.com/questions/7039811/how-do-i-process-extremely-long-lists-of-files-in-a-make-recipe
define NL


endef
#####

SEQRUNS := 1830 1833 1964 1979

# Commands to be used in individual experiments
# Make the basic ataqv session (uninformative names, individual peak files)
ataqv_session: SESSION = $(ATAQV_SESSIONS)/$(NAME)
ataqv_session:
	mkarv -r $(MKARV_REFERENCE) -t $(MKARV_TEMPLATE) $(SESSION) $(ATACSEQ_DIR)/results/ataqv/*.ataqv.json.gz
	bash $(ROOT)/src/ataqvSessionToVS.sh $(SESSION) /lab/web/data/porchard/$(NAME)

# Make the ataqv session using master peaks
ataqv_session_master_peaks: SESSION = $(ATAQV_SESSIONS)/$(NAME)-master-peaks
ataqv_session_master_peaks:
	mkarv -r $(MKARV_REFERENCE) -t $(MKARV_TEMPLATE) $(SESSION) $(foreach l,$(LIBRARIES),$(WORK)/ataqv_master_peaks/$(l).ataqv.json.gz)
	bash $(ROOT)/src/ataqvSessionToVS.sh $(SESSION) /lab/web/data/porchard/$(NAME)

ataqv_master_peaks:
	mkdir -p $(ANALYSIS)
	echo "# drmr:job working_directory=$(ANALYSIS)" > $(PIPELINE)
	$(foreach library, $(LIBRARIES), echo "ataqv --peak-file $(MASTER_PEAKS) --name $(library) --metrics-file $(ANALYSIS)/$(library).ataqv.json.gz --excluded-region-file $(DATA)/mappability/hg19.blacklist.1.bed.gz --excluded-region-file $(DATA)/mappability/hg19.blacklist.2.bed.gz --tss-file $(TSS) --ignore-read-groups human $(ATACSEQ_DIR)/results/mark_duplicates/$(library).md.bam > $(ANALYSIS)/$(library).ataqv.out" >> $(PIPELINE)$(NL))
	cd $(ANALYSIS) && drmrarray -j $(JOBNAME) -s 8 pipeline

ataqv_informative_names:
	mkdir -p $(ANALYSIS)
	$(foreach l,$(LIBRARIES),python $(ROOT)/src/rename-ataqv.py $(WORK)/ataqv_master_peaks/$(l).ataqv.json.gz "$(call LIBRARY2NAME,$(l))" | gzip -c > $(ANALYSIS)/$(l).ataqv.json.gz$(NL))

master_peaks: NEG_LOG_10_FDR=2
master_peaks: MIN_NUMBER_1X_LIBRARIES = 2
master_peaks:
	mkdir -p $(ANALYSIS)
	cat $(foreach l,$(LIBRARIES_1X),$(ATACSEQ_DIR)/results/macs2/$(l)_peaks.broadPeak.noblacklist) | awk '$$9>=$(NEG_LOG_10_FDR)' | sort -k1,1 -k2n,2 | bedtools merge -i stdin > $(ANALYSIS)/merged.bed
	cat $(foreach l,$(LIBRARIES_1X),$(ATACSEQ_DIR)/results/macs2/$(l)_peaks.broadPeak.noblacklist) | awk '$$9>=$(NEG_LOG_10_FDR)' | sort -k1,1 -k2n,2 | bedtools intersect -a $(ANALYSIS)/merged.bed -b stdin -wa -wb | sort -k1,1 -k2n,2 | cut -f1-3,7 | perl -pe 's/_peak.*//; s/___.*//' | sort -k1,1 -k2n,2 | bedtools groupby -g 1,2,3 -o count_distinct -c 4 | awk '$$4>=$(MIN_NUMBER_1X_LIBRARIES)' | cut -f1-3 > $(MASTER_PEAKS)
	rm $(ANALYSIS)/merged.bed

master_peak_counts:
	mkdir -p $(ANALYSIS)/sort
	mkdir -p $(ANALYSIS)/counts
	echo "# drmr:job working_directory=$(ANALYSIS)" > $(PIPELINE)
	$(foreach l, $(LIBRARIES),echo "ionice -c2 -n7 samtools sort -m 3G -O bam -T $(l).sort -o $(ANALYSIS)/sort/$(l).bam $(ATACSEQ_DIR)/results/prune/$(l).pruned.bam; samtools index $(ANALYSIS)/sort/$(l).bam; ionice -c2 -n7 coverageBed -counts -sorted -a $(MASTER_PEAKS) -b $(ANALYSIS)/sort/$(l).bam > $(ANALYSIS)/counts/$(l).counts.bed; rm $(ANALYSIS)/sort/$(l).bam" >> $(PIPELINE)$(NL))
	cd $(ANALYSIS) && drmrarray -s 5 -j $(JOBNAME) $(PIPELINE)

tn5_sensitive_peaks: COUNTS = $(ANALYSIS)/counts
tn5_sensitive_peaks: SESSION = $(ATAQV_SESSIONS)/$(NAME)-master-peaks
tn5_sensitive_peaks: COVARIATES := fragment_length_distance median_fragment_length short_mononucleosomal_ratio
tn5_sensitive_peaks:
	mkdir -p $(COUNTS)
	mkdir -p $(ANALYSIS)/results
	cd $(ANALYSIS) && extractAtaqvMetric --files $(SESSION)/data/*json.gz --metrics hqaa $(COVARIATES) > metrics.txt
	cd $(COUNTS) && cp $(foreach l,$(LIBRARIES),$(WORK)/master_peak_counts/counts/$(l).counts.bed) .
	echo "# drmr:job working_directory=$(ANALYSIS)/results" > $(PIPELINE)
	$(foreach cv,$(COVARIATES),echo "Rscript $(ROOT)/src/nb_glm.R --prefix covariates_replicate___$(cv)_subsample_false --size_factors hqaa --sample_info $(SAMPLE_INFO) --counts $(COUNTS) --tn5_scale log2 --covariates_file $(ANALYSIS)/metrics.txt --covariates replicate,$(cv)" >> $(PIPELINE)$(NL))
	echo "Rscript $(SRC)/nb_glm.R --prefix covariates_replicate_subsample_false --size_factors hqaa --sample_info $(SAMPLE_INFO) --counts $(COUNTS) --tn5_scale log2 --covariates_file $(ANALYSIS)/metrics.txt --covariates replicate" >> $(PIPELINE)
	cd $(ANALYSIS) && drmrarray -j $(JOBNAME) pipeline

tss_enrichment.ENCODE:
	mkdir -p $(ANALYSIS)
	zcat $(TSS) | grep -P '^chr\d+\t' > $(ANALYSIS)/tss.bed
	echo "# drmr:job working_directory=$(ANALYSIS) memory=10g" > $(PIPELINE)
	$(foreach l,$(LIBRARIES),echo "samtools sort -m 8G -n -T $(l).sort -O bam -o $(l).bam $(ATACSEQ_DIR)/results/prune/$(l).pruned.bam; python $(ROOT)/src/tss_enrichment.py --method ENCODE --read-length 36 $(l).bam $(ANALYSIS)/tss.bed $(CHROM_SIZES) > $(l).tss_coverage.txt" >> $(PIPELINE)$(NL))
	cd $(ANALYSIS) && drmrarray -s 10 -j $(JOBNAME) $(PIPELINE)

tss_enrichment.cutsite:
	mkdir -p $(ANALYSIS)
	zcat $(TSS) | grep -P '^chr\d+\t' > $(ANALYSIS)/tss.bed
	echo "# drmr:job working_directory=$(ANALYSIS) memory=10g" > $(PIPELINE)
	$(foreach l,$(LIBRARIES),echo "samtools sort -m 8G -n -T $(l).sort -O bam -o $(l).bam $(ATACSEQ_DIR)/results/prune/$(l).pruned.bam; python $(ROOT)/src/tss_enrichment.py --method ENCODE --read-length 2 $(l).bam $(ANALYSIS)/tss.bed $(CHROM_SIZES) > $(l).tss_coverage.txt" >> $(PIPELINE)$(NL))
	cd $(ANALYSIS) && drmrarray -s 10 -j $@ $(PIPELINE)

tss_enrichment.ataqv:
	mkdir -p $(ANALYSIS)
	zcat $(TSS) | grep -P '^chr\d+\t' > $(ANALYSIS)/tss.bed
	echo "# drmr:job working_directory=$(ANALYSIS) memory=10g" > $(PIPELINE)
	$(foreach l,$(LIBRARIES),echo "samtools sort -m 8G -n -T $(l).sort -O bam -o $(l).bam $(ATACSEQ_DIR)/results/prune/$(l).pruned.bam; python $(ROOT)/src/tss_enrichment.py --method ataqv $(l).bam $(ANALYSIS)/tss.bed $(CHROM_SIZES) > $(l).tss_coverage.txt" >> $(PIPELINE)$(NL))
	cd $(ANALYSIS) && drmrarray -s 10 -j $@ $(PIPELINE)

normalization:
	mkdir -p $(ANALYSIS)/bedgraph
	mkdir -p $(ANALYSIS)/bigwig
	echo "# drmr:job working_directory=$(ANALYSIS)" > $(PIPELINE)
	$(foreach l,$(LIBRARIES),echo "python $(ROOT)/src/normalize_bedgraph/src/normalize_bedgraph.py --to-number-reads 10000000 $(ATACSEQ_DIR)/results/macs2/$(l)_treat_pileup.bdg | LC_COLLATE=C sort -k1,1 -k2,2n > $(ANALYSIS)/bedgraph/$(l).normalized.bdg; bedGraphToBigWig $(ANALYSIS)/bedgraph/$(l).normalized.bdg $(CHROM_SIZES) $(ANALYSIS)/bigwig/$(l).normalized.bw" >> $(PIPELINE)$(NL))
	cd $(ANALYSIS) && drmrarray -s 5 -j $(JOBNAME) $(PIPELINE)

read_chromhmm_overlaps:
	mkdir -p $(ANALYSIS)
	# Creating the merged ChromHMM file...
	$(foreach state,$(shell ls $(DATA)/chromhmm/2013/GM12878.* | perl -pe 's@.*GM12878\.(.*)\.bed@$$1@'),cat $(DATA)/chromhmm/2013/GM12878.$(state).bed | perl -pe 's/\n/\t$(state)\n/' >> $(ANALYSIS)/chromhmm.unsorted.bed$(NL))
	cat $(ANALYSIS)/chromhmm.unsorted.bed | sort  -k1,1 -k2n,2 > $(ANALYSIS)/chromhmm.bed
	rm $(ANALYSIS)/chromhmm.unsorted.bed
	# Creating the pipeline...
	echo "# drmr:job working_directory=$(ANALYSIS) memory=10g" > $(PIPELINE)
	$(foreach l,$(LIBRARIES),echo "bedtools intersect -a $(ATACSEQ_DIR)/results/prune/$(l).pruned.bam -b $(shell ls $(DATA)/mappability/hg19.*.bed.gz) -v | bedtools intersect -wao -a stdin -b $(ANALYSIS)/chromhmm.bed -bed | cut -f4,16,17 | sort -k1,1 -k2,2 | bedtools groupby -g 1,2 -c 3 -o sum | sort -k1,1 -k3n,3 | bedtools groupby -g 1 -c 2 -o last | cut -f2 | sort | uniq -c > $(ANALYSIS)/$(l).counts.bed" >> $(PIPELINE)$(NL))
	cd $(ANALYSIS) && drmrarray -s 10 -j $(JOBNAME) $(PIPELINE)

read_tf_overlaps: CHIPSEQ_PEAKS := $(notdir $(shell ls $(DATA)/chipseq_idr/*bed.gz | perl -pe 's/.bed.gz$$//'))
read_tf_overlaps:
	rm -rf $(ANALYSIS)
	mkdir -p $(ANALYSIS)
	echo "# drmr:job working_directory=$(ANALYSIS) memory=25g" > $(PIPELINE)
	$(foreach l,$(LIBRARIES),echo "ionice -c2 -n7 bedtools intersect -a $(ATACSEQ_DIR)/results/prune/$(l).pruned.bam -b $(shell ls $(DATA)/mappability/hg19.blacklist.*.bed.gz) -v | ionice -c2 -n7 bedtools intersect -a stdin -b $(DATA)/chipseq_idr/*.bed.gz -filenames -bed -loj | cut -f4,13 | perl -pe 's@$(DATA)/chipseq_idr/(.*).bed.gz@\$$1@' | bedtools groupby -g 1 -c 2 -o distinct | bedtools groupby -g 2 -c 1 -o count_distinct | sort -S 10G -k1,1 | bedtools groupby -g 1 -c 2 -o sum > $(l).intermediate; bedtools expand -c 1 -i $(l).intermediate | sort -S 10G -k1,1 | bedtools groupby -g 1 -c 2 -o sum | awk '\$$1!=\".\"' > $(l).tf_overlap.txt; cat $(l).intermediate | perl -pe 's@\\\n@\\\t.\\\n@' | bedtools groupby -g 3 -c 2 -o sum | cut -f2 > $(l).total_reads; rm $(l).intermediate" >> $(PIPELINE)$(NL))
	cd $(ANALYSIS) && drmrarray -s 7 -j $(JOBNAME) $(PIPELINE)

heatmaps: GENES = GAPDH VCP
heatmaps:
	mkdir -p $(ANALYSIS)
	$(foreach g,$(GENES),zcat $(DATA)/gtf/hg19.gtf.gz | grep -i 'gene_name "$(g)"' | cut -f1,4,5 | sort -k1,1 | bedtools groupby -g 1 -c 2,3 -o min,max > $(ANALYSIS)/$(g).coords.bed$(NL))
	$(foreach g,$(GENES),cat $(ANALYSIS)/$(g).coords.bed | bedtools slop -i stdin -g $(CHROM_SIZES) -b 15000 | bedtools makewindows -b stdin -w 20 > $(ANALYSIS)/$(g).windows.bed$(NL))
	cat $(foreach g,$(GENES),$(ANALYSIS)/$(g).coords.bed) | bedtools slop -i stdin -g $(CHROM_SIZES) -b 15000 > $(ANALYSIS)/keep.bed
	echo "# drmr:job memory=10g" > $(PIPELINE)
	$(foreach l,$(LIBRARIES),echo "samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 30 -L $(ANALYSIS)/keep.bed $(WORK)/atacseq/results/mark_duplicates/$(l).md.bam | bedtools bamtobed | perl -lane 'if(\$$F[5]==\"-\"){\$$F[1]=\$$F[2]-1};\$$F[2]=\$$F[1]+1;print join(\"\\\t\", @F)' | sort -k1,1 -k2n,2 > $(ANALYSIS)/$(l).bed; $(foreach g,$(GENES),bedtools coverage -counts -a $(ANALYSIS)/$(g).windows.bed -b $(ANALYSIS)/$(l).bed > $(l).$(g)_coverage.bed;)" >> $(PIPELINE)$(NL))
	cd $(ANALYSIS) && drmrarray -j $@ -s 12 $(PIPELINE)

# DATA
data: tss gtf fasta bwa chrom_sizes chromhmm mappability fastq chipseq_idr

tss: ANALYSIS = $(DATA)/tss
tss: GENOMES = mm9 hg19 rn5
tss:
	mkdir -p $(ANALYSIS)
	$(foreach g,$(GENOMES),zcat $(SRC)/ataqv/data/tss/$(g).tss.refseq.housekeeping.ortho.bed.gz > $(ANALYSIS)/$(g).tss.bed$(NL))
	$(foreach g,$(GENOMES),cat $(ANALYSIS)/$(g).tss.bed | cut -f4 | sort | uniq -u > $(ANALYSIS)/$(g).single_tss_genes.bed$(NL))
	cat $(ANALYSIS)/*.single_tss_genes.bed | perl -pe 'tr/[a-z]/[A-Z]/' | sort | uniq -c | awk '$$1==3' | awk '{print($$2)}' > $(ANALYSIS)/single_in_3.txt
	$(foreach g,$(GENOMES),grep -w -i -f $(ANALYSIS)/single_in_3.txt $(ANALYSIS)/$(g).tss.bed | sort -k1,1 -k2n,2 | gzip -c > $(ANALYSIS)/$(g).tss.bed.gz$(NL))
	rm $(ANALYSIS)/*.single_tss_genes.bed $(ANALYSIS)/single_in_3.txt $(foreach g,$(GENOMES),$(ANALYSIS)/$(g).tss.bed)

gtf:
	mkdir -p $(DATA)/$@ && cd $(DATA)/$@ && wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz && ln -s gencode.v19.annotation.gtf.gz hg19.gtf.gz

fasta:
	$(foreach g,$(GENOMES),mkdir -p $(DATA)/$@/$(g)$(NL))
	$(foreach g,$(GENOMES),cd $(DATA)/$@/$(g) && wget -O $(g).tar.gz ftp://hgdownload.cse.ucsc.edu/goldenPath/$(g)/bigZips/chromFa.tar.gz && tar -xvzf $(g).tar.gz && rename 's/^/$(g)./g' chr*.fa && rename 's/$$/sta/g' $(g).*fa && cat $(g).chr* > $(g).fa && rm $(g)*.fasta $(g).tar.gz$(NL))

bwa:  ANALYSIS = $(DATA)/$@
bwa: 
	mkdir -p $(DATA)/$@
	echo "# drmr:job working_directory=$(ANALYSIS) memory=10g" > $(PIPELINE)
	$(foreach g,$(GENOMES),mkdir -p $(ANALYSIS)/$(g)$(NL))
	$(foreach g,$(GENOMES),echo "cd $(ANALYSIS)/$(g) && ln -s $(DATA)/fasta/$(g)/$(g).fa $(g) && bwa index $(g)" >> $(PIPELINE)$(NL))
	cd $(ANALYSIS) && drmrarray -j $(ROOT_NAME)-$@ $(PIPELINE)

chipseq_idr: ANALYSIS = $(DATA)/$@
chipseq_idr:
	mkdir -p $(ANALYSIS)
	cat $(ROOT)/manuscript/TableS3/TableS3.csv | perl -pe 's/,/\t/g' | awk 'NR>1' | awk '{print($$4, $$3, $$1)}' | perl -pe 's/(.*) (.*) (ENCFF.*)$$/wget $$1 -O $$2.$$3.bed.gz/' >> $(PIPELINE)
	cd $(ANALYSIS) && bash $(PIPELINE)


chrom_sizes:
	mkdir -p $(DATA)/$@
	cd $(DATA)/$@ && $(foreach g,$(GENOMES),wget http://hgdownload.cse.ucsc.edu/goldenPath/$(g)/bigZips/$(g).chrom.sizes;)

chromhmm:
	@mkdir -p data/chromhmm/2013
	@wget https://research.nhgri.nih.gov/manuscripts/Collins/islet_chromatin/hg19/ChromHMM/GM12878_chromHMM.bb
	@bigBedToBed GM12878_chromHMM.bb GM12878_chromHMM.bed
	@cat GM12878_chromHMM.bed | perl -pe 's/255,252,4/weak_enhancer/' > GM12878_chromHMM.bed.tmp; mv GM12878_chromHMM.bed.tmp GM12878_chromHMM.bed
	@cat GM12878_chromHMM.bed | perl -pe 's/250,202,0/strong_enhancer/' > GM12878_chromHMM.bed.tmp; mv GM12878_chromHMM.bed.tmp GM12878_chromHMM.bed
	@cat GM12878_chromHMM.bed | perl -pe 's/0,176,80/transcribed/' > GM12878_chromHMM.bed.tmp; mv GM12878_chromHMM.bed.tmp GM12878_chromHMM.bed
	@cat GM12878_chromHMM.bed | perl -pe 's/10,190,254/insulator/' > GM12878_chromHMM.bed.tmp; mv GM12878_chromHMM.bed.tmp GM12878_chromHMM.bed
	@cat GM12878_chromHMM.bed | perl -pe 's/255,105,105/weak_promoter/' > GM12878_chromHMM.bed.tmp; mv GM12878_chromHMM.bed.tmp GM12878_chromHMM.bed
	@cat GM12878_chromHMM.bed | perl -pe 's/255,0,0/active_promoter/' > GM12878_chromHMM.bed.tmp; mv GM12878_chromHMM.bed.tmp GM12878_chromHMM.bed
	@cat GM12878_chromHMM.bed | perl -pe 's/207,11,198/poised_promoter/' > GM12878_chromHMM.bed.tmp; mv GM12878_chromHMM.bed.tmp GM12878_chromHMM.bed
	@cat GM12878_chromHMM.bed | perl -pe 's/127,127,127/repressed/' > GM12878_chromHMM.bed.tmp; mv GM12878_chromHMM.bed.tmp GM12878_chromHMM.bed
	@cat GM12878_chromHMM.bed | perl -pe 's/255,255,255/low_signal/' > GM12878_chromHMM.bed.tmp; mv GM12878_chromHMM.bed.tmp GM12878_chromHMM.bed
	@for i in `cut -f9 GM12878_chromHMM.bed | sort | uniq`; do grep $$i GM12878_chromHMM.bed | cut -f1-3 > data/chromhmm/2013/GM12878.$$i.bed; done
	@rm GM12878_chromHMM.bed
	@mkdir -p data/chromhmm/tracks/2013
	@mv GM12878_chromHMM.bb data/chromhmm/tracks/2013

mappability:
	mkdir -p $(DATA)/$@
	cd $(DATA)/$@ && wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz && wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz && ln -s wgEncodeDacMapabilityConsensusExcludable.bed.gz hg19.blacklist.1.bed.gz && ln -s wgEncodeDukeMapabilityRegionsExcludable.bed.gz hg19.blacklist.2.bed.gz
	cd $(DATA)/$@ && wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm9-mouse/mm9-blacklist.bed.gz && ln -s mm9-blacklist.bed.gz mm9.blacklist.1.bed.gz

fastq: PUBLIC_RUNS = $(shell cut -f1 $(ROOT)/public_survey/sample_info/sample_info.txt | awk 'NR>1' | sort | uniq)
fastq:
	$(foreach r,1964 1979 1833 1830,mkdir -p $(DATA)/$@/run_$(r)$(NL))
	cd $(DATA)/$@/run_1964 && ln -s /lab/data/seqcore/Run_1964/fastq/* . && rename 's/(^\d+)___(\d+)___L(\d+)___GM12878.*.([12]).fq.gz/$$1___L$$3___1964.$$4.fastq.gz/' *
	cd $(DATA)/$@/run_1979 && ln -s /lab/data/seqcore/Run_1979/fastq/* . && rename 's/(^\d+)___(\d+)___L(\d+)___GM12878.*.([12]).fq.gz/$$1___L$$3___1979.$$4.fastq.gz/' *
	cd $(DATA)/$@/run_1830 && ln -s /lab/data/seqcore/Run_1830/fastq/* . && rename 's/^.*low___(\d+)___(L\d+)___GM12878.*\.([12]).fq.gz/$$1___$$2___1830.$$3.fastq.gz/' *
	cd $(DATA)/$@/run_1833 && ln -s /lab/data/seqcore/Run_1833/fastq/* . && rename 's/^.*high___(\d+)___(L\d+)___GM12878.*\.([12]).fq.gz/$$1___$$2___1833.$$3.fastq.gz/' *
	mkdir -p $(DATA)/$@/public
	echo "# drmr:job memory=2g working_directory=$(DATA)/$@/public" > $(DATA)/$@/public/pipeline
	$(foreach r,$(PUBLIC_RUNS),echo "fastq-dump -F --gzip --split-files $(r) && mv $(r)_1.fastq.gz $(r).1.fastq.gz && mv $(r)_2.fastq.gz $(r).2.fastq.gz" >> $(DATA)/$@/public/pipeline$(NL))
	cd $(DATA)/$@/public && drmrarray -s 10 -j $@ $(PIPELINE)

compare-tss-lists: HOUSEKEEPING_TSS = $(ROOT)/src/ataqv/data/tss/hg19.tss.refseq.housekeeping.all.bed.gz
compare-tss-lists: NONHOUSEKEEPING_TSS = $(ROOT)/src/ataqv/data/tss/hg19.tss.refseq.bed.gz
compare-tss-lists:
	mkdir -p $(ANALYSIS)
	echo "# drmr:job working_directory=$(ANALYSIS)" > $(PIPELINE)
	$(foreach l,$(LIBRARIES),echo "ataqv --peak-file $(WORK)/atacseq/results/macs2/$(l)_peaks.broadPeak --name $(l)-housekeeping --metrics-file $(ANALYSIS)/$(l)-housekeeping.ataqv.json.gz --excluded-region-file $(DATA)/mappability/hg19.blacklist.2.bed.gz --excluded-region-file $(DATA)/mappability/hg19.blacklist.1.bed.gz --tss-file $(HOUSEKEEPING_TSS) --ignore-read-groups human $(WORK)/atacseq/results/mark_duplicates/$(l).md.bam > $(ANALYSIS)/$(l)-housekeeping.ataqv.out" >> $(PIPELINE)$(NL))
	$(foreach l,$(LIBRARIES),echo "ataqv --peak-file $(WORK)/atacseq/results/macs2/$(l)_peaks.broadPeak --name $(l)-all --metrics-file $(ANALYSIS)/$(l)-all.ataqv.json.gz --excluded-region-file $(DATA)/mappability/hg19.blacklist.2.bed.gz --excluded-region-file $(DATA)/mappability/hg19.blacklist.1.bed.gz --tss-file $(NONHOUSEKEEPING_TSS) --ignore-read-groups human $(WORK)/atacseq/results/mark_duplicates/$(l).md.bam > $(ANALYSIS)/$(l)-all.ataqv.out" >> $(PIPELINE)$(NL))
	cd $(ANALYSIS) && drmrarray -j $@ $(PIPELINE)


# Figures
# For the manuscript
EXAMPLE_PROMOTER_PEAK=chr1:40348769-40349427
EXAMPLE_ENHANCER_PEAK=chr12:125389250-125389658

###### MAIN FIGURES
figures.public: ANALYSIS = $(FIGURES)/public
figures.public: SESSION = $(ATAQV_SESSIONS)/$(ROOT_NAME)-public-bulk
figures.public:
	rm -rf $(ANALYSIS)
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && python2 $(SRC)/extractFragmentLengths.py --files $(SESSION)/data/*json.gz > fragment_length_distributions.txt
	cd $(ANALYSIS) && python2 $(SRC)/extractTssCoverage.py --files $(SESSION)/data/*json.gz > tss_coverage.txt
	cd $(ANALYSIS) && python2 $(SRC)/extractAtaqvMetric.py --metrics total_reads duplicate_reads hqaa percent_mitochondrial median_fragment_length total_peaks fragment_length_distance hqaa_overlapping_peaks_percent tss_enrichment --files $(SESSION)/data/*.json.gz > metrics.txt
	cd $(ANALYSIS) && python2 $(SRC)/extractAtaqvMetric.py --metrics total_reads duplicate_reads hqaa percent_mitochondrial median_fragment_length total_peaks fragment_length_distance hqaa_overlapping_peaks_percent tss_enrichment --files $(ROOT)/public_survey/work/ataqv_sessions/multi-flowcell/data/*.json.gz > metrics_flowcell.txt


figures.tn5: ANALYSIS = $(FIGURES)/tn5
figures.tn5: SESSION = $(ATAQV_SESSIONS)/$(TN5_NAME)-master-peaks
figures.tn5: GB_SESSION = ataqv-manuscript-2018-07-18
figures.tn5: CHIPSEQ = $(patsubst %.bed.gz,%,$(notdir $(wildcard $(DATA)/chipseq_idr/*.bed.gz)))
figures.tn5: MASTER_PEAKS=$(ROOT)/tn5/work/master_peaks/master_peaks.bed
figures.tn5:
	mkdir -p $(ANALYSIS) && cd $(ANALYSIS) && python $(SRC)/extractFragmentLengths.py --files $(SESSION)/data/* > fragment_length_distributions.txt
	mkdir -p $(ANALYSIS) && cd $(ANALYSIS) && python $(SRC)/extractTssCoverage.py --files $(SESSION)/data/* > tss_coverage.txt
	cd $(ANALYSIS) && python $(SRC)/extractAtaqvMetric.py --metrics total_reads duplicate_reads hqaa percent_mitochondrial median_fragment_length total_peaks short_mononucleosomal_ratio fragment_length_distance hqaa_overlapping_peaks_percent tss_enrichment --files $(SESSION)/data/* > metrics.txt
	cat $(DATA)/qpcr/2017-08-10.tsv | awk 'NR>36' | perl -pe 's@2/3@0.66@g; s@1/2@0.5@g; s@1/5@0.2@g' > $(ANALYSIS)/2017-08-10.txt
	cat $(DATA)/qpcr/2017-08-14.tsv | awk 'NR>33' | perl -pe 's@2/3@0.66@g; s@1/2@0.5@g; s@1/5@0.2@g' > $(ANALYSIS)/2017-08-14.txt
	mkdir -p $(ANALYSIS) && cd $(ANALYSIS) && $(foreach chip,$(CHIPSEQ),bedtools intersect -c -a $(MASTER_PEAKS) -b $(DATA)/chipseq_idr/$(chip).bed.gz > $(chip).peak_overlap.bed;)
	mkarv $(ROOT)/tn5/work/compare-tss-lists/session $(ROOT)/tn5/work/compare-tss-lists/*.json.gz
	python $(SRC)/extractTssCoverage.py --files $(ROOT)/tn5/work/compare-tss-lists/session/data/* > $(ANALYSIS)/tss-list-comparison.tss_coverage.txt
	python $(SRC)/extractAtaqvMetric.py --metrics tss_enrichment --files $(ROOT)/tn5/work/compare-tss-lists/session/data/* > $(ANALYSIS)/tss-list-comparison.tss_enrichment.txt


figures.cd: ANALYSIS = $(FIGURES)/cd
figures.cd: SESSION = $(ATAQV_SESSIONS)/$(CD_NAME)-master-peaks
figures.cd:
	rm -rf $(ANALYSIS)
	mkdir -p $(ANALYSIS) && cd $(ANALYSIS) && python $(SRC)/extractFragmentLengths.py --files $(SESSION)/data/* > fragment_length_distributions.txt
	mkdir -p $(ANALYSIS) && cd $(ANALYSIS) && python $(SRC)/extractTssCoverage.py --files $(SESSION)/data/* > tss_coverage.txt
	cd $(ANALYSIS) && python $(SRC)/extractAtaqvMetric.py --metrics total_reads duplicate_reads hqaa percent_mitochondrial median_fragment_length total_peaks short_mononucleosomal_ratio fragment_length_distance hqaa_overlapping_peaks_percent tss_enrichment --files $(SESSION)/data/* > metrics.txt

figures.supplement.gb_screenshots: ANALYSIS = $(FIGURES)/gb
figures.supplement.gb_screenshots: REPS := $(shell seq 1 6)
figures.supplement.gb_screenshots:
	rm -rf $(ANALYSIS)
	mkdir -p $(ANALYSIS)
	rm -rf $(ANALYSIS)/*.pdf
	cd $(ANALYSIS) && echo "$(EXAMPLE_PROMOTER_PEAK)" | perl -pe 's/[:-]/\t/g' > region_of_interest.bed
	$(foreach rep, $(REPS), cd $(ANALYSIS) && python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session ataqv-rep-$(rep) --bed region_of_interest.bed --prefix promoter_rep$(rep) --viewer-width 400 --label-width 14 --extension 1700 --text-size 18$(NL))
	cd $(ANALYSIS) && echo "$(EXAMPLE_ENHANCER_PEAK)" | perl -pe 's/[:-]/\t/g' > region_of_interest.bed
	$(foreach rep, $(REPS), cd $(ANALYSIS) && python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session ataqv-rep-$(rep) --bed region_of_interest.bed --prefix enhancer_rep$(rep) --viewer-width 400 --label-width 14 --extension 1000 --text-size 18$(NL))
	rm $(ANALYSIS)/region_of_interest.bed
	cd $(ANALYSIS) && $(foreach rep,$(REPS),mv promoter_rep$(rep).*pdf $(ROOT)/knitr/promoter_rep$(rep).pdf;)
	cd $(ANALYSIS) && $(foreach rep,$(REPS),mv enhancer_rep$(rep).*pdf $(ROOT)/knitr/enhancer_rep$(rep).pdf;)
	
knit: figures.public figures.tn5 figures.cd figures.supplement.gb_screenshot
	cd knitr && Rscript -e "library(knitr); knit('main.Rnw')" && pdflatex main.tex && pdflatex main.tex

TableS2: SAMPLE_INFO = $(ROOT)/public_survey/sample_info/sample_info.txt
TableS2: ANALYSIS = $(FIGURES)/manuscript/$@
TableS2:
	@echo "Creating $@"
	@rm -rf $(ANALYSIS)
	@mkdir -p $(ANALYSIS)
	@echo "library,run,organism,is_single_cell" > $(ANALYSIS)/$@.csv
	@cat $(SAMPLE_INFO) | grep -v is_single_cell | awk -F '\t' '{print($$3"\t"$$1"\t"$$5"\t"$$8)}' | perl -pe 's/\t/,/g; s/TRUE/yes/; s/FALSE/no/' >> $(ANALYSIS)/$@.csv
