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
WORK = $(LOCAL_ROOT)/work
CONTROL = $(LOCAL_ROOT)/control
ATACSEQ_DIR = $(WORK)/atacseq
SAMPLE_INFO = $(LOCAL_ROOT)/sample_info/sample_info.txt
JOBNAME = $(NAME)-$@
MKARV_REFERENCE=SRR891268
MKARV_TEMPLATE=$(ROOT)/src/ataqv/src/web

.PHONY: data sample_info

##### https://stackoverflow.com/questions/7039811/how-do-i-process-extremely-long-lists-of-files-in-a-make-recipe
define NL


endef
#####

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
	cd $(ANALYSIS) && bash $(PIPELINE)

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
	cd $(DATA)/$@/public && bash $(PIPELINE)



figure-1-panels: ANALYSIS=$(ROOT)/knitr/$@
figure-1-panels:
	rm -rf $(ANALYSIS)
	mkdir -p $(ANALYSIS)
	ln -s $(ROOT)/public_survey/work/downstream/results/figures/all-corr-heatmap.pdf $(ANALYSIS)/public-ataqv-metric-correlations.pdf
	ln -s $(ROOT)/public_survey/work/downstream/results/figures/public_meta.pdf $(ANALYSIS)/
	ln -s $(ROOT)/public_survey/work/downstream/results/figures/public_median_fragment_length_vs_tss_enrichment.pdf $(ANALYSIS)/
	ln -s $(ROOT)/public_survey/work/downstream/results/figures/GAPDH-heatmap.pdf $(ANALYSIS)/
	ln -s $(ROOT)/public_survey/work/downstream/results/figures/VCP-heatmap.pdf $(ANALYSIS)/
	ln -s $(ROOT)/public_survey/work/downstream/results/figures/fig1-hqaa-vs-max-fraction-color-by-chrom.pdf $(ANALYSIS)/
	ln -s $(ROOT)/public_survey/work/downstream/results/figures/compare-chromosome-coverages.pdf $(ANALYSIS)/
	ln -s $(ROOT)/public_survey/work/downstream/results/pc1-vs-ataqv-metrics/PRJNA259243-hg19.cor_barplot.pdf $(ANALYSIS)/

figure-2-panels: ANALYSIS=$(ROOT)/knitr/$@
figure-2-panels:
	rm -rf $(ANALYSIS)
	mkdir -p $(ANALYSIS)
	ln -s $(ROOT)/tn5/work/downstream/results/figures/tn5-GAPDH-heatmap.pdf $(ANALYSIS)/
	ln -s $(ROOT)/tn5/work/downstream/results/figures/tn5_chromatin_state_overlap.pdf $(ANALYSIS)/
	ln -s $(ROOT)/tn5/work/downstream/results/figures/tn5_hqaa_overlapping_peaks_percent.pdf $(ANALYSIS)/
	ln -s $(ROOT)/tn5/work/downstream/results/figures/tn5_tss_enrichment.pdf $(ANALYSIS)/
	ln -s $(ROOT)/tn5/work/downstream/results/figures/tn5_fragment_length_distributions.pdf $(ANALYSIS)/
	ln -s $(ROOT)/knitr/gb-screenshots/promoter_rep6.pdf $(ANALYSIS)/
	ln -s $(ROOT)/knitr/gb-screenshots/enhancer_rep6.pdf $(ANALYSIS)/
	ln -s $(ROOT)/knitr/manual/tn5_cartoon.pdf $(ANALYSIS)/

	
knit: ANALYSIS=$(ROOT)/knitr
knit:
	# link in completed plots
	ln -sf $(ROOT)/knitr/manual/* $(ANALYSIS)/
	ln -sf $(ROOT)/knitr/gb-screenshots/* $(ANALYSIS)/
	ln -sf $(ROOT)/public_survey/work/downstream/results/figures/example_fld_calculation_density.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/public_survey/work/downstream/results/figures/example_fld_calculation_cumulative.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/public_survey/work/downstream/results/figures/public_intrastudy_heterogeneity_fragment_length_distributions.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/public_survey/work/downstream/results/figures/public_intrastudy_heterogeneity_tss_enrichment.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/public_survey/work/downstream/results/max-fraction-from-single-autosome/bulk.hg19.hqaa-vs-max-fraction-group-by-cell-type.pdf $(ANALYSIS)/public-bulk-hg19-max-fraction.pdf
	ln -sf $(ROOT)/public_survey/work/downstream/results/max-fraction-from-single-autosome/bulk.mm9.hqaa-vs-max-fraction-group-by-cell-type.pdf $(ANALYSIS)/public-bulk-mm9-max-fraction.pdf
	ln -sf $(ROOT)/public_survey/work/downstream/results/pc1-vs-ataqv-metrics/PRJNA259243-hg19.tss_enrichment_scatter.pdf $(ROOT)/knitr/
	ln -sf $(ROOT)/public_survey/work/downstream/results/su/differential-peak-analysis/su-et-al-qqplot.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/public_survey/work/downstream/results/su/differential-peak-analysis/su_et_al_median_fragment_length_vs_batch.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/public_survey/work/downstream/results/max-fraction-from-single-autosome/public-bulk-hg19.hqaa-vs-max-fraction-group-by-cell-type.pdf $(ROOT)/knitr/public-bulk-hg19-hqaa-vs-max-fraction.pdf
	ln -sf $(ROOT)/public_survey/work/downstream/results/max-fraction-from-single-autosome/public-bulk-mm9.hqaa-vs-max-fraction-group-by-cell-type.pdf $(ROOT)/knitr/public-bulk-mm9-hqaa-vs-max-fraction.pdf
	ln -sf $(ROOT)/cluster_density/work/downstream/results/figures/cd_high_vs_low_median_fragment_length.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/cluster_density/work/downstream/results/figures/cd_high_vs_low_tss_enrichment.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/tn5/work/downstream/results/figures/tn5_median_fragment_length_all_six.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/tn5/work/downstream/results/figures/tn5_tss_enrichment_all_six.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/tn5/work/downstream/results/figures/tn5_mitochondrial_percent.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/tn5/work/downstream/results/figures/tn5_hqaa_overlapping_peaks_percent.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/tn5/work/downstream/results/figures/tn5_proportion_duplicate.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/tn5/work/downstream/results/figures/tn5_proportion_hqaa.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/cluster_density/work/downstream/results/figures/tn5_median_fragment_length_all_six.pdf $(ANALYSIS)/cd_median_fragment_length.pdf
	ln -sf $(ROOT)/cluster_density/work/downstream/results/figures/tn5_tss_enrichment_all_six.pdf $(ANALYSIS)/cd_tss_enrichment.pdf
	ln -sf $(ROOT)/cluster_density/work/downstream/results/figures/tn5_mitochondrial_percent.pdf $(ANALYSIS)/cd_mitochondrial_percent.pdf
	ln -sf $(ROOT)/cluster_density/work/downstream/results/figures/tn5_hqaa_overlapping_peaks_percent.pdf $(ANALYSIS)/cd_hqaa_overlapping_peaks_percent.pdf
	ln -sf $(ROOT)/cluster_density/work/downstream/results/figures/tn5_proportion_duplicate.pdf $(ANALYSIS)/cd_proportion_duplicate.pdf
	ln -sf $(ROOT)/cluster_density/work/downstream/results/figures/tn5_proportion_hqaa.pdf $(ANALYSIS)/cd_proportion_hqaa.pdf
	ln -sf $(ROOT)/tn5/work/downstream/results/figures/number-of-peaks-full.pdf $(ANALYSIS)/pcr-constant-number-of-peaks.pdf
	ln -sf $(ROOT)/tn5/work/downstream/results/figures/number-of-peaks-subsampled.pdf $(ANALYSIS)/pcr-constant-number-of-peaks-after-subsampling-bams.pdf
	ln -sf $(ROOT)/tn5/work/downstream/results/figures/pairwise-jaccard.mean-pairwise-jaccard-heatmap.pdf $(ANALYSIS)/pcr-constant-pairwise-peak-jaccard.pdf
	ln -sf $(ROOT)/tn5/work/downstream/results/subsample/figures/pairwise-jaccard.subsample.mean-pairwise-jaccard-heatmap.pdf $(ANALYSIS)/pcr-constant-pairwise-peak-jaccard-after-subsampling-bams.pdf
	ln -sf $(ROOT)/tn5/work/downstream/results/figures/tn5_pca.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/tn5/work/downstream/results/figures/tn5_modeling.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/cluster_density/work/downstream/results/figures/tn5_modeling.pdf $(ANALYSIS)/cd_modeling.pdf
	ln -sf $(ROOT)/cluster_density/work/downstream/results/figures/tn5_chromatin_state_overlap.pdf $(ANALYSIS)/cd_chromatin_state_overlap.pdf
	ln -sf $(ROOT)/tn5/work/downstream/results/figures/tn5_tf_overlap.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/tn5/work/downstream/results/figures/tn5_sensitivity_prob_by_chipseq_overlap_and_peak_signal.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/tn5/work/downstream/results/figures/tss_coverage_different_methods.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/tn5/work/downstream/results/figures/tss_coverage_positional_stability.pdf $(ANALYSIS)/
	ln -sf $(ROOT)/tn5/work/downstream/results/figures/tss_enrichment_housekeeping_vs_all.pdf $(ANALYSIS)/
	# make addition plots
	cd $(ANALYSIS) && Rscript plot-nuclear-isolation-efficiency.R $(DATA)/nuclei_isolation/nuclei_isolation.csv
	# plot qPCR curves
	cat $(DATA)/qpcr/2017-08-10.tsv | awk 'NR>36' | perl -pe 's@2/3@0.66@g; s@1/2@0.5@g; s@1/5@0.2@g' > $(ANALYSIS)/2017-08-10.txt
	cat $(DATA)/qpcr/2017-08-14.tsv | awk 'NR>33' | perl -pe 's@2/3@0.66@g; s@1/2@0.5@g; s@1/5@0.2@g' > $(ANALYSIS)/2017-08-14.txt
	cd $(ANALYSIS) && Rscript plot-qPCR-curves.R && rm 2017-08-10.txt 2017-08-14.txt
	# plot PCR cycles for PCR-variable experiment
	cd $(ANALYSIS) && Rscript plot-pcr-cycles-vs-tn5.R $(ROOT)/cluster_density/sample_info/sample_info.txt
	cd $(ANALYSIS) && pdflatex main.tex && pdflatex main.tex 
