#!/usr/bin/env nextflow

IONICE = 'ionice -c2 -n7'

// Generic data
AUTOSOMAL_REFERENCES = ['hg19': (1..22).collect({it -> 'chr' + it}),
			'hg38': (1..22).collect({it -> 'chr' + it}),
			'rn5': (1..20).collect({it -> 'chr' + it}),
			'rn6': (1..20).collect({it -> 'chr' + it}),
			'mm9': (1..19).collect({it -> 'chr' + it}),
			'mm10': (1..19).collect({it -> 'chr' + it})
]

ORGANISMS = ['hg19': 'human', 
		'hg38': 'human',
		'rn5': 'rat',
		'rn6': 'rat',
		'mm9': 'mouse',
		'mm10': 'mouse']

libraries = params.libraries.keySet()

make_excluded_regions_arg = {
	genome ->
	return params.blacklist[genome].collect({'--excluded-region-file ' + it}).join(' ')
}


get_genome = {
	library ->
	params.libraries[library].genome
}


get_tss = {
	genome ->
	params.tss[genome]
}


get_organism = {
	genome ->
	ORGANISMS[genome]
}


get_chrom_sizes = {
	genome ->
	params.chrom_sizes[genome]
}


get_gtf = {
	genome ->
	params.gtf[genome]
}

get_peaks = {
	library ->
	params.atacseq.peaks + '/' + library + '_peaks.broadPeak.noblacklist'
}

get_bam = {
	library ->
	params.atacseq.pruned + '/' + library + '.pruned.bam'
}

get_md_bam = {
	library ->
	params.atacseq.md + '/' + library + '.md.bam'
}

get_project = {
	library ->
	params.libraries[library].project
}

get_runs = {
	library ->
	params.libraries[library].runs
}

get_class = {
	library ->
	params.libraries[library]['class']
}

is_bulk = {
	library ->
	params.libraries[library]['class'] == 'bulk'
}

// get flowcell info...
flowcell_in = []

for (library in libraries) {
	if (get_project(library) == 'PRJNA323617.PRJNA341508') {
		for (run in get_runs(library)) {
			flowcell_in << run
		}
	}
}

/*
process wait_2_seconds {

	maxForks 1

	input:
	val(srr) from Channel.from(flowcell_in)

	output:
	val(srr) into wait_2_seconds_out

	"""
	sleep 2
	""'

}

//flowcell_fastq_in = Channel.fromSRA(flowcell_in, apiKey: 'fba9b758a379c9f16486e5513169979f6407').map({it -> [it[0], it[1][0]]})
*/

process flowcell_fastq {

	publishDir "${params.results}/fastq"

	maxForks 2
	cpus 1
	errorStrategy 'retry'
	maxRetries 2

	input:
	val(srr) from Channel.from(flowcell_in)

	output:
	set val(srr), file("${srr}_1.fastq") into flowcell_fastq_out

	"""
	fasterq-dump ${srr} -S -e 1
	"""
	
}

process flowcell {

	publishDir "${params.results}/flowcell"
	time '4h'
	maxForks 3
	cpus 1
	errorStrategy 'retry'
	maxRetries 1

	input:
	set val(srr), file(fastq) from flowcell_fastq_out
	
	output:
	file("${srr}.flowcells.txt") into flowcell_out_chan

	"""
	cat $fastq | awk '(NR-1)%4==0' | perl -pe 's/ /\t/g' | cut -f2 | perl -pe 's/:.*//' | uniq | sort | uniq  > ${srr}.flowcells.txt
	"""

}




// run ataqv with the TSS file we want
ataqv_in = []

for (library in libraries) {
	ataqv_in << [library, file(get_md_bam(library)), file(get_md_bam(library) + '.bai'), file(get_peaks(library)), file(get_tss(get_genome(library)))]
}

process ataqv {

	publishDir "${params.results}/ataqv"
	time '4h'
	maxForks 40
	errorStrategy 'retry'
	maxRetries 2

	input:
	set val(library), file(md_bam), file(index), file(peaks), file(tss) from Channel.from(ataqv_in)

	output:
	set val(library), file("${library}.ataqv.json.gz"), file("${library}.ataqv.out") into ataqv_out

	"""
	ataqv --peak-file $peaks --name $library --metrics-file ${library}.ataqv.json.gz ${make_excluded_regions_arg(get_genome(library))} --tss-file $tss --ignore-read-groups ${get_organism(get_genome(library))} $md_bam > ${library}.ataqv.out
	"""

}

ataqv_out.into{ ataqv_out_1; ataqv_out_2 }

// make the ataqv sessions
// want to do one session per project
// ataqv_session_in_chan = ataqv_out_1.map({it -> [[get_project(it[0]), get_class(it[0])], [it[1], it[1]]]}).transpose().groupTuple()
ataqv_session_in_chan = ataqv_out_1.map({it -> [get_project(it[0]), it[1]]}).groupTuple()

process ataqv_session {

	publishDir "${params.results}/ataqv-sessions"
	errorStrategy 'retry'
	maxRetries 3
	time '4h'

	input:
	set val(project), file(ataqv) from ataqv_session_in_chan

	output:
	set val(project), file("${project}") into ataqv_metrics_in

	"""
	mkarv -r SRR891268 $project *.json.gz
	"""
}

process ataqv_metrics {
	publishDir "${params.results}/ataqv-metrics"
	
	input:
	set val(project), file(session) from ataqv_metrics_in

	output:
	set val(project), file("${project}.fragment_length_distributions.txt"), file("${project}.tss_coverage.txt"), file("${project}.chrom_counts.txt"), file("${project}.metrics.txt") into ataqv_metrics_out
	set val(project), file("${project}.metrics.txt") into bulk_session
	set file("${project}.fragment_length_distributions.txt"), file("${project}.tss_coverage.txt"), file("${project}.chrom_counts.txt"), file("${project}.metrics.txt") into split_metrics 

	"""
	extractFragmentLengths.py --files ${session}/data/*.json.gz | perl -pe 's@.*/data/(.*).json.gz@\$1@' > ${project}.fragment_length_distributions.txt
	extractTssCoverage.py --files ${session}/data/*.json.gz | perl -pe 's@.*/data/(.*).json.gz@\$1@' > ${project}.tss_coverage.txt
	extractChromosomeCounts.py --files ${session}/data/*.json.gz | perl -pe 's@.*/data/(.*).json.gz@\$1@' > ${project}.chrom_counts.txt
	extractAtaqvMetric.py --files ${session}/data/*.json.gz --metrics percent_hqaa percent_properly_paired_and_mapped percent_autosomal_duplicate short_mononucleosomal_ratio tss_enrichment duplicate_fraction_in_peaks duplicate_fraction_not_in_peaks peak_duplicate_ratio hqaa_overlapping_peaks_percent total_peak_territory total_reads percent_secondary percent_supplementary percent_duplicate mean_mapq median_mapq percent_unmapped percent_unmapped_mate percent_qcfailed percent_unpaired percent_mapq_0 percent_rf percent_ff percent_rr percent_mate_separate_chromosome percent_mate_too_distant percent_improperly_paired percent_autosomal percent_mitochondrial percent_mitochondrial_duplicate total_peaks fragment_length_distance hqaa median_fragment_length max_fraction_reads_from_single_autosome | perl -pe 's@.*/data/(.*).json.gz@\$1@' > ${project}.metrics.txt
	"""
}

split_metrics_in =  split_metrics.flatten().toSortedList()

classes = ['bulk', 'single-cell']
species = ['Homo sapiens', 'Mus musculus']

process split_metrics_by_class_and_species {

	input:
	file(ataqv_output) from split_metrics_in
	each cl from classes
	each sp from species

	output:
	set val(cl), val(sp), file("${cl}.${genome}.fragment_length_distributions.txt"), file("${cl}.${genome}.tss_coverage.txt"), file("${cl}.${genome}.chrom_counts.txt"), file("${cl}.${genome}.metrics.txt")
	set val(cl), file("${cl}.${genome}.metrics.txt") into split_metrics_bulk
	set val(cl), file("${cl}.${genome}.metrics.txt") into public_mfl_vs_tss_enrichment_in
	set val(cl), val(genome), file("${cl}.${genome}.fragment_length_distributions.txt") into example_fldd_calculation_in
	set val(cl), val(genome), file("${cl}.${genome}.metrics.txt"), file("${cl}.${genome}.chrom_counts.txt") into split_metrics_max_fraction_in
	set val(cl), val(genome), file("${cl}.${genome}.metrics.txt"), file("${cl}.${genome}.chrom_counts.txt") into split_metrics_max_fraction_fig_1_in
	set val(cl), file("${cl}.${genome}.fragment_length_distributions.txt"), file("${cl}.${genome}.tss_coverage.txt") into intrastudy_heterogeneity_in
	
	script:
	genome = (sp == 'Homo sapiens') ? 'hg19' : 'mm9'

	"""
	cat *.fragment_length_distributions.txt > flds.txt
	cat *.tss_coverage.txt > tss_cov.txt
	cat *.metrics.txt > metrics.txt
	cat *.chrom_counts.txt > chrom_counts.txt
	split-metrics-by-class-and-species.R ${params.sample_info} '$sp' '$cl' flds.txt ${cl}.${genome}.fragment_length_distributions.txt
	split-metrics-by-class-and-species.R ${params.sample_info} '$sp' '$cl' tss_cov.txt ${cl}.${genome}.tss_coverage.txt
	split-metrics-by-class-and-species.R ${params.sample_info} '$sp' '$cl' metrics.txt ${cl}.${genome}.metrics.txt
	split-metrics-by-class-and-species.R ${params.sample_info} '$sp' '$cl' chrom_counts.txt ${cl}.${genome}.chrom_counts.txt
	"""

}

ataqv_metrics_out.into{ ataqv_metrics_out_2; ataqv_metrics_out_3 }
bulk_session = split_metrics_bulk.filter({it[0].toString() == 'bulk'}).map({it -> it[1]}).toSortedList()

process explore_ataqv_metrics {

	publishDir "${params.results}/figures"

	input:
	file(metrics) from bulk_session

	output:
	file('all-corr-heatmap.pdf')

	"""
	cat ${metrics.join(' ')} > metrics.txt
	explore-ataqv-metrics.py --ataqv-metrics metrics.txt --sample-info ${params.sample_info}
	"""

}

// run PCA
// first, create master peaks for each project
// then gather the peak counts
// then run pca
master_peaks_in = []
for (library in libraries) {
	// want project, genome, peaks
	if (is_bulk(library)) {
		master_peaks_in << [get_project(library), get_genome(library), file(get_peaks(library))]
	}
}

master_peaks_in_chan = Channel.from(master_peaks_in).groupTuple(by: [0, 1]).filter({it -> it[2].size() > 1})

process master_peaks {

	publishDir "${params.results}/master-peaks"
	
	input:
	set val(project), val(genome), file(peaks) from master_peaks_in_chan

	output:
	set val(project), val(genome), file("${project}-${genome}.bed") into master_peaks_out_chan

	"""
	make-master-peaks.py --fdr 0.05 --min-samples 2 --peak-files $peaks | sort -k1V,1 -k2n,2 > ${project}-${genome}.bed
	"""

}

peak_counts_in = []
for (library in libraries) {
	if (is_bulk(library)) {
		peak_counts_in << [get_project(library), get_genome(library), library, file(get_bam(library))]
	}
}

peak_counts_in_chan = Channel.from(peak_counts_in).combine(master_peaks_out_chan, by: [0, 1])

process peak_counts {

	publishDir "${params.results}/peak-counts"
	maxForks 5
	validExitStatus 0,141
	

	input:
	set val(project), val(genome), val(library), file(bam), file(master_peaks) from peak_counts_in_chan

	output:
	set val(project), val(genome), file("${library}.bed") into peak_counts_out_chan

	"""
	coverageBed -counts -sorted -a $master_peaks -b $bam > ${library}.bed
	"""

}

pca_in_chan = peak_counts_out_chan.groupTuple(by: [0, 1])

process pca {
	publishDir "${params.results}/pca"
	
	input:
	set val(project), val(genome), file(counts) from pca_in_chan

	output:
	set val(project), val(genome), file("${project}-${genome}.pca_scores.txt") into pca_out_chan

	"""
	pca.R --counts-dir . --log2 --pseudocount 1 --prefix ${project}-${genome}.pca_
	"""

}

// plot top PC vs QC metrics
pc1_vs_ataqv_metrics_in = pca_out_chan.map({it -> [it[0].toString(), it[1], it[2]]}).combine(ataqv_metrics_out_3.map({it -> [it[0].toString(), it[4]]}), by: 0)

process pc1_vs_ataqv_metrics {

	publishDir "${params.results}/pc1-vs-ataqv-metrics"

	input:
	set val(project), val(genome), file(pc_scores), file(metrics) from pc1_vs_ataqv_metrics_in

	output:
	set file("${project}-${genome}.cor_barplot.pdf"), file("${project}-${genome}.all_scatter.pdf"), file("${project}-${genome}.tss_enrichment_scatter.pdf")

	"""
	plot-pc1-vs-qc-metrics.R $pc_scores $metrics ${project}-${genome}. $project ${params.project_colors}
	"""

}


// make heatmaps
// need GTF file
heatmaps_in = []
for (library in libraries) {
	if (is_bulk(library) && get_genome(library) == 'hg19') {
		heatmaps_in << [library, get_genome(library), file(get_md_bam(library)), file(get_md_bam(library) + '.bai')]
	}
}

heatmaps_in_chan = Channel.from(heatmaps_in)

genes = ['GAPDH', 'VCP']
process heatmap_windows {

	input:
	each gene from genes

	output:
	set val(gene), file('windows.bed') into heatmap_windows_out

	"""
	zcat ${get_gtf('hg19')} | grep -i 'gene_name \"$gene\"' | cut -f1,4,5 | sort -k1,1 | bedtools groupby -g 1 -c 2,3 -o min,max > coords.bed
	cat coords.bed | bedtools slop -i stdin -g ${get_chrom_sizes('hg19')} -b 15000 | bedtools makewindows -b stdin -w 20 > windows.bed
	"""

}

heatmaps_in_chan = heatmaps_in_chan.combine(heatmap_windows_out)

process heatmaps {

	publishDir "${params.results}/heatmaps"
	errorStrategy 'retry'
	maxRetries 1
	memory { 20.GB * task.attempt }
	maxForks 5

	input:
	set val(library), val(genome), file(md_bam), file(index), val(gene), file(windows) from heatmaps_in_chan

	output:
	file("${library}.${gene}_coverage.bed") into heatmaps_out_chan

	"""
	samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 30 -L $windows $md_bam | bedtools bamtobed | perl -lane 'if(\$F[5]=="-"){\$F[1]=\$F[2]-1};\$F[2]=\$F[1]+1;print join("\\t", @F)' | sort -k1,1 -k2n,2 > ${library}.bed
	bedtools coverage -counts -a windows.bed -b ${library}.bed > ${library}.${gene}_coverage.bed
	"""

}

prep_gtf_in = Channel.from(file(get_gtf('hg19')))

process prep_gtf_for_python {

	memory '20 GB'

	input:
	file('my_gtf.gtf.gz') from prep_gtf_in

	output:
	file('gtf.db') into gtf_db

	"""
	#!/usr/bin/env python
	import gffutils
	GTF = 'my_gtf.gtf.gz'
	db = gffutils.create_db(GTF, dbfn='gtf.db', force=True, keep_order=False, merge_strategy='merge', sort_attribute_values=False, disable_infer_transcripts=True, disable_infer_genes=True)
	"""

}

plot_heatmaps_coverage_chan = heatmaps_out_chan.toSortedList()
plot_heatmaps_ataqv_chan = ataqv_out_2.map({it -> it[1]}).toSortedList()

process plot_heatmaps {
	publishDir "${params.results}/figures"
	memory '75 GB'
	
	input:
	file("*") from plot_heatmaps_coverage_chan
	file("*") from plot_heatmaps_ataqv_chan
	file('gtf.db') from gtf_db

	output:
	set file("GAPDH-heatmap.pdf"), file("VCP-heatmap.pdf")

	"""
	public-heatmap.py --sample-info ${params.sample_info} --gtf-db gtf.db --heatmap-files *.GAPDH_coverage.bed --chrom chr12 --start 6642500 --end 6645000 --ataqv-files *.json.gz --out GAPDH-heatmap.pdf --width 5.5 --height 4 --project-colors ${params.project_colors}
	public-heatmap.py --sample-info ${params.sample_info} --gtf-db gtf.db --heatmap-files *.VCP_coverage.bed --chrom chr9 --start 35071000 --end 35074000 --ataqv-files *.json.gz --out VCP-heatmap.pdf --width 5.5 --height 4 --heatmap-and-gene-model-only --project-colors ${params.project_colors}
	"""
}

// Run Su et al differential peak analysis
// First, get the master peaks
// then get the master peak counts
// then run the analysis (make the covariate files and run the R script)

SU_ET_AL_DIFFERENTIAL_PEAK_LIBRARIES = ["SRX1804207", "SRX1804208", "SRX1804209", "SRX1804210", "SRX1804211", "SRX1804212", "SRX1804213", "SRX1804214"]
su_e0_e1_peaks = []
su_e0_e1_bams = []
for (library in SU_ET_AL_DIFFERENTIAL_PEAK_LIBRARIES) {
	su_e0_e1_peaks << file(get_peaks(library))
	su_e0_e1_bams << [library, file(get_bam(library))]
}

process su_master_peaks {

	publishDir "${params.results}/su/master-peaks"
	
	input:
	file(peaks) from Channel.from(su_e0_e1_peaks).toSortedList()

	output:
	file('master-peaks.bed') into su_master_peaks_out_chan

	"""
	make-master-peaks.py --fdr 0.01 --min-samples 2 --peak-files ${peaks.join(' ')} | sort -k1V,1 -k2n,2 > master-peaks.bed
	"""

}

su_master_peak_counts_in_chan = Channel.from(su_e0_e1_bams).combine(su_master_peaks_out_chan)

process su_master_peak_counts {

	publishDir "${params.results}/su/master-peak-counts"
	
	input:
	set val(library), file(bam), file(master_peaks) from su_master_peak_counts_in_chan

	output:
	file("${library}.counts.bed") into su_master_peak_counts_out_chan

	"""
	coverageBed -counts -sorted -a $master_peaks -b $bam > ${library}.counts.bed
	"""
	
}

counts_to_library = {
        (full, library) = (it =~ /.*\/(.*).counts.bed/)[0]
        return library
}

//su_master_peak_counts_for_differential_analysis = su_master_peak_counts_out_chan.filter({it -> SU_ET_AL_DIFFERENTIAL_PEAK_LIBRARIES.indexOf(counts_to_library(it)) != -1}).toSortedList()
//differential_peak_analysis_in = ataqv_metrics_out_2.filter({it[0].toString() == 'PRJNA323617.PRJNA341508'}).map({it -> it[4]}).combine(su_master_peak_counts_for_differential_analysis).combine(flowcell_out_chan.collect())
su_metrics = ataqv_metrics_out_2.filter({it[0].toString() == 'PRJNA323617.PRJNA341508'}).map({it -> it[4]})
//su_master_peak_counts_out_chan.toSortedList().subscribe({println it})
//flowcell_out_chan.toSortedList().subscribe({println it})

process su_differential_peak_analysis {

	publishDir "${params.results}/su/differential-peak-analysis"

	input:
	file("*") from su_master_peak_counts_out_chan.toSortedList() 
	file(metrics) from su_metrics
	file("*") from flowcell_out_chan.toSortedList()

	output:
	set file('su_et_al_median_fragment_length_vs_batch.pdf'), file('su-et-al-qqplot.pdf'), file('median-fragment-length-covariate-p-value-distribution.pdf')

	"""
	make-su-et-al-sample-matrix.R ${params.sample_info} $metrics
	deseq.R --counts . --out no-covariate.results.txt --covariates covariates.no_median_fragment_length.txt
	deseq.R --counts . --out with-covariate.results.txt --covariates covariates.with_median_fragment_length.txt
	su-et-al-qqplot.py no-covariate.results.txt with-covariate.results.txt
	plot-su-et-al-median-fragment-length.R ${params.sample_info} $metrics
	"""

}



plot_max_fraction_in = split_metrics_max_fraction_in.filter({it[0] == 'single-cell' || it[0] == 'bulk'})
process plot_max_fraction_from_single_autosome {

	publishDir "${params.results}/max-fraction-from-single-autosome"

	input:
	set val(category), val(genome), file(metrics), file(chrom_counts) from plot_max_fraction_in

	output:
	file("${project}*")

	script:
	project = "${category}.${genome}"

	"""
	max-fraction.R ${params.sample_info} $metrics $chrom_counts ${project}. $category ${params.project_colors}
	"""

}

process fig1_max_fraction_from_single_autosome {

	publishDir "${params.results}/figures"

	input:
	set val(category), val(genome), file(metrics), file(chrom_counts) from split_metrics_max_fraction_fig_1_in

	output:
	set file("fig1-hqaa-vs-max-fraction-color-by-chrom.pdf"), file('nonoutliers.txt')
	file('nonoutliers.txt') into fig1_chromosome_coverage

	when:
	genome == 'hg19' && category == 'single-cell'

	"""
	Rscript ${params.src}/max-fraction-figure-1.R ${params.sample_info} $metrics $chrom_counts
	"""

}

make_chrom_coverage_in = []
for (library in libraries) {
	if (!is_bulk(library)) {
		make_chrom_coverage_in << [library, file(get_bam(library))]
	}
}

process make_chromosome_coverage {

	memory '5 GB'
	time '1h'
	errorStrategy 'retry'
	maxRetries 2

	publishDir "${params.results}/chromosome-coverage", mode: 'rellink', overwrite: true

	input:
	set val(library), file(bam) from Channel.from(make_chrom_coverage_in)

	output:
	file("${name}.chromosome-windows.bed") into make_chrom_coverage_out

	script:
	window_size = 2000000
	sliding_window_size = window_size / 2
	name = "${library}.${window_size}"

	"""
	bedtools makewindows -g ${get_chrom_sizes(get_genome(library))} -w $window_size -s $sliding_window_size > windows.bed
	bedtools intersect -a $bam -b ${params.blacklist[get_genome(library)].join(' ')} -v -bed | bedtools coverage -counts -a windows.bed -b stdin > ${name}.chromosome-windows.bed
	plot-coverage-across-chromosomes.R ${name}.chromosome-windows.bed ${name}.coverage-across-chromosomes.pdf 
	"""

}

process fig1_chrom_coverage_plot {

	publishDir "${params.results}/figures"

	input:
	file(covs) from make_chrom_coverage_out.toSortedList()
	file(nonoutliers) from fig1_chromosome_coverage

	output:
	file('compare-chromosome-coverages.pdf')

	"""
	compare-chromosome-coverages.R $nonoutliers ${covs.join(' ')}
	"""

}

// plot public metadata
process plot_public_metadata {

	publishDir "${params.results}/figures"

	output:
	file('public_meta.pdf')

	"""
	plot-public-metadata.R ${params.sample_info} ${params.project_colors}
	"""

}

process example_fldd_calculation {

	publishDir "${params.results}/figures"

	input:
	set val(cl), val(genome), file(flds) from example_fldd_calculation_in

	output:
	set file('example_fld_calculation_cumulative.pdf'), file('example_fld_calculation_density.pdf')

	when:
	cl == 'bulk' && genome == 'hg19'

	"""
	example-fldd-calculation.R $flds
	"""

}

bulk_metrics_public_mfl_vs_tss_enrichment_in = public_mfl_vs_tss_enrichment_in.filter({it[0].toString() == 'bulk'}).map({it -> it[1]}).flatten().toSortedList()
process public_mfl_vs_tss_enrichment {

	publishDir "${params.results}/figures"

	input:
	file(metrics) from bulk_metrics_public_mfl_vs_tss_enrichment_in

	output:
	file('public_median_fragment_length_vs_tss_enrichment.pdf')

	"""
	cat ${metrics.join(' ')} > metrics.txt
	plot-public-mfl-vs-tss-enrichment.R ${params.sample_info} metrics.txt ${params.project_colors}
	"""

}

intrastudy_heterogeneity_in_chan = intrastudy_heterogeneity_in.filter({it[0].toString() == 'bulk'}).map({it -> [it[1], it[2]]}).flatten().toSortedList()

process intrastudy_heterogeneity {

	publishDir "${params.results}/figures"

	input:
	file("*") from intrastudy_heterogeneity_in_chan

	output:
	set file('public_intrastudy_heterogeneity_fragment_length_distributions.pdf'), file('public_intrastudy_heterogeneity_tss_enrichment.pdf')

	"""
	cat *.fragment_length_distributions.txt > flds.txt
	cat *.tss_coverage.txt > tss_cov.txt
	plot-intrastudy-heterogeneity.R ${params.sample_info} flds.txt tss_cov.txt
	"""

}
