#!/usr/bin/env nextflow

// TODO: capture numeric values we want for the manuscript

IONICE = 'ionice -c2 -n7'

// Generic data and functions
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
low_cd_libraries = params.low_cd_libraries.keySet()

make_excluded_regions_arg = {
	genome ->
	return params.blacklist[genome].collect({'--excluded-region-file ' + it}).join(' ')
}

get_cluster_density = {
	library ->
	if (params.libraries.containsKey(library)) {
		return 'high'
	} else {
		return 'low'
	}
}

get_genome = {
	library ->
	if (get_cluster_density(library) == 'high') {
		params.libraries[library].genome
	} else {
		params.low_cd_libraries[library].genome
	}
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

get_ataqv = {
	library ->
	params.atacseq.ataqv + '/' + library + '.ataqv.json.gz'
}

get_peaks = {
	library ->
	params.atacseq.peaks + '/' + library + '_peaks.broadPeak.noblacklist'
}

get_bedgraph = {
	library ->
	params.atacseq.peaks + '/' + library + '_treat_pileup.bdg'
}

get_bam = {
	library ->
	params.atacseq.pruned + '/' + library + '.pruned.bam'
}

get_md_bam = {
	library ->
	params.atacseq.md + '/' + library + '.md.bam'
}

get_tn5 = {
	library ->
	if (get_cluster_density(library) == 'high') {
		params.libraries[library].tn5
	} else {
		params.low_cd_libraries[library].tn5
	}
}

get_replicate = {
	library ->
	if (get_cluster_density(library) == 'high') {
		params.libraries[library].replicate
	} else {
		params.low_cd_libraries[library].replicate
	}
}


// use low and high for master peaks and for ataqv
// use just high for everything else

// Processing starts here
// First, create master peaks
master_peaks_in = []

for (library in libraries) {
	if (get_tn5(library) == '1X') {
		master_peaks_in << [file(get_peaks(library))]
	}
}

for (library in low_cd_libraries) {
	if (get_tn5(library) == '1X') {
		master_peaks_in << [file(get_peaks(library))]
	}
}


master_peaks_in_chan = Channel.from(master_peaks_in).collect()

process master_peaks {

	publishDir "${params.results}/master-peaks"

	input:
	file("peaks-*.bed") from master_peaks_in_chan

	output:
	file('master_peaks.bed') into master_peaks_out_chan
	file('master_peaks.bed') into master_peaks_tf_overlap_in

	"""
	cat peaks-*.bed | awk '\$9>=2' | sort -k1,1 -k2n,2 | bedtools merge -i stdin > merged.bed
	cat peaks-*.bed | awk '\$9>=2' | sort -k1,1 -k2n,2 | bedtools intersect -a merged.bed -b stdin -wa -wb | sort -k1,1 -k2n,2 | cut -f1-3,7 | perl -pe 's/_peak.*//; s/___.*//' | sort -k1,1 -k2n,2 | bedtools groupby -g 1,2,3 -o count_distinct -c 4 | awk '\$4>=2' | cut -f1-3 > master_peaks.bed
	"""

}

sort_pruned_in = []
normalized_bigwig_in = []
tss_enrichment_in = []

for (library in libraries) {
	sort_pruned_in << [library, file(get_bam(library))]
	normalized_bigwig_in << [library, file(get_bedgraph(library))]
	tss_enrichment_in << [library, file(get_bam(library))]
}

process sort_pruned {

	maxForks 10

	input:
	set val(library), file(bam) from Channel.from(sort_pruned_in)

	output:
	set val(library), file("${library}.sorted.pruned.bam"), file("${library}.sorted.pruned.bam.bai") into sort_pruned_out_chan

	"""
	${IONICE} samtools sort -m 3G -O bam -T ${library}.sort -o ${library}.sorted.pruned.bam $bam
	${IONICE} samtools index ${library}.sorted.pruned.bam
	"""

}

master_peaks_out_chan.into{master_peak_counts_in_chan; ataqv_master_peaks_in}
master_peak_counts_in_chan = sort_pruned_out_chan.combine(master_peak_counts_in_chan)

process master_peak_counts {
	
	maxForks 10
	publishDir "${params.results}/master_peak_counts"

	input:
	set val(library), file(bam), file(index), file(master_peaks) from master_peak_counts_in_chan

	output:
	set val(library), file("${library}.counts.bed") into master_peak_counts_out_chan
	file("${library}.counts.bed") into tn5_sensitive_peaks_counts_in
	file("${library}.counts.bed") into pca_in
	file("${library}.counts.bed") into tn5_sensitivity_vs_peak_signal_in

	"""
	${IONICE} coverageBed -counts -sorted -a $master_peaks -b $bam > ${library}.counts.bed	
	"""
}


process normalized_bigwig {

	publishDir "${params.results}/bigwig"
	maxForks 2
	errorStrategy 'ignore'

	input:
	set val(library), file(bedgraph) from Channel.from(normalized_bigwig_in)

	output:
	file("${library}.normalized.bw")

	"""
	normalize_bedgraph.py --to-number-reads 10000000 $bedgraph | LC_COLLATE=C sort -k1,1 -k2,2n > normalized.bdg
	bedGraphToBigWig normalized.bdg ${get_chrom_sizes(get_genome(library))} ${library}.normalized.bw
	"""

}

process tss_enrichment {
	
	maxForks 10
	publishDir "${params.results}/tss-enrichment"
	memory '11 GB'

	input:
	set val(library), file(bam) from Channel.from(tss_enrichment_in)

	output:
	set file("${library}.encode.txt"), file("${library}.cutsite.txt"), file("${library}.ataqv.txt") into tss_enrichment_out

	"""
	${IONICE} samtools sort -m 8G -n -T ${library}.sort -O bam -o name_sorted.bam $bam
	zcat ${get_tss('hg19')} | grep -P '^chr\\d+\\t' > tss.bed
	tss_enrichment.py --method ENCODE --read-length 36 name_sorted.bam tss.bed ${get_chrom_sizes('hg19')} > ${library}.encode.txt
	tss_enrichment.py --method ENCODE --read-length 2 name_sorted.bam tss.bed ${get_chrom_sizes('hg19')} > ${library}.cutsite.txt
	tss_enrichment.py --method ataqv name_sorted.bam tss.bed ${get_chrom_sizes('hg19')} > ${library}.ataqv.txt
	"""
}

process plot_compare_tss_enrichment_methods {

	publishDir "${params.results}/figures"
	
	input:
	file(covs) from tss_enrichment_out.flatten().toSortedList()

	output:
	set file('tss_coverage_different_methods.pdf'), file('tss_coverage_positional_stability.pdf')

	"""
	tss-enrichment-method-comparison.R ${params.sample_info} ${covs.join(' ')}
	"""

}

// Now, make the ataqv sessions
// Want to do one with individual peak calls
// Want to do one with master peaks as well
ataqv_in = []

for (library in libraries) {
	ataqv_in << file(get_ataqv(library))
}
for (library in low_cd_libraries) {
	ataqv_in << file(get_ataqv(library))
}

ataqv_session_in_chan = Channel.from(ataqv_in).collect()

process ataqv_session {

	publishDir "${params.results}/ataqv-sessions"
	errorStrategy 'retry'
	maxRetries 3
	time '4h'

	input:
	file(ataqv) from ataqv_session_in_chan

	output:
	file("${params.name}")

	"""
	mkarv -r SRR891268 ${params.name} ${ataqv.join(' ')}
	"""
}

// run ataqv using master peaks
ataqv_in = []
for (library in libraries) {
	ataqv_in << [library, file(get_md_bam(library)), file(get_md_bam(library) + '.bai')]
}
for (library in low_cd_libraries) {
	ataqv_in << [library, file(get_md_bam(library)), file(get_md_bam(library) + '.bai')]
}


ataqv_in_chan = Channel.from(ataqv_in).combine(ataqv_master_peaks_in)

process ataqv_master_peaks {

	publishDir "${params.results}/ataqv_master_peaks"
	time '2h'

	input:
	set val(library), file(bam), file(index), file(peaks) from ataqv_in_chan

	output:
	file("${library}.ataqv.json.gz") into ataqv_master_peaks_out_chan
	set val(library), file("${library}.ataqv.json.gz") into ataqv_informative_names_in
	file("${library}.ataqv.json.gz") into plot_heatmaps_ataqv_in

	"""
	ataqv --peak-file $peaks --name $library --metrics-file ${library}.ataqv.json.gz ${make_excluded_regions_arg(get_genome(library))} --tss-file ${get_tss(get_genome(library))} --ignore-read-groups human $bam > ${library}.ataqv.out
	"""

}

process ataqv_informative_names {

	publishDir "${params.results}/ataqv_informative_names"
	time '1h'

	input:
	set val(library), file('old_ataqv.json.gz') from ataqv_informative_names_in

	output:
	file("${library}.ataqv.json.gz") into ataqv_informative_names_out

	script:
	new_name = get_tn5(library) + ' Tn5, ' + get_replicate(library) + ', ' + get_cluster_density(library) + ' cluster density'	

	"""
	rename-ataqv.py old_ataqv.json.gz '${new_name}' | gzip -c > ${library}.ataqv.json.gz
	"""

}

process ataqv_informative_names_session {

	publishDir "${params.results}/ataqv-sessions"

	input:
	file(x) from ataqv_informative_names_out.toSortedList()

	output:
	file("ataqv-tn5-series-pcr-controlled-informative-names")

	"""
	mkarv -r SRR891268 ataqv-tn5-series-pcr-controlled-informative-names ${x.join(' ')}
	"""

}


process ataqv_master_peaks_session {

	publishDir "${params.results}/ataqv-sessions"
	
	input:
	file(ataqv) from ataqv_master_peaks_out_chan.toSortedList()

	output:
	file("${params.name}-master-peaks") into ataqv_master_peaks_session_out

	"""
	mkarv -r SRR891268 ${params.name}-master-peaks ${ataqv.join(' ')}
	"""
}

process ataqv_metrics {

	publishDir "${params.results}/ataqv-metrics"

	input:
	file(session) from ataqv_master_peaks_session_out

	output:
        set file("fragment_length_distributions.txt"), file("tss_coverage.txt"), file("metrics.txt") into ataqv_metrics_out
        set file("high-cd-fragment_length_distributions.txt"), file("high-cd-tss_coverage.txt") into plot_fld_and_tss_enrichment_in
        file("high-cd-metrics.txt") into plot_ataqv_metrics_in
        file("high-cd-metrics.txt") into tn5_sensitive_peaks_metrics_in
        file("metrics.txt") into compare_high_and_low_in

	"""
	extractFragmentLengths.py --files ${session}/data/*.json.gz | perl -pe 's@.*/data/(.*).json.gz@\$1@' > fragment_length_distributions.txt
        extractTssCoverage.py --files ${session}/data/*.json.gz | perl -pe 's@.*/data/(.*).json.gz@\$1@' > tss_coverage.txt
        extractAtaqvMetric.py --files ${session}/data/*.json.gz --metrics percent_hqaa percent_properly_paired_and_mapped percent_autosomal_duplicate short_mononucleosomal_ratio tss_enrichment duplicate_fraction_in_peaks duplicate_fraction_not_in_peaks peak_duplicate_ratio hqaa_overlapping_peaks_percent total_peak_territory total_reads percent_secondary percent_supplementary percent_duplicate mean_mapq median_mapq percent_unmapped percent_unmapped_mate percent_qcfailed percent_unpaired percent_mapq_0 percent_rf percent_ff percent_rr percent_mate_separate_chromosome percent_mate_too_distant percent_improperly_paired percent_autosomal percent_mitochondrial percent_mitochondrial_duplicate total_peaks fragment_length_distance hqaa median_fragment_length | perl -pe 's@.*/data/(.*).json.gz@\$1@' > metrics.txt
	grep -v '___1830' metrics.txt > high-cd-metrics.txt
	grep -v '___1830' fragment_length_distributions.txt > high-cd-fragment_length_distributions.txt
	grep -v '___1830' tss_coverage.txt > high-cd-tss_coverage.txt
	"""

}

process compare_high_and_low {
	
	publishDir "${params.results}/figures"

	input:
	file(metrics) from compare_high_and_low_in

	output:
	file("*.pdf")

	"""
	plot-high-vs-low.R ${params.sample_info_high_and_low} $metrics
	"""

}


process plot_ataqv_metrics {

	publishDir "${params.results}/figures"

	input:
	file(metrics) from plot_ataqv_metrics_in

	output:
	set file("tn5_tss_enrichment_all_six.pdf"), file("tn5_median_fragment_length_all_six.pdf"), file("tn5_hqaa_overlapping_peaks_percent.pdf"), file("tn5_mitochondrial_percent.pdf"), file("tn5_proportion_hqaa.pdf"), file("tn5_proportion_duplicate.pdf")

	"""
	plot-metric-vs-tn5.R --sample_info ${params.sample_info} --metrics $metrics --plot_metric tss_enrichment --ylab 'TSS enrichment' --out tn5_tss_enrichment_all_six.pdf
	plot-metric-vs-tn5.R --sample_info ${params.sample_info} --metrics $metrics --plot_metric median_fragment_length --ylab 'Median fragment length (bp)' --out tn5_median_fragment_length_all_six.pdf
	plot-metric-vs-tn5.R --sample_info ${params.sample_info} --metrics $metrics --plot_metric hqaa_overlapping_peaks_percent --ylab '% reads overlapping peaks' --out tn5_hqaa_overlapping_peaks_percent.pdf
	plot-metric-vs-tn5.R --sample_info ${params.sample_info} --metrics $metrics --plot_metric percent_mitochondrial --ylab '% mitochondrial reads' --out tn5_mitochondrial_percent.pdf
	plot-metric-vs-tn5.R --sample_info ${params.sample_info} --metrics $metrics --plot_metric percent_hqaa --ylab '% reads remaining after\nfiltering' --out tn5_proportion_hqaa.pdf
	plot-metric-vs-tn5.R --sample_info ${params.sample_info} --metrics $metrics --plot_metric percent_duplicate --ylab '% duplicate reads' --out tn5_proportion_duplicate.pdf
	"""

}


process plot_fld_and_tss_enrichment {

	publishDir "${params.results}/figures"

	input:
	set file(flds), file(tss_coverage) from plot_fld_and_tss_enrichment_in

	output:
	set file('tn5_fragment_length_distributions.pdf'), file('tn5_tss_enrichment.pdf')

	"""
	fig2-fld-and-tss-enrichment.R --sample_info ${params.sample_info} --fld $flds --tss_coverage $tss_coverage
	"""

}

// make heatmaps
heatmaps_in = []
for (library in libraries) {
	heatmaps_in << [library, get_genome(library), file(get_md_bam(library)), file(get_md_bam(library) + '.bai')]
}

genes = ['GAPDH']

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

heatmaps_in_chan = Channel.from(heatmaps_in).combine(heatmap_windows_out)

process heatmaps {

	publishDir "${params.results}/heatmaps"

	input:
	set val(library), val(genome), file(md_bam), file(index), val(gene), file(windows) from heatmaps_in_chan

	output:
	file("${library}.${gene}_coverage.bed") into heatmaps_out_chan


	"""
	samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 30 -L $windows $md_bam | bedtools bamtobed | perl -lane 'if(\$F[5]=="-"){\$F[1]=\$F[2]-1};\$F[2]=\$F[1]+1;print join("\\t", @F)' | sort -k1,1 -k2n,2 > ${library}.bed
	bedtools coverage -counts -a $windows -b ${library}.bed > ${library}.${gene}_coverage.bed
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

process plot_heatmaps {

	publishDir "${params.results}/figures"

	input:
	file(coverages) from heatmaps_out_chan.toSortedList()
	file(ataqv) from plot_heatmaps_ataqv_in.toSortedList()
	file(gtf) from gtf_db

	output:
	file("tn5-GAPDH-heatmap.pdf")

	"""
	tn5-heatmap.py --sample-info ${params.sample_info} --gtf-db $gtf --heatmap-files ${coverages.join(' ')} --chrom chr12 --start 6642500 --end 6645000 --ataqv-files ${ataqv.join(' ')} --out tn5-GAPDH-heatmap.pdf --width 5.5 --height 4
	"""

}

subsample_in = []
pairwise_jaccard_in = []

for (library in libraries) {
	subsample_in << file(get_bam(library))
	pairwise_jaccard_in << file(get_peaks(library))
}


process subsample {

	publishDir "${params.results}/subsample/bams"

	input:
	file(bams) from Channel.from(subsample_in).toSortedList()

	output:
	file("*.pruned.subsampled.bam") into subsample_peaks_in mode flatten

	"""
	subsampleBams.py --same-depth *.bam
	"""
}


process subsample_peaks {
	
	publishDir "${params.results}/subsample/peaks"

	input:
	file(bam) from subsample_peaks_in

	output:
	set val(library), file("${library}_peaks.broadPeak.noblacklist") into subsample_peaks_out
	set val('subsampled'), file("${library}_peaks.broadPeak.noblacklist") into subsample_plot_number_of_peaks_in

	script:
	library = bam.getName().replaceAll('.pruned.subsampled.bam', '')
	
	"""
	macs2 callpeak -t $bam --outdir . -f BAM -n $library -g hs --nomodel --shift -100 --seed 762873 --extsize 200 -B --broad --keep-dup all
        bedtools intersect -a ${library}_peaks.broadPeak -b ${params.blacklist['hg19'].join(' ')} -v > ${library}_peaks.broadPeak.noblacklist
	"""
}

subsample_pairwise_jaccard_in = subsample_peaks_out.map({it -> it[1]}).toSortedList()

process subsample_pairwise_jaccard {

	publishDir "${params.results}/subsample/figures"

	input:
	file(peaks) from subsample_pairwise_jaccard_in

	output:
	set file('binary-jaccard-indices.txt'), file('*.pdf')

	"""
	peak-binary-jaccard.py binary-jaccard-indices.txt ${peaks.join(' ')}
	plot-pairwise-jaccard.py ${params.sample_info} binary-jaccard-indices.txt pairwise-jaccard.subsample.
	"""
}


process pairwise_jaccard {

	publishDir "${params.results}/figures"

	input:
	file(peaks) from Channel.from(pairwise_jaccard_in).toSortedList() 

	output:
	set file('binary-jaccard-indices.txt'), file('*.pdf')

	"""
	peak-binary-jaccard.py binary-jaccard-indices.txt ${peaks.join(' ')}
	plot-pairwise-jaccard.py ${params.sample_info} binary-jaccard-indices.txt pairwise-jaccard.
	"""
}


nonsubsampled_peaks_in = []

for (library in libraries) {
	nonsubsampled_peaks_in << ["full", file(get_peaks(library))]
}

plot_number_of_peaks_in = subsample_plot_number_of_peaks_in.mix(Channel.from(nonsubsampled_peaks_in)).groupTuple(by: 0)

process plot_number_of_peaks {

	publishDir "${params.results}/figures"

	input:
	set val(g), file(peaks) from plot_number_of_peaks_in

	output:
	file("number-of-peaks-${g}.pdf")

	"""
	plot-number-of-peaks.py --out number-of-peaks-${g}.pdf --sample-info ${params.sample_info} --peaks ${peaks.join(' ')}
	"""

}


process chromhmm {

	input:
	file(x) from Channel.fromPath(params.chromhmm_glob)

	output:
	file("${state}.bed") into chromhmm_out

	script:
	state_match = x.toString() =~ /GM12878.(.*).bed/
	state = state_match[0][1]

	"""
	cat $x | perl -pe 's/\$/\t${state}/' > ${state}.bed	
	"""
}

concat_chromhmm_in = chromhmm_out.toSortedList() 


process concat_chromhmm {

	input:
	file(x) from concat_chromhmm_in

	output:
	file("chromhmm.bed") into concat_chromhmm_out

	"""
	cat ${x.join(' ')} | sort -k1,1 -k2n,2 > chromhmm.bed
	"""
}


blacklist_filter_in = []

for (library in libraries) {
	blacklist_filter_in << [library, file(get_bam(library))]
}

process blacklist_filter_reads {

	memory '5 GB'
	maxForks 10

	input:
	set val(library), file(bam) from Channel.from(blacklist_filter_in)

	output:
	set val(library), file("${library}.reads.bed") into read_tf_overlaps_in
	set val(library), file("${library}.reads.bed") into read_chromhmm_overlaps_in
	file("${library}.total_reads.txt") into blacklist_filtered_reads_out

	"""
	${IONICE} bedtools intersect -a $bam -b ${params.blacklist['hg19'].join(' ')} -v -bed > ${library}.reads.bed
	wc -l ${library}.reads.bed | perl -pe 's/^\\s+//; s/\\s/\\t/' > ${library}.total_reads.txt
	"""

}


blacklist_filtered_read_counts_in = blacklist_filtered_reads_out.toSortedList()

process blacklist_filtered_read_counts {

	publishDir "${params.results}/read-tf-overlaps"

	input:
	file(x) from blacklist_filtered_read_counts_in

	output:
	file("blacklist-filtered-read-counts.txt") into blacklist_filtered_read_counts_out

	"""
	cat ${x.join(' ')} | perl -pe 's/.reads.bed\$//' > blacklist-filtered-read-counts.txt
	"""

}


process read_chromhmm_overlaps {

	publishDir "${params.results}/read_chromhmm_overlaps"
	memory '10 GB'
	maxForks 10

	input:
	set val(library), file(bed), file(chromhmm) from read_chromhmm_overlaps_in.combine(concat_chromhmm_out)

	output:
	file("${library}.counts.bed") into chromhmm_overlaps_out

	"""
	bedtools intersect -wao -a $bed -b $chromhmm | cut -f4,16,17 | sort -k1,1 -k2,2 | bedtools groupby -g 1,2 -c 3 -o sum | sort -k1,1 -k3n,3 | bedtools groupby -g 1 -c 2 -o last | cut -f2 | sort | uniq -c > ${library}.counts.bed
	"""	

}

process plot_chromatin_state_overlap {

	cache false
	publishDir "${params.results}/figures"

	input:
	file("*") from chromhmm_overlaps_out.toSortedList()

	output:
	file('tn5_chromatin_state_overlap.pdf')

	"""
	plot-chromatin-state-overlap.R --sample_info ${params.sample_info} --out tn5_chromatin_state_overlap.pdf
	"""

}


process read_tf_overlaps {

	memory '10 GB'
	maxForks 20
	validExitStatus 0,141

	input:
	set val(library), file(reads), file(chipseq) from read_tf_overlaps_in.combine(Channel.fromPath(params.chipseq_glob))

	output:
	file("${library}.${tf}.overlaps.txt") into concat_read_tf_overlaps_in

	script:
	tf_match = chipseq.getName() =~ /(.*).bed.gz/
	tf = tf_match[0][1]

	"""
	${IONICE} bedtools intersect -a $reads -b $chipseq -u | wc -l | perl -pe 's/^/${library}\\t${tf}\\t/' > ${library}.${tf}.overlaps.txt
	"""

}

process concat_read_tf_overlaps {
	
	publishDir "${params.results}/read-tf-overlaps"

	input:
	file(x) from concat_read_tf_overlaps_in.toSortedList()

	output:
	file('read-tf-overlaps.txt') into concat_read_tf_overlaps_out

	"""
	cat ${x.join(' ')} | sort | uniq > read-tf-overlaps.txt
	"""

}


process plot_read_tf_overlaps {

	publishDir "${params.results}/figures"

	input:
	file(total_counts) from blacklist_filtered_read_counts_out
	file(overlaps) from concat_read_tf_overlaps_out

	output:
	file('tn5_tf_overlap.pdf')

	"""
	plot-read-tf-overlap.R ${params.sample_info} $overlaps $total_counts
	"""

}


tn5_sensitive_peak_covariates = ['fragment_length_distance', 'median_fragment_length', 'short_mononucleosomal_ratio', 'none']

process tn5_sensitive_peaks {

	publishDir "${params.results}/tn5-sensitive-peaks"
	
	input:
	file(metrics) from tn5_sensitive_peaks_metrics_in
	file(counts) from tn5_sensitive_peaks_counts_in.toSortedList()
	each covariate from tn5_sensitive_peak_covariates

	output:
	file("${prefix_arg}.results.txt") into tn5_sensitive_peaks_out

	script:
	prefix_arg = 'covariates_replicate_subsample_false'
	covariates_arg = 'replicate'
	grep_string = "-w -e hqaa"
	if (covariate != 'none') {
		prefix_arg = "covariates_replicate___${covariate}_subsample_false"
		covariates_arg = 'replicate,' + covariate
		grep_string = '-w -e "hqaa" -e "' + covariate + '"'
	}

	"""
	grep ${grep_string} $metrics > covariates-file.txt
	mkdir counts_dir
	mv ${counts.join(' ')} counts_dir/
	nb_glm.R --prefix $prefix_arg --size_factors hqaa --sample_info ${params.sample_info} --counts counts_dir --tn5_scale log2 --covariates_file covariates-file.txt --covariates $covariates_arg
	"""

}


process analyze_tn5_sensitive_modeling {

	publishDir "${params.results}/figures"

	input:
	file(res) from tn5_sensitive_peaks_out.toSortedList()

	output:
	set file('tn5-sensitivity.txt'), file('tn5_modeling.pdf')
	file('tn5-sensitivity.txt') into tn5_sensitivity

	"""
	analyze-tn5-sensitive-modeling.R
	echo 'hello'
	"""

}


process master_peaks_tf_overlap {

	validExitStatus 0,141

	input:
	set file(master_peaks), file(chipseq) from master_peaks_tf_overlap_in.combine(Channel.fromPath(params.chipseq_glob))

	output:
	file("${tf}.master-peak-overlap.txt") into master_peaks_tf_overlap_out

	script:	
	tf_match = chipseq.getName() =~ /(.*).bed.gz/
	tf = tf_match[0][1]

	"""
	bedtools intersect -c -a $master_peaks -b $chipseq | perl -pe 's/\$/\\t${tf}/' > ${tf}.master-peak-overlap.txt
	"""
	
}

process tn5_sensitivity_vs_tf_overlap_and_peak_size {

	input:
	file(tn5_sensitivity) from tn5_sensitivity
	file(x) from master_peaks_tf_overlap_out.toSortedList()
	file(peak_counts) from tn5_sensitivity_vs_peak_signal_in.toSortedList()

	"""
	cat ${x.join(' ')} > master-peak-overlap.txt
	"""

}


process pca {

	publishDir "${params.results}/pca"

	input:
	file(counts) from pca_in.toSortedList()

	output:
	file('pca.Rda') into pca_out

	"""
	pca.R --counts-dir . --scale --out pca.Rda
	"""

}

process plot_pca {

	publishDir "${params.results}/figures"

	input:
	file(pca_rda) from pca_out

	output:
	file('tn5_pca.pdf')

	"""
	plot-pca.R --sample-info ${params.sample_info} --rda $pca_rda --out tn5_pca.pdf
	"""

}

// compare TSS enrichment using different TSS annotations
// no need for peaks, since we are only interested in the TSS enrichment

ataqv_different_tss_annotations = [params.refseq_tss, params.refseq_housekeeping_tss]

ataqv_in = []
for (library in libraries) {
	ataqv_in << [library, file(get_md_bam(library)), file(get_md_bam(library) + '.bai')]
}

ataqv_different_tss_annotations_bam = Channel.from(ataqv_in)

// compare TSS enrichment using different TSS annotations
process ataqv_different_tss {

	maxForks 10
	validExitStatus 0,141

	input:
	set val(library), file(bam), file(ind) from ataqv_different_tss_annotations_bam
	each tss_annot from ataqv_different_tss_annotations

	output:
	file("${library}-${tss_list}.ataqv.json.gz") into ataqv_different_tss_out

	script:
	tss_list =  (tss_annot == params.refseq_tss) ? 'all_refseq' : 'housekeeping_refseq'

	"""
	ataqv --name ${library}-${tss_list} --metrics-file ${library}-${tss_list}.ataqv.json.gz ${make_excluded_regions_arg(get_genome(library))} --tss-file ${tss_annot} --ignore-read-groups human $bam > ${library}-${tss_list}.ataqv.out
	"""

}

process plot_ataqv_different_tss {

	publishDir "${params.results}/figures"

	input:
	file(ataqv_json) from ataqv_different_tss_out.toSortedList()

	output:
	file("tss_enrichment_housekeeping_vs_all.pdf")

	"""
	mkarv tmp-session ${ataqv_json.join(' ')}
	extractAtaqvMetric.py --files tmp-session/data/*.json.gz --metrics tss_enrichment > tss_enrichment.txt
	plot-tss-enrichment-different-annotations.R ${params.sample_info} tss_enrichment.txt
	"""

}
