import os
import re

# assess the genome annotations ...
ReferenceFasta = config["genome_fasta"]
GenomeGFF  = config["genome_annot"]
SAMPLES=config["samples"]

rule all:
  input:
    expand("Results/Salmon/{sample}/quant.sf", sample=SAMPLES),
    "Results/Pinfish/clustered_transcripts_collapsed.gff",
    "Results/GffCompare/nanopore.combined.gtf"

rule dump_versions:
    output:
        ver = "versions.txt"
    conda:
        "environment.yml"
    shell:
        "conda list > {output.ver}"

# build known splice junction bed file
rule SpliceJunctionIndex:
    input:
        "ReferenceData/"+GenomeGFF
    output:
        junc_bed = "ReferenceData/junctions.bed"
    shell:
        "paftools.js gff2bed -j {output.junc_bed} {input}"

# build minimap2 index
rule Minimap2Index:
    input:
        genome = "ReferenceData/"+ReferenceFasta
    output:
        index = "Results/Minimap2/"+ReferenceFasta+".mmi"
    params:
        opts = config["minimap_index_opts"]
    threads: config["threads"]
    shell:
      "minimap2 -t {threads} {params.opts} -I 1000G -d {output.index} {input.genome}"

rule FilterReads:
    input:
        "RawData/{sample}.fastq"
    output:
        "FilteredData/{sample}.fastq",
    params:
        config["min_mean_q"]
    shell:
        "filtlong --min_mean_q {params} {input} > {output}"

rule Pychopper:
  input:
    "FilteredData/{sample}.fastq"
  output:
    pdf = "Results/Pychopper/{sample}.pychopper_report.pdf",
    fastq = "Results/Pychopper/{sample}.pychop.fastq",
    stats = "Results/Pychopper/{sample}.pychop.stats",
    scores = "Results/Pychopper/{sample}.pychop.scores",
    unclass = "Results/Pychopper/{sample}.unclassified.fastq",
  params:
    opts = config["porechop_heu_stringency"]
  run:
    shell("cdna_classifier.py -b ReferenceData/cdna_barcodes.fas -x -r {output.pdf} -S {output.stats} -A {output.scores} -u {output.unclass} -l {params.opts} {input} {output.fastq}")

rule Minimap2: ## map reads using minimap2
    input:
       index = rules.Minimap2Index.output.index,
       fastq = expand("Results/Pychopper/{sample}.pychop.fastq", sample=SAMPLES) if (config["pychopper"]==True) else expand("RawData/{sample}.fastq", sample=SAMPLES)
    output:
       bam = "Results/Minimap2/merged.mapping.bam"
    params:
        opts = config["minimap2_opts"],
        min_mq = config["minimum_mapping_quality"],
    threads: config["threads"]
    shell:"""
    minimap2 -t {threads} -ax splice {params.opts} {input.index} {input.fastq}\
    | samtools view -q {params.min_mq} -F 2304 -Sb | samtools sort -@ {threads} - -o {output.bam};
    samtools index {output.bam}
    """

rule PinfishRawBAM2GFF: ## convert BAM to GFF
    input:
        bam = rules.Minimap2.output.bam
    output:
        raw_gff = "Results/Pinfish/raw_transcripts.gff"
    params:
        opts = config["spliced_bam2gff_opts"]
    threads: config["threads"]
    shell:
        "spliced_bam2gff {params.opts} -t {threads} -M {input.bam} > {output.raw_gff}"


rule PinfishClusterGFF: ## cluster transcripts in GFF
    input:
        raw_gff = rules.PinfishRawBAM2GFF.output.raw_gff
    output:
        cls_gff = "Results/Pinfish/clustered_pol_transcripts.gff",
        cls_tab = "Results/Pinfish/cluster_memberships.tsv",
    params:
        c = config["minimum_cluster_size"],
        d = config["exon_boundary_tolerance"],
        e = config["terminal_exon_boundary_tolerance"],
        min_iso_frac = config["minimum_isoform_percent"],
    threads: config["threads"]
    shell:
        "cluster_gff -p {params.min_iso_frac} -t {threads} -c {params.c} -d {params.d} -e {params.e} -a {output.cls_tab} {input.raw_gff} > {output.cls_gff}"


rule PinfishCollapseRawPartials: ## collapse clustered read artifacts
    input:
        cls_gff = rules.PinfishClusterGFF.output.cls_gff
    output:
        cls_gff_col = "Results/Pinfish/clustered_transcripts_collapsed.gff"
    params:
        d = config["collapse_internal_tol"],
        e = config["collapse_three_tol"],
        f = config["collapse_five_tol"],
    shell:
       "collapse_partials -d {params.d} -e {params.e} -f {params.f} {input.cls_gff} > {output.cls_gff_col}"


rule PinfishPolishClusters: ## polish read clusters
    input:
        cls_gff = rules.PinfishClusterGFF.output.cls_gff,
        cls_tab = rules.PinfishClusterGFF.output.cls_tab,
        bam = rules.Minimap2.output.bam
    output:
        pol_trs = "Results/Pinfish/polished_transcripts.fas"
    params:
        c = config["minimum_cluster_size"]
    threads: config["threads"]
    shell:
        "polish_clusters -t {threads} -a {input.cls_tab} -c {params.c} -o {output.pol_trs} {input.bam}"


rule MinimapPolishedClusters: ## map polished transcripts to genome
    input:
       index = rules.Minimap2Index.output.index,
       fasta = rules.PinfishPolishClusters.output.pol_trs,
    output:
       pol_bam = "Results/Minimap2/polished_reads_aln_sorted.bam"
    params:
        extra = config["minimap2_opts_polished"]
    threads: config["threads"]
    shell:"""
    minimap2 -t {threads} {params.extra} -ax splice {input.index} {input.fasta}\
    | samtools view -Sb -F 2304 | samtools sort -@ {threads} - -o {output.pol_bam};
    samtools index {output.pol_bam}
    """


rule PinfishPolishedBAM2GFF: ## convert BAM of polished transcripts to GFF
    input:
        bam = rules.MinimapPolishedClusters.output.pol_bam
    output:
        pol_gff = "Results/Pinfish/polished_transcripts.gff"
    params:
        extra = config["spliced_bam2gff_opts_pol"]
    threads: config["threads"]
    shell:
        "spliced_bam2gff {params.extra} -t {threads} -M {input.bam} > {output.pol_gff}"


rule PinfishCollapsePolishedPartials: ## collapse polished read artifacts
    input:
        pol_gff = rules.PinfishPolishedBAM2GFF.output.pol_gff
    output:
        pol_gff_col = "Results/Pinfish/polished_transcripts_collapsed.gff"
    params:
        d = config["collapse_internal_tol"],
        e = config["collapse_three_tol"],
        f = config["collapse_five_tol"],
    shell:
        "collapse_partials -d {params.d} -e {params.e} -f {params.f} {input.pol_gff} > {output.pol_gff_col}"


rule PrepareCorrectedTranscriptomeFasta: ## Generate corrected transcriptome.
    input:
        genome = "ReferenceData/"+ReferenceFasta,
        gff = rules.PinfishCollapsePolishedPartials.output.pol_gff_col,
    output:
        fasta = "Results/Pinfish/corrected_transcriptome_polished_collapsed.fas"
    shell:"""
    gffread -g {input.genome} -w {output.fasta} {input.gff}
    """


rule GffCompare:
    input:
        reference = "ReferenceData/"+GenomeGFF,
        exptgff = rules.PinfishCollapsePolishedPartials.output.pol_gff_col
    output:
        "Results/GffCompare/nanopore.combined.gtf",
        "Results/GffCompare/nanopore.loci",
        "Results/GffCompare/nanopore.redundant.gtf",
        "Results/GffCompare/nanopore.stats",
        "Results/GffCompare/nanopore.tracking"
    shell:
        "gffcompare -r {input.reference} -R -M -C -K -o Results/GffCompare/nanopore {input.exptgff}"

#map to transcriptome
rule Map2Transcriptome:
    input:
       target = rules.PrepareCorrectedTranscriptomeFasta.output.fasta,
       fastq = "FilteredData/{sample}.fastq"
    output:
       sam = "Results/Quantification/{sample}.sam"
    threads: config["threads"]
    priority: 10
    shell:
      "minimap2 -a -t {threads} {input.target} {input.fastq} > {output.sam}"

rule SalmonQuantifyMapped:
    input:
        target = rules.PrepareCorrectedTranscriptomeFasta.output.fasta,
        sam = "Results/Quantification/{sample}.sam"
    output:
        quant = "Results/Salmon/{sample}/quant.sf",
    log:
        'Results/Salmon/{sample}.log'
    shell:
        "salmon quant -l U -a {input.sam} -t {input.target} -o Results/Quantification/{wildcards.sample} --noErrorModel --writeUnmappedNames > {log} "
