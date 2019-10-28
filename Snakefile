
ReferenceFasta = config["genome_fasta"]
GenomeGFF  = config["genome_annot"]
SAMPLES=config["samples"]

rule all:
  input:
    "Results/Quantification/all_counts.txt",
    expand("IGV/{sample}.genome.bam", sample=SAMPLES),
    "Results/Pinfish/clustered_transcripts_collapsed.gff",

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
        "paftools.js gff2bed {input} > {output.junc_bed}"

# build minimap2 index
rule Minimap2Index:
    input:
        genome = "ReferenceData/"+ReferenceFasta
    output:
        index = "Results/Minimap2/"+ReferenceFasta+".mmi"
    params:
        opts = config["minimap2_index_opts"]
    threads: config["threads"]
    shell:
      "minimap2 -t {threads} {params.opts} -I 1000G -d {output.index} {input.genome}"

rule FilterReads:
    input:
        lambda wildcards: "RawData/"+SAMPLES[wildcards.sample]+".fastq"
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
    shell("cdna_classifier.py -x True -b ReferenceData/cdna_barcodes.fas -r {output.pdf} -S {output.stats} -A {output.scores} -u {output.unclass} -l {params.opts} {input} {output.fastq}")

rule Minimap2Pinfish: ## map reads using minimap2
    input:
       index = rules.Minimap2Index.output.index,
       fastq = expand("Results/Pychopper/{sample}.pychop.fastq", sample=SAMPLES) if (config["pychopper"]==True) else expand("FilteredData/{sample}.fastq", sample=SAMPLES),
       use_junc = rules.SpliceJunctionIndex.output.junc_bed if config["minimap2_opts_junction"] else ""
    output:
       bam = "Results/Minimap2/merged.mapping.bam"
    params:
        opts = config["minimap2_opts"],
        min_mq = config["minimum_mapping_quality"],
    threads: config["threads"]
    shell:"""
    minimap2 -t {threads} -ax splice {params.opts} --junc-bed {input.use_junc} {input.index} {input.fastq}\
    | samtools view -q {params.min_mq} -F 2304 -Sb | samtools sort -@ {threads} - -o {output.bam};
    samtools index {output.bam}
    """

rule PinfishRawBAM2GFF: ## convert BAM to GFF
    input:
        bam = rules.Minimap2Pinfish.output.bam
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
        bam = rules.Minimap2Pinfish.output.bam
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
        fasta = "Results/Pinfish/corrected_transcriptome_polished_collapsed.fas",
    shell:
        "gffread -g {input.genome} -w {output.fasta} {input.gff}"

# build minimap2 index for transcriptome
rule Transcriptome2Index:
    input:
        genome = rules.PrepareCorrectedTranscriptomeFasta.output.fasta
    output:
        index = "Results/Minimap2/Transcriptome.mmi"
    threads: config["threads"]
    shell:
      "minimap2 -t {threads} -I 1000G -d {output.index} {input.genome}"

rule GffCompare:
    input:
        reference = "ReferenceData/"+GenomeGFF,
        exptgff = rules.PinfishCollapsePolishedPartials.output.pol_gff_col
    output:
        "Results/GffCompare/nanopore.combined.gtf",
        "Results/GffCompare/nanopore.loci",
        "Results/GffCompare/nanopore.redundant.gtf",
        "Results/GffCompare/nanopore.stats",
        "Results/GffCompare/nanopore.tracking",
    shell:
        "gffcompare -r {input.reference} -R -A -K -o Results/GffCompare/nanopore {input.exptgff}"

rule Minimap2Genome: ## map reads using minimap2
    input:
       index = rules.Minimap2Index.output.index,
       fastq = "FilteredData/{sample}.fastq",
       use_junc = rules.SpliceJunctionIndex.output.junc_bed if config["minimap2_opts_junction"] else ""
    output:
       bam = "IGV/{sample}.genome.bam"
    params:
        opts = config["minimap2_opts"],
    threads: config["threads"]
    shell:"""
    minimap2 -t {threads} -ax splice {params.opts} --junc-bed {input.use_junc} {input.index} {input.fastq}\
    | samtools view -F 260 -Sb | samtools sort -@ {threads} - -o {output.bam};
    samtools index {output.bam}
    """

#map to transcriptome
rule Map2Transcriptome:
    input:
       index = rules.Transcriptome2Index.output.index,
       fastq = "FilteredData/{sample}.fastq"
    output:
       bam = "Results/Quantification/{sample}.bam",
       sbam = "Results/Quantification/{sample}.sorted.bam"
    threads: config["threads"]
    params:
        msec = config["maximum_secondary"],
        psec = config["secondary_score_ratio"]
    priority: 10
    shell:"""
        minimap2 -ax map-ont -t {threads} -p {params.psec} -N {params.msec} {input.index} {input.fastq} | samtools view -Sb > {output.bam};
        samtools sort -@ {threads} {output.bam} -o {output.sbam};
        samtools index {output.sbam};
        """

rule ExtractPrimaryMapping:
    input:
        "Results/Quantification/{sample}.sorted.bam"
    output:
        "IGV/{sample}.transcriptome.bam"
    shell:"""
        samtools view -b -F 260 {input} > {output};
        samtools index {output}
        """

rule CountTranscripts:
    input:
        bam = "IGV/{sample}.transcriptome.bam"
    output:
        quant = "Results/Quantification/{sample}.counts"
    shell:"""
        echo -e "counts\ttranscript" > {output.quant};
        samtools view {input.bam} | cut -f3 | uniq -c | grep -v "*" | sed -e 's/^[ \t]*//' | sed 's/ /\t/' >> {output.quant}
        """

rule CalculateTPM:
    input:
        annotation = "ReferenceData/"+GenomeGFF,
        transcriptome = "Results/GffCompare/nanopore.combined.gtf",
        counts = expand("Results/Quantification/{sample}.counts", sample=SAMPLES),
    output:
        counts = "Results/Quantification/all_counts.txt"
    script:
        "scripts/concatenateCounts.py"

rule SalmonQuantifyMapped:
    input:
        target = rules.PrepareCorrectedTranscriptomeFasta.output.fasta,
        bam = "Results/Quantification/{sample}.bam"
    output:
        quant = "Results/Salmon/{sample}/quant.sf"
    params:
        dir = "Results/Salmon/{sample}"
    threads: config["threads"]
    shell:
        "salmon quant --noErrorModel -p {threads} -l U -a {input.bam} -t {input.target} -o {params.dir}"

rule SpliceJunction:
    input:
        "IGV/{sample}.genome.bam"
    output:
        "Results/Splicing/{sample}_junction.bed"
    shell:
        "regtools junctions extract -o {output} -s 0 {input}"

rule Strandness:
  input:
    "FilteredData/{sample}.fastq"
  output:
    pdf = "Results/Strand/{sample}.pychopper_report.pdf",
    fastq = "Results/Strand/{sample}.pychop.fastq",
    stats = "Results/Strand/{sample}.pychop.stats",
    scores = "Results/Strand/{sample}.pychop.scores",
    unclass = "Results/Strand/{sample}.unclassified.fastq",
  run:
    shell("cdna_classifier.py -m edlib -b ReferenceData/cdna_barcodes.fas -r {output.pdf} -S {output.stats} -A {output.scores} -u {output.unclass} {input} {output.fastq}")

rule CalculatePsi:
    input:
        expand("Results/Splicing/{sample}_junction.bed", sample=SAMPLES),
        expand("Results/Strand/{sample}.pychopper_report.pdf", sample=SAMPLES)
