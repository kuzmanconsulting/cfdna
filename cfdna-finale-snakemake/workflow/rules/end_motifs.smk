rule filter_file:
    """Filter markdup BAM by MAPQ + fragment length class.

    filter-file auto-enforces: properly-paired, read1, not secondary/supplementary,
    not QC-flagged (so MarkDuplicates-flagged dups are dropped here).
    """
    input:
        bam="work/{sample}/{sample}.md.bam",
        bai="work/{sample}/{sample}.md.bam.bai",
    output:
        bam="work/{sample}/filter/{sample}.{cls}.bam",
    params:
        mapq=config["params"]["mapq"],
        lflags=lambda wc: length_flags(wc.cls),
    threads: config["params"]["threads"]["finaletoolkit"]
    conda:
        "../envs/finaletoolkit.yaml"
    log:
        "work/logs/filter_file_{sample}_{cls}.log",
    shell:
        "finaletoolkit filter-file "
        "-q {params.mapq} {params.lflags} "
        "-w {threads} "
        "-o {output.bam} {input.bam} > {log} 2>&1"


rule samtools_index_filtered:
    input:
        "work/{sample}/filter/{sample}.{cls}.bam",
    output:
        "work/{sample}/filter/{sample}.{cls}.bam.bai",
    threads: config["params"]["threads"]["finaletoolkit"]
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -@ {threads} {input}"


rule end_motifs:
    input:
        bam="work/{sample}/filter/{sample}.{cls}.bam",
        bai="work/{sample}/filter/{sample}.{cls}.bam.bai",
        ref2bit="work/refs/hg19_nochr.2bit",
    output:
        "results/{sample}/end_motifs/{sample}.{cls}.tsv",
    params:
        kmer=config["params"]["kmer"],
        mapq=config["params"]["mapq"],
    threads: config["params"]["threads"]["finaletoolkit"]
    conda:
        "../envs/finaletoolkit.yaml"
    log:
        "work/logs/end_motifs_{sample}_{cls}.log",
    shell:
        "finaletoolkit end-motifs "
        "-k {params.kmer} -q {params.mapq} "
        "-w {threads} "
        "-o {output} "
        "{input.bam} {input.ref2bit} > {log} 2>&1"


rule mds:
    input:
        "results/{sample}/end_motifs/{sample}.{cls}.tsv",
    output:
        "results/{sample}/end_motifs/{sample}.{cls}.mds.txt",
    conda:
        "../envs/finaletoolkit.yaml"
    log:
        "work/logs/mds_{sample}_{cls}.log",
    shell:
        "finaletoolkit mds {input} > {output} 2> {log}"
