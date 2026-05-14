def windows_bed(wildcards):
    return f"work/windows/genome.{config['params']['window_size_bp']}bp.bed"


rule interval_end_motifs:
    input:
        bam="work/{sample}/filter/{sample}.{cls}.bam",
        bai="work/{sample}/filter/{sample}.{cls}.bam.bai",
        ref2bit="work/refs/hg19_nochr.2bit",
        intervals=windows_bed,
    output:
        "results/{sample}/interval_end_motifs/{sample}.{cls}.interval_end_motifs.tsv",
    params:
        kmer=config["params"]["kmer"],
        mapq=config["params"]["mapq"],
    threads: config["params"]["threads"]["finaletoolkit"]
    conda:
        "../envs/finaletoolkit.yaml"
    log:
        "work/logs/interval_end_motifs_{sample}_{cls}.log",
    shell:
        "finaletoolkit interval-end-motifs "
        "-k {params.kmer} -q {params.mapq} "
        "-w {threads} "
        "-o {output} "
        "{input.bam} {input.ref2bit} {input.intervals} > {log} 2>&1"


rule interval_mds:
    input:
        "results/{sample}/interval_end_motifs/{sample}.{cls}.interval_end_motifs.tsv",
    output:
        "results/{sample}/interval_mds/{sample}.{cls}.interval_mds.bed",
    conda:
        "../envs/finaletoolkit.yaml"
    log:
        "work/logs/interval_mds_{sample}_{cls}.log",
    shell:
        "finaletoolkit interval-mds {input} -o {output} > {log} 2>&1"
