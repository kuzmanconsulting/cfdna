rule frag_length_bins:
    input:
        bam="work/{sample}/{sample}.md.bam",
        bai="work/{sample}/{sample}.md.bam.bai",
    output:
        tsv="results/{sample}/frag_lengths/{sample}.frag_length_bins.tsv",
    params:
        mapq=config["params"]["mapq"],
    threads: config["params"]["threads"]["finaletoolkit"]
    conda:
        "../envs/finaletoolkit.yaml"
    log:
        "work/logs/frag_length_bins_{sample}.log",
    shell:
        "finaletoolkit frag-length-bins "
        "-q {params.mapq} --bin-size 1 "
        "-o {output.tsv} "
        "{input.bam} > {log} 2>&1"


rule plot_frag_lengths:
    input:
        tsv="results/{sample}/frag_lengths/{sample}.frag_length_bins.tsv",
    output:
        png="results/{sample}/frag_lengths/{sample}.frag_length_dist.png",
    conda:
        "../envs/finaletoolkit.yaml"
    log:
        "work/logs/plot_frag_lengths_{sample}.log",
    script:
        "../scripts/plot_frag_lengths.py"
