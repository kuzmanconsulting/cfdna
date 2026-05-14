rule samtools_sort:
    input:
        input_bam,
    output:
        bam=temp("work/{sample}/{sample}.sorted.bam"),
    threads: config["params"]["threads"]["preprocess"]
    conda:
        "../envs/samtools.yaml"
    log:
        "work/logs/sort_{sample}.log",
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input} > {log} 2>&1"


rule mark_duplicates:
    input:
        "work/{sample}/{sample}.sorted.bam",
    output:
        bam="work/{sample}/{sample}.md.bam",
        metrics="results/{sample}/metrics/{sample}.dup_metrics.txt",
    resources:
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
    conda:
        "../envs/picard.yaml"
    log:
        "work/logs/markdup_{sample}.log",
    shell:
        """
        picard -Xms{resources.java_min_gb}g -Xmx{resources.java_max_gb}g \
            MarkDuplicates \
            I={input} \
            O={output.bam} \
            M={output.metrics} \
            CREATE_INDEX=false \
            VALIDATION_STRINGENCY=LENIENT \
            TMP_DIR=work/tmp \
            > {log} 2>&1
        """


rule samtools_index_md:
    input:
        "work/{sample}/{sample}.md.bam",
    output:
        "work/{sample}/{sample}.md.bam.bai",
    threads: config["params"]["threads"]["preprocess"]
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -@ {threads} {input}"
