rule make_windows:
    input:
        chrom_sizes="work/refs/hg19_nochr.chrom.sizes",
    output:
        "work/windows/genome.{w}bp.bed",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools makewindows -g {input.chrom_sizes} -w {wildcards.w} > {output}"
