rule make_windows:
    input:
        chrom_sizes=config["refs"]["chrom_sizes"],
    output:
        "work/windows/genome.{w}bp.bed",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools makewindows -g {input.chrom_sizes} -w {wildcards.w} > {output}"
