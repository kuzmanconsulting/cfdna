"""One-time reference conversion: strip 'chr' prefix to match Snyder BAM contig names."""


rule strip_chr_chrom_sizes:
    """Remove chr prefix from chrom.sizes so window BED contigs match BAM headers."""
    input:
        config["refs"]["chrom_sizes"],
    output:
        "work/refs/hg19_nochr.chrom.sizes",
    shell:
        "sed 's/^chr//' {input} > {output}"


rule strip_chr_2bit:
    """Convert chr-prefixed hg19.2bit to bare-contig-name version for FinaleToolkit."""
    input:
        config["refs"]["genome_2bit"],
    output:
        "work/refs/hg19_nochr.2bit",
    conda:
        "../envs/ucsc.yaml"
    shell:
        "twoBitToFa {input} /dev/stdout | sed 's/>chr/>/' | faToTwoBit /dev/stdin {output}"
