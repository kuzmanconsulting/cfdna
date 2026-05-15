"""One-time reference conversion: restrict refs to the BAM's @SQ contigs.

The Snyder BAMs are b37/1000G-style (bare contig names, mito = MT, plus
hs37d5/NC_007605/GL000*.1 decoys). UCSC hg19 refs use chr-prefixed names,
mito = chrM, and carry alt-haplotype contigs (e.g. chr6_ssto_hap7). Building
FinaleToolkit intervals from ref contigs absent in the BAM crashes with
`ValueError: invalid contig`. Intersecting against the BAM @SQ list keeps
exactly the contigs present in both (1..22, X, Y) and drops mito + decoys +
alt contigs, which is the desired set for cfDNA fragmentomics anyway.
"""


rule bam_contigs:
    """Authoritative contig list = BAM @SQ headers."""
    input:
        bam=os.path.join(bam_dir, f"{srr_for(SAMPLES[0])}.bam"),
    output:
        "work/refs/bam_contigs.txt",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -H {input.bam} | awk '$1==\"@SQ\"' "
        "| sed 's/.*SN://; s/\\t.*//' > {output}"


rule strip_chr_chrom_sizes:
    """chrom.sizes restricted to contigs present in the BAM header."""
    input:
        sizes=config["refs"]["chrom_sizes"],
        contigs="work/refs/bam_contigs.txt",
    output:
        "work/refs/hg19_nochr.chrom.sizes",
    shell:
        "sed 's/^chr//' {input.sizes} | grep -Fwf {input.contigs} - > {output}"


rule strip_chr_2bit:
    """2bit subset to the same contigs (2bit ⊆ chrom.sizes ⊆ BAM @SQ)."""
    input:
        twobit=config["refs"]["genome_2bit"],
        sizes="work/refs/hg19_nochr.chrom.sizes",
    output:
        "work/refs/hg19_nochr.2bit",
    conda:
        "../envs/ucsc.yaml"
    shell:
        "cut -f1 {input.sizes} | sed 's/^/chr/' "
        "| twoBitToFa -seqList=/dev/stdin {input.twobit} /dev/stdout "
        "| sed 's/>chr/>/' | faToTwoBit /dev/stdin {output}"
