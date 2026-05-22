def _chrs_for_read_counter(cfg):
    """Convert R-style 'c(1:22)' to comma-separated list for readCounter CLI."""
    raw = cfg["params"]["chrs"]
    # Handle c(1:22) pattern
    import re
    m = re.match(r"c\((\d+):(\d+)\)", raw.strip())
    if m:
        return ",".join(str(i) for i in range(int(m.group(1)), int(m.group(2)) + 1))
    # Fallback: strip c() and return as-is
    return raw.strip().lstrip("c(").rstrip(")")


CHRS_RC = _chrs_for_read_counter(config)


rule read_counter:
    input:
        bam=input_bam,
    output:
        wig="results/{sample}/{sample}.wig",
    params:
        bin=config["params"]["bin_size"],
        mapq=config["params"]["mapq"],
        chrs=CHRS_RC,
    conda:
        "../envs/hmmcopy.yaml"
    log:
        "work/logs/readCounter_{sample}.log",
    shell:
        "readCounter"
        " --window {params.bin}"
        " --quality {params.mapq}"
        " --chromosome {params.chrs}"
        " {input.bam}"
        " > {output.wig}"
        " 2> {log}"
