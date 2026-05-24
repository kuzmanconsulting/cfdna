def _ichorcna_script():
    """Path to the bundled runIchorCNA.R script."""
    return os.path.join(workflow.basedir, "workflow", "scripts", "runIchorCNA.R")


def _gc_wig():
    bin_kb = config["params"]["bin_size"] // 1000
    genome = config["refs"]["genome"]
    return f"$(Rscript -e 'cat(system.file(\"extdata\", \"gc_{genome}_{bin_kb}kb.wig\", package=\"ichorCNA\"))')"


def _map_wig():
    bin_kb = config["params"]["bin_size"] // 1000
    genome = config["refs"]["genome"]
    return f"$(Rscript -e 'cat(system.file(\"extdata\", \"map_{genome}_{bin_kb}kb.wig\", package=\"ichorCNA\"))')"


rule gather_tumor_fractions:
    input:
        expand("results/{sample}/{sample}.params.txt", sample=SAMPLES),
    output:
        csv="results/tumor_fractions.csv",
    run:
        import csv, re

        rows = []
        for params_file in input:
            sample = None
            tf = ploidy = gc_mad = None
            with open(params_file) as f:
                for line in f:
                    line = line.rstrip()
                    if re.match(r"^[A-Za-z0-9_-]+$", line) and sample is None:
                        # second bare identifier block is the sample name
                        pass
                    m = re.match(r"^(\S.*?)\s*:\s*(.+)$", line)
                    if m:
                        key, val = m.group(1).strip(), m.group(2).strip()
                        if key == "Tumor Fraction":
                            tf = val
                        elif key == "Ploidy":
                            ploidy = val
                        elif key == "GC-Map correction MAD":
                            gc_mad = val
            # derive sample name from file path
            sample = os.path.basename(os.path.dirname(params_file))
            rows.append({"sample": sample, "tumor_frac": tf,
                         "ploidy": ploidy, "gc_mad": gc_mad})

        rows.sort(key=lambda r: r["sample"])
        with open(output.csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["sample", "tumor_frac", "ploidy", "gc_mad"])
            writer.writeheader()
            writer.writerows(rows)


rule plot_tumor_fractions:
    input:
        csv="results/tumor_fractions.csv",
    output:
        pdf="results/tumor_fractions.pdf",
    conda:
        "../envs/plots.yaml"
    log:
        "work/logs/plot_tumor_fractions.log",
    shell:
        "python3 {workflow.basedir}/workflow/scripts/plot_tumor_fractions.py"
        " {input.csv} {output.pdf} > {log} 2>&1"


rule ichorcna:
    input:
        wig="results/{sample}/{sample}.wig",
    output:
        params_txt="results/{sample}/{sample}.params.txt",
        seg="results/{sample}/{sample}.seg",
    params:
        outdir=lambda wc: f"results/{wc.sample}",
        script=_ichorcna_script(),
        gc_wig=_gc_wig(),
        map_wig=_map_wig(),
        genome=config["refs"]["genome"],
        chrs=config["params"]["chrs"],
        normal=config["params"]["normal_states"],
        ploidy=config["params"]["ploidy"],
        max_cn=config["params"]["max_cn"],
        include_homd=config["params"]["include_homd"],
        txne=config["params"]["txne"],
        txn_strength=config["params"]["txn_strength"],
        plot_type=config["params"]["plot_file_type"],
    resources:
        mem_mb=config["resources"]["ichorcna_mem_mb"],
    log:
        "work/logs/ichorcna_{sample}.log",
    shell:
        """
        {{
            GC_WIG={params.gc_wig}
            MAP_WIG={params.map_wig}
            Rscript {params.script} \\
                --id {wildcards.sample} \\
                --WIG {input.wig} \\
                --gcWig "$GC_WIG" \\
                --mapWig "$MAP_WIG" \\
                --genomeBuild {params.genome} \\
                --chrs '{params.chrs}' \\
                --chrTrain '{params.chrs}' \\
                --normal '{params.normal}' \\
                --ploidy '{params.ploidy}' \\
                --maxCN {params.max_cn} \\
                --includeHOMD {params.include_homd} \\
                --txnE {params.txne} \\
                --txnStrength {params.txn_strength} \\
                --plotFileType {params.plot_type} \\
                --outDir {params.outdir}
        }} > {log} 2>&1
        """
