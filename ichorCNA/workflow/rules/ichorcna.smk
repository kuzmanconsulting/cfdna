def _ichorcna_script():
    """Resolve the runIchorCNA.R path from inside the conda env at runtime."""
    return "$(Rscript -e 'cat(system.file(\"extdata\", \"runIchorCNA.R\", package=\"ichorCNA\"))')"


def _gc_wig():
    bin_kb = config["params"]["bin_size"] // 1000
    genome = config["refs"]["genome"]
    return f"$(Rscript -e 'cat(system.file(\"extdata\", \"gc_{genome}_{bin_kb}kb.wig\", package=\"ichorCNA\"))')"


def _map_wig():
    bin_kb = config["params"]["bin_size"] // 1000
    genome = config["refs"]["genome"]
    # ichorCNA ships mappability WIGs with a bowtie suffix; try the common name
    return (
        f"$(Rscript -e 'f <- system.file(\"extdata\","
        f" \"map_{genome}_{bin_kb}kb_150mer_{genome}_bowtie.wig\","
        f" package=\"ichorCNA\"); if (nchar(f)==0) stop(\"map WIG not found\") else cat(f)')"
    )


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
            --outDir {params.outdir} \\
            > {log} 2>&1
        """
