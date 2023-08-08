rule multiqc:
    input:
        rules.run_metrics.input
    output:
        multiqcdir + "/QC.html"
    log:
        multiqcdir + "/logs/multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    params:
        input_path = [fastpdir, markdupdir, metricsdir],
        output_path = multiqcdir,
        output_name = "QC.html",
        config = "config/multiqc.yaml"
    shell:
        "(set -x; multiqc -f "
        "  -o {params.output_path} "
        "  -n {params.output_name} "
        "  -c {params.config} "
        "  {params.input_path} "
        ") &> {log}"

rule run_multiqc:
    input: rules.multiqc.output
