#!/usr/bin/env nextflow

params.input_structures = 'data/BILELIB19_spectralID_and_SMILES.tsv'

params.publishdir = "nf_output"
TOOL_FOLDER = "$baseDir/bin"

process calculate {
    //errorStrategy 'ignore'
    //time '4h'
    //memory { 12.GB }

    conda "$TOOL_FOLDER/conda_env.yml"

    publishDir "$params.publishdir", mode: 'copy'
    
    input:
    file(smiles_file) from Channel.fromPath(params.input_structures)
    each x from Channel.from( 0..2000 )

    output:
    file "*output.json" optional true into _query_results_ch

    """
    echo $x
    python $TOOL_FOLDER/run_comparisons.py \
        "$smiles_file" \
        "${x}_output.json" \
        --node_current ${x} \
        --node_total 2001
    """
}