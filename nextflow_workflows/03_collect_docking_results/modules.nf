process CALCULATE_RMSD{
    publishDir "${params.dockedLigandRMSDs}/${name}", mode: 'copy', overwrite: true
    conda "${params.asap}"
    tag "calculate-rmsds ${uuid}"
    label 'cpushort'

    input:
    tuple val(name), val(uuid), path(docked_dir), path(ligand_file_3d)

    output:
    path("*.csv"), emit: rmsd_csv

    script:
    """
    python3 "${params.scripts}"/calculate_rmsd_from_docking_results.py \
    -d "${docked_dir}" \
    -l "${ligand_file_3d}" \
    -o "${name}_${uuid}_rmsd_results.csv"
    """
}
process COMBINE_AND_PROCESS_RESULTS {
    publishDir "${params.combinedDockingResultsPath}", mode: 'copy', overwrite: true
    conda "${params.drugforge}"
    tag "combine-and-process-results ${name}"
    label 'cpushort'
    memory 128.GB

    input:
    path(dockedLigandRMSDs)
    val(name)

    output:
    path("*")

    script:
    """
    python3 "${params.scripts}"/combine_and_process_results.py \
    ${dockedLigandRMSDs.join(' ')} \
    --tc-data "${params.chemicalSimilarityData}/tanimoto_combo/tanimoto_combo.csv" \
    --ecfp-data "${params.chemicalSimilarityData}/ecfp_tanimoto/fingerprint_similarities.csv" \
    --mcs-data "${params.chemicalSimilarityData}/mcs_tanimoto/mcs_tanimoto.csv" \
    --date-dict "${params.dateDictPath}" \
    --structure-cmpd-dict "${params.dataPath}/cmpd_date_dict/structure_to_cmpd_dict.json" \
    --scaffold-data "${params.genericScaffoldPath}" \
    --output-file-prefix "${name}" \
    --protein-rmsd-data "${params.combinedProteinRmsdPath}" \
    --deduplicate
    """
}

process CONVERT_TO_DOCKING_DATA_MODEL {
    publishDir "${params.combinedDockingResultsPath}", mode: 'copy', overwrite: true
    conda "${params.harbor}"
    tag "combine-and-process-results ${method}"
    errorStrategy = { task.exitStatus in [137,140,143,247] ? 'retry' : 'finish' }
    maxRetries 3
    // Dynamic memory allocation
    memory { task.attempt > 1 ? (2 ** (task.attempt - 1)) * 8.GB : 8.GB }
    // Dynamic time allocation
    time { task.attempt > 1 ? (2 ** (task.attempt - 1)) * 2.h : 2.h }

    input:
    path(input_csv)
    val(method)

    output:
    path("*.parquet"), emit: docking_data_model_dataframe
    path("*.json"), emit: docking_data_model_schema

    script:
    """
    python3 "${params.scripts}"/convert_to_docking_data_model.py \
    --input-csv ${input_csv} \
    --output-file-prefix "${method}"_combined_results \
    """
}