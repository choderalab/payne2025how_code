process GENERATE_DATE_DICTIONARY {
    publishDir "${params.dataPath}", mode: 'copy', overwrite: true
    conda "${params.asap}"
    tag "generate-date-dictionary"

    output:
    path "cmpd_date_dict"
    path "cmpd_date_dict/date_dict.json", emit: structure_to_date_dict
    path "cmpd_date_dict/structure_to_cmpd_dict.json", emit: structure_to_cmpd_dict

    script:
    """
    python3 "${params.scripts}"/generate_date_dict.py --fragalysis-dir "${params.curatedFragalysis}" --output-dir cmpd_date_dict
    """
}
process CALCULATE_ECFP_TANIMOTO {
    publishDir "${params.chemicalSimilarityData}", mode: 'copy', overwrite: true
    conda "${params.asap}"
    tag "calculate-ecfp-tanimoto"

    input:
    path(ligand_file_3d)

    output:
    path("ecfp_tanimoto"), emit: ecfp_tanimoto

    script:
    """
    python3 "${params.scripts}"/calculate_ecfp_tanimoto.py --ref-ligand-sdf "${ligand_file_3d}" --output-dir ecfp_tanimoto
    """
}
process CALCULATE_MCS_TANIMOTO {
    publishDir "${params.chemicalSimilarityData}", mode: 'copy', overwrite: true
    conda "${params.asap}"
    tag "calculate-mcs-tanimoto"
    clusterOptions '--partition "cpu" --time=24:00:00 --mem=64GB --cpus-per-task=32'

    input:
    path(ligand_file_3d)

    output:
    path("mcs_tanimoto"), emit: mcs_tanimoto

    script:
    """
    python3 "${params.scripts}"/calculate_mcs_tanimoto.py --ref-ligand-sdf "${ligand_file_3d}" --output-dir mcs_tanimoto --ncpus 32
    """
}
process CALCULATE_TANIMOTO_COMBO {
    publishDir "${params.chemicalSimilarityData}", mode: 'copy', overwrite: true
    conda "${params.asap}"
    tag "calculate-tanimoto-combo"
    clusterOptions '--partition "cpu" --time=06:00:00 --mem=64GB --cpus-per-task=32'

    input:
    path(ligand_file_3d)

    output:
    path("tanimoto_combo"), emit: tanimoto_combo

    script:
    """
    python3 "${params.scripts}"/calculate_tanimoto_combo.py --ref-ligand-sdf "${ligand_file_3d}" --output-dir tanimoto_combo
    """
}
process COMBINE_CHEMICAL_SIMILARITY_DATA {
    publishDir "${params.chemicalSimilarityData}", mode: 'copy', overwrite: true
    conda "${params.asap}"
    tag "combine-chemical-similarity-data"

    input:
    path csv_files

    output:
    path "combined_chemical_similarity_data.csv", emit: combined_chemical_similarity_data

    script:
    """
    python3 "${params.scripts}/combine_chemical_similarity_data.py" ${csv_files.join(' ')}
    """
}
process CALCULATE_PROTEIN_RMSD_FULL {
    publishDir "${params.proteinRmsdData}/full_protein", mode: 'copy', overwrite: true
    conda "${params.drugforge}"
    tag "full|${ref_id}"

    input:
    tuple val(ref_id), path(ref_json), path(cache_dir, stageAs: 'cache_dir')

    output:
    path("protein_rmsd_full_${ref_id}.csv"), emit: protein_rmsd_full

    script:
    """
    python3 "${params.scripts}"/calculate_protein_RMSD.py \
        --ref_json "${ref_json}" \
        --cache_dir "cache_dir" \
        --output_csv "protein_rmsd_full_${ref_id}.csv"
    """
}
process CALCULATE_PROTEIN_RMSD_BINDING_SITE {
    publishDir "${params.proteinRmsdData}/binding_site", mode: 'copy', overwrite: true
    conda "${params.drugforge}"
    tag "binding_site|${ref_id}"

    input:
    tuple val(ref_id), path(ref_json), path(cache_dir, stageAs: 'cache_dir')

    output:
    path("protein_rmsd_binding_site_${ref_id}.csv"), emit: protein_rmsd_binding_site

    script:
    """
    python3 "${params.scripts}"/calculate_protein_RMSD.py \
        --ref_json "${ref_json}" \
        --cache_dir "cache_dir" \
        --binding_site \
        --output_csv "protein_rmsd_binding_site_${ref_id}.csv"
    """
}
process COMBINE_PROTEIN_RMSD_DATA {
    publishDir "${params.proteinRmsdData}", mode: 'copy', overwrite: true
    conda "${params.asap}"
    tag "combine-protein-rmsd-data"

    input:
    path csv_files

    output:
    path "combined_protein_rmsd.csv", emit: combined_protein_rmsd

    script:
    """
    head -n 1 \$(ls *.csv | head -1) > combined_protein_rmsd.csv
    for f in *.csv; do tail -n +2 "\$f" >> combined_protein_rmsd.csv; done
    """
}
process RUN_BEMIS_MURCKO_CLUSTERING {
    publishDir "${params.chemicalSimilarityData}", mode: 'copy', overwrite: true
    conda "${params.asap}"
    tag "run-bemis-murcko-clustering"

    input:
    path(ligand_file_2d)

    output:
    path "${params.scaffoldDataName}", emit: bemis_murcko_clustering

    script:
    """
    python "${params.scripts}"/run_bemis_murcko_clustering.py --sdf-2d ${ligand_file_2d} --output-dir "${params.scaffoldDataName}"
    """
}