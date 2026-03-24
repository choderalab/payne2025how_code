#!/usr/bin/env nextflow
include {
    GENERATE_DATE_DICTIONARY
    CALCULATE_ECFP_TANIMOTO
    CALCULATE_MCS_TANIMOTO
    CALCULATE_TANIMOTO_COMBO
    COMBINE_CHEMICAL_SIMILARITY_DATA
    RUN_BEMIS_MURCKO_CLUSTERING
    CALCULATE_PROTEIN_RMSD_FULL
    CALCULATE_PROTEIN_RMSD_BINDING_SITE
    COMBINE_PROTEIN_RMSD_DATA
} from "./modules.nf"

workflow {
    // Your input channel for ligand file
    ligand_file_3d = Channel.fromPath("${params.ligandFiles}/${params.ligandFile3d}")
    ligand_file_2d = Channel.fromPath("${params.ligandFiles}/${params.ligandFile2d}")

    // Run calculation processes and collect their CSV outputs
    ecfp_results = CALCULATE_ECFP_TANIMOTO(ligand_file_3d)
        .ecfp_tanimoto
        .map { dir -> dir.listFiles().find { it.name.endsWith('.csv') } }

    mcs_results = CALCULATE_MCS_TANIMOTO(ligand_file_3d)
        .mcs_tanimoto
        .map { dir -> dir.listFiles().find { it.name.endsWith('.csv') } }

    tanimoto_combo_results = CALCULATE_TANIMOTO_COMBO(ligand_file_3d)
        .tanimoto_combo
        .map { dir -> dir.listFiles().find { it.name.endsWith('.csv') } }

    // Combine all CSV files into a single channel
    all_csv_files = ecfp_results
        .mix(mcs_results)
        .mix(tanimoto_combo_results)
        .collect()

    // Pass collected CSV files to combine process
    COMBINE_CHEMICAL_SIMILARITY_DATA(all_csv_files)

    // Generate date dictionary
    GENERATE_DATE_DICTIONARY()

    // Run scaffolding
    RUN_BEMIS_MURCKO_CLUSTERING(ligand_file_2d)
}

// Entry point: Calculate ECFP Tanimoto similarity only
workflow ECFP_ANALYSIS {
    ligand_file_3d = Channel.fromPath("${params.ligandFiles}/${params.ligandFile3d}")
    CALCULATE_ECFP_TANIMOTO(ligand_file_3d)
}

// Entry point: Combine existing CSV files (assumes CSV files already exist)
workflow COMBINE_SIMILARITY_DATA {
    // This assumes CSV files already exist in the expected locations
    // You might need to adjust paths based on your directory structure
    csv_files = Channel.fromPath("${params.chemicalSimilarityData}/*/*.csv").collect()
    COMBINE_CHEMICAL_SIMILARITY_DATA(csv_files)
}

// Entry point: Calculate pairwise protein RMSD for all structures in the fragalysis cache.
// One job is launched per reference structure; each job compares that ref against the
// entire cache directory (--cache_dir), so Nextflow fans out N jobs (one per structure)
// rather than N² jobs.  Both full chain-A and binding-site RMSD are run in parallel.
workflow PROTEIN_RMSD_ANALYSIS {
    cache_dir = Channel.fromPath("${params.fixedFragalysisCachePath}", type: 'dir')

    // One channel item per reference JSON: (ref_id, ref_json, cache_dir)
    ref_inputs = Channel
        .fromPath("${params.fixedFragalysisCachePath}/**/*.json")
        .map { f -> tuple(f.baseName, f) }
        .combine(cache_dir)   // attach the single cache_dir value to every ref

    // Full chain-A RMSD
    full_csvs = CALCULATE_PROTEIN_RMSD_FULL(ref_inputs).protein_rmsd_full

    // Binding-site-only RMSD
    bs_csvs = CALCULATE_PROTEIN_RMSD_BINDING_SITE(ref_inputs).protein_rmsd_binding_site

    // Combine all per-ref CSVs (full + binding-site) into one output.
    // The Binding_Site_Only column in each row distinguishes the two modes.
    all_csvs = full_csvs.mix(bs_csvs).collect()
    COMBINE_PROTEIN_RMSD_DATA(all_csvs)
}

