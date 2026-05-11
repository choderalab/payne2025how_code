#!/usr/bin/env nextflow
include {
    CREATE_EVALUATOR_FACTORY_SETTINGS
    CREATE_EVALUATOR_FACTORY_SETTINGS_CUTOFF
    CREATE_EVALUATOR_FACTORY_SETTINGS_PLIF
    CREATE_EVALUATORS
    RUN_EVALUATORS
    RUN_EVALUATORS_LIGHTWEIGHT
    COMBINE_EVALUATIONS
    CREATE_MULTIPOSE_EVALUATORS
    CREATE_EVALUATORS_MODULAR
    CALCULATE_PLIF_RECALL
    COMBINE_PLIF_RECALL
    MERGE_PLIF_RECALL
} from "./modules.nf"
params.K = 4

workflow RUN_DOCKING_ANALYSIS {
    take:
        name
        docking_results_parquet
        docking_results_json
        evaluator_settings

    main:
        CREATE_EVALUATORS(
            name,
            evaluator_settings,
            docking_results_parquet,
            docking_results_json
        )

        // Create channel from JSON files only after evaluator creation
        eval_inputs_ch = CREATE_EVALUATORS.output.evaluator_json_directory
            .flatMap { dir -> file("${dir}/*.json") }
            .buffer(size: params.K)

        RUN_EVALUATORS(
            name,
            docking_results_parquet,
            docking_results_json,
            eval_inputs_ch,
        )

        // Collect all evaluator results before combining
        all_results = RUN_EVALUATORS.output.evaluator_results
            .flatten()
            .collect()

        COMBINE_EVALUATIONS(
            name,
            all_results
        )
}
dataset_names = [
"posit_single_pose": "ALL_1_poses",
"fred_single_pose": "FRED_1_poses",
"posit_multipose": "ALL_50_poses",
"fred_multipose": "FRED_50_poses",
]

def results = [:]
dataset_names.each { label, name ->
    def parquet = "${params.combinedDockingResultsPath}/${name}.parquet"
    def json = "${params.combinedDockingResultsPath}/${name}.json"

    results[label] = [  // Store directly in map with label as key
        name: name,
        docking_results_parquet: parquet,
        docking_results_json: json
    ]
}

// Print dataset definitions
println "\nDataset Definitions:"
println "===================="
results.each { label, data ->
    println """
    Name: ${label}
    Label: ${data.name}
    Parquet: ${data.docking_results_parquet}
    JSON: ${data.docking_results_json}
    -------------------"""
}

settings_map = [
"datesplit": "reference_split_comparison.yaml",
"x_to_x": "x_to_x_scaffold_split.yaml",
"x_to_x_5": "x_to_x_scaffold_split_5_refs.yaml",
"x_to_y": "x_to_y_scaffold_split.yaml",
"x_to_y_5": "x_to_y_scaffold_split_5_refs.yaml",
"x_to_not_x": "x_to_not_x_scaffold_split.yaml",
"not_x_to_x": "not_x_to_x_scaffold_split.yaml",
"not_x_to_x_5": "not_x_to_x_scaffold_split_5_refs.yaml",
"ecfp4": "increasing_similarity_ecfp4.yaml",
"mcs": "increasing_similarity_mcs.yaml",
"tc": "increasing_similarity_tanimoto_combo_aligned.yaml",]

def settings = [:]
settings_map.each { label, filename ->
    settings[label] = [label: label, filename: "${params.evaluator_configs}/${filename}"]
}

workflow CREATE_EVALUATOR_FACTORY_SETTINGS_WORKFLOW {
    CREATE_EVALUATOR_FACTORY_SETTINGS()
}

workflow RUN_ANALYSIS {
    take:
        result
        setting
    main:
        RUN_DOCKING_ANALYSIS(
            "${result.name}_${setting.label}",
            result.docking_results_parquet,
            result.docking_results_json,
            setting.filename,
        )
}

workflow DATESPLIT_POSIT {RUN_ANALYSIS(results.posit_single_pose, settings.datesplit)}
workflow DATESPLIT_FRED {RUN_ANALYSIS(results.fred_single_pose, settings.datesplit)}
workflow X_TO_X_POSIT {RUN_ANALYSIS(results.posit_single_pose, settings.x_to_x)}
workflow X_TO_X_POSIT_5_REFS {RUN_ANALYSIS(results.posit_single_pose, settings.x_to_x_5)}
workflow X_TO_Y_POSIT {RUN_ANALYSIS(results.posit_single_pose, settings.x_to_y)}
workflow X_TO_Y_POSIT_5_REFS {RUN_ANALYSIS(results.posit_single_pose, settings.x_to_y_5)}
workflow NOT_X_TO_X_POSIT {RUN_ANALYSIS(results.posit_single_pose, settings.not_x_to_x)}
workflow NOT_X_TO_X_POSIT_5_REFS {RUN_ANALYSIS(results.posit_single_pose, settings.not_x_to_x_5)}
workflow X_TO_NOT_X_POSIT {RUN_ANALYSIS(results.posit_single_pose, settings.x_to_not_x)}
workflow INCREASING_SIMILARITY_TC_ALIGNED_POSIT {
    RUN_ANALYSIS(results.posit_single_pose, settings.tc)
}
workflow INCREASING_SIMILARITY_TC_ALIGNED_FRED {
    RUN_ANALYSIS(results.fred_single_pose, settings.tc)
}
workflow INCREASING_SIMILARITY_MCS_POSIT {
    RUN_ANALYSIS(results.posit_single_pose, settings.mcs)
}
workflow INCREASING_SIMILARITY_MCS_FRED {
    RUN_ANALYSIS(results.fred_single_pose, settings.mcs)
}
workflow INCREASING_SIMILARITY_ECFP4_POSIT {
    RUN_ANALYSIS(results.posit_single_pose, settings.ecfp4)
}
workflow INCREASING_SIMILARITY_ECFP4_FRED {
    RUN_ANALYSIS(results.fred_single_pose, settings.ecfp4)
}
workflow INCREASING_SIMILARITY_ECFP4_POSIT_MULTIPOSE {
    RUN_ANALYSIS(results.posit_single_pose, settings.ecfp4)
}
workflow posit_scaffold_splits {
    X_TO_X_POSIT()
    X_TO_X_POSIT_5_REFS()
    X_TO_Y_POSIT()
    X_TO_Y_POSIT_5_REFS()
    NOT_X_TO_X_POSIT()
    NOT_X_TO_X_POSIT_5_REFS()
    X_TO_NOT_X_POSIT()
}
workflow analyze_posit {
    DATESPLIT_POSIT()
    X_TO_X_POSIT()
    X_TO_X_POSIT_5_REFS()
    X_TO_Y_POSIT()
    X_TO_Y_POSIT_5_REFS()
    NOT_X_TO_X_POSIT()
    NOT_X_TO_X_POSIT_5_REFS()
    X_TO_NOT_X_POSIT()
    INCREASING_SIMILARITY_TC_ALIGNED_POSIT()
    INCREASING_SIMILARITY_MCS_POSIT()
    INCREASING_SIMILARITY_ECFP4_POSIT()
    INCREASING_SIMILARITY_ECFP4_POSIT_MULTIPOSE()
}
workflow analyze_fred {
    DATESPLIT_FRED()
    INCREASING_SIMILARITY_TC_ALIGNED_FRED()
    INCREASING_SIMILARITY_MCS_FRED()
    INCREASING_SIMILARITY_ECFP4_FRED()
}
workflow {
    analyze_posit()
    analyze_fred()
}

workflow POSIT_MULTIPOSE_ANALYSIS {
    name = "posit_multipose_analysis"
    CREATE_MULTIPOSE_EVALUATORS(name)

    // Create channel from JSON files only after evaluator creation
    eval_inputs_ch = CREATE_MULTIPOSE_EVALUATORS.output.evaluator_json_directory
        .flatMap { dir -> file("${dir}/*.json") }
        .buffer(size: 1)

    RUN_EVALUATORS(
        name,
        results.posit_multipose.docking_results_parquet,
        results.posit_multipose.docking_results_json,
        eval_inputs_ch,
    )

    // Collect all evaluator results before combining
    all_results = RUN_EVALUATORS.output.evaluator_results
        .flatten()
        .collect()

    COMBINE_EVALUATIONS(
        name,
        all_results
    )
}
// ── RMSD cutoff sensitivity workflows (R3.7) ────────────────────────────────
// Re-run the core POSIT analyses at 1.5 Å and 2.5 Å success thresholds.
// Settings files are generated with _rmsd1.5 / _rmsd2.5 name suffixes so they
// sit alongside the existing 2.0 Å configs in params.evaluator_configs.

def make_settings_map_for_cutoff(String cutoff_label) {
    def suffix = "_rmsd${cutoff_label}"
    return [
        "datesplit"   : "${params.evaluator_configs}/reference_split_comparison${suffix}.yaml",
        "x_to_x"      : "${params.evaluator_configs}/x_to_x_scaffold_split${suffix}.yaml",
        "x_to_x_5"    : "${params.evaluator_configs}/x_to_x_scaffold_split_5_refs${suffix}.yaml",
        "x_to_y"      : "${params.evaluator_configs}/x_to_y_scaffold_split${suffix}.yaml",
        "x_to_y_5"    : "${params.evaluator_configs}/x_to_y_scaffold_split_5_refs${suffix}.yaml",
        "x_to_not_x"  : "${params.evaluator_configs}/x_to_not_x_scaffold_split${suffix}.yaml",
        "not_x_to_x"  : "${params.evaluator_configs}/not_x_to_x_scaffold_split${suffix}.yaml",
        "not_x_to_x_5": "${params.evaluator_configs}/not_x_to_x_scaffold_split_5_refs${suffix}.yaml",
        "ecfp4"        : "${params.evaluator_configs}/increasing_similarity_ecfp4${suffix}.yaml",
        "mcs"          : "${params.evaluator_configs}/increasing_similarity_mcs${suffix}.yaml",
        "tc"           : "${params.evaluator_configs}/increasing_similarity_tanimoto_combo_aligned${suffix}.yaml",
    ]
}

workflow GENERATE_SETTINGS_RMSD_1P5 {
    CREATE_EVALUATOR_FACTORY_SETTINGS_CUTOFF(Channel.value(1.5))
}
workflow GENERATE_SETTINGS_RMSD_2P5 {
    CREATE_EVALUATOR_FACTORY_SETTINGS_CUTOFF(Channel.value(2.5))
}

// 1.5 Å — individual workflows (DSL2 requires one process instance per workflow scope)
def s1p5 = make_settings_map_for_cutoff("1.5")
workflow DATESPLIT_POSIT_RMSD_1P5     { RUN_ANALYSIS(results.posit_single_pose, [label: "datesplit_rmsd1p5",    filename: s1p5.datesplit]) }
workflow X_TO_X_POSIT_RMSD_1P5        { RUN_ANALYSIS(results.posit_single_pose, [label: "x_to_x_rmsd1p5",       filename: s1p5.x_to_x]) }
workflow X_TO_X_5_POSIT_RMSD_1P5      { RUN_ANALYSIS(results.posit_single_pose, [label: "x_to_x_5_rmsd1p5",     filename: s1p5.x_to_x_5]) }
workflow X_TO_Y_POSIT_RMSD_1P5        { RUN_ANALYSIS(results.posit_single_pose, [label: "x_to_y_rmsd1p5",       filename: s1p5.x_to_y]) }
workflow X_TO_Y_5_POSIT_RMSD_1P5      { RUN_ANALYSIS(results.posit_single_pose, [label: "x_to_y_5_rmsd1p5",     filename: s1p5.x_to_y_5]) }
workflow X_TO_NOT_X_POSIT_RMSD_1P5    { RUN_ANALYSIS(results.posit_single_pose, [label: "x_to_not_x_rmsd1p5",   filename: s1p5.x_to_not_x]) }
workflow NOT_X_TO_X_POSIT_RMSD_1P5    { RUN_ANALYSIS(results.posit_single_pose, [label: "not_x_to_x_rmsd1p5",   filename: s1p5.not_x_to_x]) }
workflow NOT_X_TO_X_5_POSIT_RMSD_1P5  { RUN_ANALYSIS(results.posit_single_pose, [label: "not_x_to_x_5_rmsd1p5", filename: s1p5.not_x_to_x_5]) }
workflow ECFP4_POSIT_RMSD_1P5         { RUN_ANALYSIS(results.posit_single_pose, [label: "ecfp4_rmsd1p5",         filename: s1p5.ecfp4]) }
workflow MCS_POSIT_RMSD_1P5           { RUN_ANALYSIS(results.posit_single_pose, [label: "mcs_rmsd1p5",           filename: s1p5.mcs]) }
workflow TC_POSIT_RMSD_1P5            { RUN_ANALYSIS(results.posit_single_pose, [label: "tc_rmsd1p5",            filename: s1p5.tc]) }

workflow analyze_posit_rmsd_1p5 {
    DATESPLIT_POSIT_RMSD_1P5()
    X_TO_X_POSIT_RMSD_1P5()
    X_TO_X_5_POSIT_RMSD_1P5()
    X_TO_Y_POSIT_RMSD_1P5()
    X_TO_Y_5_POSIT_RMSD_1P5()
    X_TO_NOT_X_POSIT_RMSD_1P5()
    NOT_X_TO_X_POSIT_RMSD_1P5()
    NOT_X_TO_X_5_POSIT_RMSD_1P5()
    ECFP4_POSIT_RMSD_1P5()
    MCS_POSIT_RMSD_1P5()
    TC_POSIT_RMSD_1P5()
}

// 2.5 Å — individual workflows
def s2p5 = make_settings_map_for_cutoff("2.5")
workflow DATESPLIT_POSIT_RMSD_2P5     { RUN_ANALYSIS(results.posit_single_pose, [label: "datesplit_rmsd2p5",    filename: s2p5.datesplit]) }
workflow X_TO_X_POSIT_RMSD_2P5        { RUN_ANALYSIS(results.posit_single_pose, [label: "x_to_x_rmsd2p5",       filename: s2p5.x_to_x]) }
workflow X_TO_X_5_POSIT_RMSD_2P5      { RUN_ANALYSIS(results.posit_single_pose, [label: "x_to_x_5_rmsd2p5",     filename: s2p5.x_to_x_5]) }
workflow X_TO_Y_POSIT_RMSD_2P5        { RUN_ANALYSIS(results.posit_single_pose, [label: "x_to_y_rmsd2p5",       filename: s2p5.x_to_y]) }
workflow X_TO_Y_5_POSIT_RMSD_2P5      { RUN_ANALYSIS(results.posit_single_pose, [label: "x_to_y_5_rmsd2p5",     filename: s2p5.x_to_y_5]) }
workflow X_TO_NOT_X_POSIT_RMSD_2P5    { RUN_ANALYSIS(results.posit_single_pose, [label: "x_to_not_x_rmsd2p5",   filename: s2p5.x_to_not_x]) }
workflow NOT_X_TO_X_POSIT_RMSD_2P5    { RUN_ANALYSIS(results.posit_single_pose, [label: "not_x_to_x_rmsd2p5",   filename: s2p5.not_x_to_x]) }
workflow NOT_X_TO_X_5_POSIT_RMSD_2P5  { RUN_ANALYSIS(results.posit_single_pose, [label: "not_x_to_x_5_rmsd2p5", filename: s2p5.not_x_to_x_5]) }
workflow ECFP4_POSIT_RMSD_2P5         { RUN_ANALYSIS(results.posit_single_pose, [label: "ecfp4_rmsd2p5",         filename: s2p5.ecfp4]) }
workflow MCS_POSIT_RMSD_2P5           { RUN_ANALYSIS(results.posit_single_pose, [label: "mcs_rmsd2p5",           filename: s2p5.mcs]) }
workflow TC_POSIT_RMSD_2P5            { RUN_ANALYSIS(results.posit_single_pose, [label: "tc_rmsd2p5",            filename: s2p5.tc]) }

workflow analyze_posit_rmsd_2p5 {
    DATESPLIT_POSIT_RMSD_2P5()
    X_TO_X_POSIT_RMSD_2P5()
    X_TO_X_5_POSIT_RMSD_2P5()
    X_TO_Y_POSIT_RMSD_2P5()
    X_TO_Y_5_POSIT_RMSD_2P5()
    X_TO_NOT_X_POSIT_RMSD_2P5()
    NOT_X_TO_X_POSIT_RMSD_2P5()
    NOT_X_TO_X_5_POSIT_RMSD_2P5()
    ECFP4_POSIT_RMSD_2P5()
    MCS_POSIT_RMSD_2P5()
    TC_POSIT_RMSD_2P5()
}
// ── end RMSD cutoff sensitivity ──────────────────────────────────────────────

workflow FRED_MULTIPOSE_ANALYSIS {
    name = "fred_multipose_analysis"
    CREATE_MULTIPOSE_EVALUATORS(name)

    // Create channel from JSON files only after evaluator creation
    eval_inputs_ch = CREATE_MULTIPOSE_EVALUATORS.output.evaluator_json_directory
        .flatMap { dir -> file("${dir}/*.json") }
        .buffer(size: 1)

    RUN_EVALUATORS(
        name,
        results.fred_multipose.docking_results_parquet,
        results.fred_multipose.docking_results_json,
        eval_inputs_ch,
    )

    // Collect all evaluator results before combining
    all_results = RUN_EVALUATORS.output.evaluator_results
        .flatten()
        .collect()

    COMBINE_EVALUATIONS(
        name,
        all_results
    )
}
workflow SCAFFOLD_DATE_SPLIT {
    name = "posit_scaffold_date_split"
    script_path = "${params.scripts}/create_evaluators_scaffold_datesplit.py"

    CREATE_EVALUATORS_MODULAR(
        name,
        script_path,
        results.posit_single_pose.docking_results_parquet,
        results.posit_single_pose.docking_results_json,
    )

    // Create channel from JSON files only after evaluator creation
    eval_inputs_ch = CREATE_EVALUATORS_MODULAR.output.evaluator_json_directory
        .flatMap { dir -> file("${dir}/*.json") }
        .buffer(size: 1)

    RUN_EVALUATORS_LIGHTWEIGHT(
        name,
        results.posit_single_pose.docking_results_parquet,
        results.posit_single_pose.docking_results_json,
        eval_inputs_ch,
    )

    // Collect all evaluator results before combining
    all_results = RUN_EVALUATORS_LIGHTWEIGHT.output.evaluator_results
        .flatten()
        .collect()

    COMBINE_EVALUATIONS(
        name,
        all_results
    )
}
workflow POSIT_REVERSE_SIMILARITY_SPLIT {
    name = "posit_reverse_similarity_split"
    script_path = "${params.scripts}/create_reverse_similarity_split_evaluators.py"

    CREATE_EVALUATORS_MODULAR(
        name,
        script_path,
        results.posit_single_pose.docking_results_parquet,
        results.posit_single_pose.docking_results_json,
    )

    // Create channel from JSON files only after evaluator creation
    eval_inputs_ch = CREATE_EVALUATORS_MODULAR.output.evaluator_json_directory
        .flatMap { dir -> file("${dir}/*.json") }
        .buffer(size: 1)

    RUN_EVALUATORS_LIGHTWEIGHT(
        name,
        results.posit_single_pose.docking_results_parquet,
        results.posit_single_pose.docking_results_json,
        eval_inputs_ch,
    )

    // Collect all evaluator results before combining
    all_results = RUN_EVALUATORS_LIGHTWEIGHT.output.evaluator_results
        .flatten()
        .collect()

    COMBINE_EVALUATIONS(
        name,
        all_results
    )
}
// ── PLIF Recall analysis (R1.2, R2.2) ───────────────────────────────────────
workflow PLIF_RECALL_POSIT {
    sdfs_ch = Channel.fromPath(
        "${params.dockedFiles}/ALL_1_poses/*/docking_results.sdf"
    ).map { sdf -> tuple(sdf.parent.name, sdf) }

    CALCULATE_PLIF_RECALL(sdfs_ch)

    COMBINE_PLIF_RECALL(
        CALCULATE_PLIF_RECALL.output.plif_recall_csv
            .collect()
    )
}
// ── end PLIF Recall ──────────────────────────────────────────────────────────

// ── PLIF Recall analysis — evaluator pipeline (runs after PLIF_RECALL_POSIT) ─
// Merges plif_recall_combined.csv into the POSIT single-pose parquet, generates
// evaluator configs, and runs the full analysis with PLIF recall as the success
// metric.  Run -entry GENERATE_SETTINGS_PLIF first, then PLIF_RECALL_ANALYSIS.

def make_plif_settings_map(String cutoff_label) {
    def suffix = "_plif${cutoff_label}"
    return [
        "datesplit"   : "${params.evaluator_configs}/reference_split_comparison${suffix}.yaml",
        "x_to_x"      : "${params.evaluator_configs}/x_to_x_scaffold_split${suffix}.yaml",
        "x_to_x_5"    : "${params.evaluator_configs}/x_to_x_scaffold_split_5_refs${suffix}.yaml",
        "x_to_y"      : "${params.evaluator_configs}/x_to_y_scaffold_split${suffix}.yaml",
        "x_to_y_5"    : "${params.evaluator_configs}/x_to_y_scaffold_split_5_refs${suffix}.yaml",
        "x_to_not_x"  : "${params.evaluator_configs}/x_to_not_x_scaffold_split${suffix}.yaml",
        "not_x_to_x"  : "${params.evaluator_configs}/not_x_to_x_scaffold_split${suffix}.yaml",
        "not_x_to_x_5": "${params.evaluator_configs}/not_x_to_x_scaffold_split_5_refs${suffix}.yaml",
        "ecfp4"        : "${params.evaluator_configs}/increasing_similarity_ecfp4${suffix}.yaml",
        "mcs"          : "${params.evaluator_configs}/increasing_similarity_mcs${suffix}.yaml",
        "tc"           : "${params.evaluator_configs}/increasing_similarity_tanimoto_combo_aligned${suffix}.yaml",
    ]
}

workflow GENERATE_SETTINGS_PLIF {
    CREATE_EVALUATOR_FACTORY_SETTINGS_PLIF(Channel.value(0.5))
}

def sp5 = make_plif_settings_map("0.5")

workflow PLIF_MERGE_POSIT {
    plif_csv = Channel.fromPath("${params.evaluationResults}/plif_recall_combined.csv")
    MERGE_PLIF_RECALL(
        "ALL_1_poses",
        results.posit_single_pose.docking_results_parquet,
        results.posit_single_pose.docking_results_json,
        plif_csv,
    )
}

def plif_result = [
    name: "ALL_1_poses_plif",
    docking_results_parquet: "${params.combinedDockingResultsPath}/ALL_1_poses_plif.parquet",
    docking_results_json:    "${params.combinedDockingResultsPath}/ALL_1_poses_plif.json",
]

workflow PLIF_DATESPLIT_POSIT     { RUN_ANALYSIS(plif_result, [label: "datesplit_plif0.5",    filename: sp5.datesplit]) }
workflow PLIF_X_TO_X_POSIT        { RUN_ANALYSIS(plif_result, [label: "x_to_x_plif0.5",        filename: sp5.x_to_x]) }
workflow PLIF_X_TO_X_5_POSIT      { RUN_ANALYSIS(plif_result, [label: "x_to_x_5_plif0.5",      filename: sp5.x_to_x_5]) }
workflow PLIF_X_TO_Y_POSIT        { RUN_ANALYSIS(plif_result, [label: "x_to_y_plif0.5",        filename: sp5.x_to_y]) }
workflow PLIF_X_TO_Y_5_POSIT      { RUN_ANALYSIS(plif_result, [label: "x_to_y_5_plif0.5",      filename: sp5.x_to_y_5]) }
workflow PLIF_NOT_X_TO_X_POSIT    { RUN_ANALYSIS(plif_result, [label: "not_x_to_x_plif0.5",    filename: sp5.not_x_to_x]) }
workflow PLIF_NOT_X_TO_X_5_POSIT  { RUN_ANALYSIS(plif_result, [label: "not_x_to_x_5_plif0.5",  filename: sp5.not_x_to_x_5]) }
workflow PLIF_X_TO_NOT_X_POSIT    { RUN_ANALYSIS(plif_result, [label: "x_to_not_x_plif0.5",    filename: sp5.x_to_not_x]) }
workflow PLIF_TC_POSIT            { RUN_ANALYSIS(plif_result, [label: "tc_plif0.5",            filename: sp5.tc]) }
workflow PLIF_MCS_POSIT           { RUN_ANALYSIS(plif_result, [label: "mcs_plif0.5",           filename: sp5.mcs]) }
workflow PLIF_ECFP4_POSIT         { RUN_ANALYSIS(plif_result, [label: "ecfp4_plif0.5",         filename: sp5.ecfp4]) }

workflow analyze_posit_plif {
    PLIF_DATESPLIT_POSIT()
    PLIF_X_TO_X_POSIT()
    PLIF_X_TO_X_5_POSIT()
    PLIF_X_TO_Y_POSIT()
    PLIF_X_TO_Y_5_POSIT()
    PLIF_NOT_X_TO_X_POSIT()
    PLIF_NOT_X_TO_X_5_POSIT()
    PLIF_X_TO_NOT_X_POSIT()
    PLIF_TC_POSIT()
    PLIF_MCS_POSIT()
    PLIF_ECFP4_POSIT()
}
// ── end PLIF Recall analysis ─────────────────────────────────────────────────
