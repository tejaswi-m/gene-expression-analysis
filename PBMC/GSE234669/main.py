from pathlib import Path
import sys

def _find_project_root(start: Path) -> Path:
    current = start.resolve()
    for parent in [current, *current.parents]:
        if (parent / "preprocessing.py").exists():
            return parent
    raise RuntimeError("Could not locate project root containing preprocessing.py")


PROJECT_ROOT = _find_project_root(Path(__file__).parent)
TEST_DIR = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from preprocessing import run_preprocessor
from analyze_depression_expression import run_expression_analysis
from analyze_cell_population_proxies import run_cell_proxy_analysis
from run_qc_normalize_r import run_qc_normalization

mapping_text = """
GSM7473078	S_control_1
GSM7473079	S_control_2
GSM7473080	S_control_3
GSM7473081	S_control_4
GSM7473082	S_control_5
GSM7473083	S_control_6
GSM7473084	S_control_7
GSM7473085	S_control_8_2
GSM7473086	S_control_9_2
GSM7473087	S_control_10
GSM7473088	S_control_11
GSM7473089	S_control_12_2
GSM7473090	S_ect_ref_1
GSM7473091	S_ect_ref_2
GSM7473092	S_ect_ref_3
GSM7473093	S_ect_ref_4
GSM7473094	S_ect_test_1
GSM7473095	S_ect_test_2
GSM7473096	S_ect_test_3
GSM7473097	S_ect_test_4
GSM7473098	S_ect_test_5
GSM7473099	S_kat_ref_1_2
GSM7473100	S_kat_ref_2
GSM7473101	S_kat_ref_3
GSM7473102	S_kat_ref_4
GSM7473103	S_kat_test_1
GSM7473104	S_kat_test_2
GSM7473105	S_kat_test_3
GSM7473106	S_kat_test_4
GSM7473107	S_sri_ref_1
GSM7473108	S_sri_ref_2
GSM7473109	S_sri_test_1
GSM7473110	S_sri_test_2
GSM7473111	S_sri_test_3
GSM7473112	S_sri_test_4
"""


expression_csv, metadata_csv = run_preprocessor(
    mapping_text,
    TEST_DIR / "original_data/GSE234669_series_matrix.txt",
    TEST_DIR / "processed_data",
    supplementary_tar_path=TEST_DIR / "original_data" / "GSE234669_RAW.tar",
    supplementary_value_format="star_gene_counts",
    supplementary_annotation_path=PROJECT_ROOT / "Ensembl" / "ensembl_gene_id_to_symbol.csv",
    mapping_strip_prefix=None,
    sample_title_prefix=None,
    phenotype_column_name="depression_status",
    phenotype_characteristic_name="treatment",
    phenotype_value_map={
        "None": "Control",
        "electroconvulsive therapy": "Case",
        "ketamine": "Case",
        "ssri": "Case",
    },
)

qc_dir = TEST_DIR / "qc_normalized"
qc_normalized_csv = run_qc_normalization(
    expression_file=TEST_DIR / "processed_data" / "annotated_renamed_unlogged.csv",
    metadata_file=metadata_csv,
    output_dir=qc_dir,
    metadata_sample_column="analysis_sample",
    phenotype_column="depression_status",
    batch_column="batch",
)

results_dir = TEST_DIR / "expression_analysis_results"
run_expression_analysis(
    input_file=qc_normalized_csv,
    output_dir=results_dir,
    metadata_file=TEST_DIR / "processed_data" / "sample_metadata.csv",
    metadata_sample_column="analysis_sample",
    metadata_phenotype_column="depression_status",
    control_label_keywords="control",
    case_label_keywords="case",
    control_label="Control",
    case_label="Case",
)

run_cell_proxy_analysis(
    input_file=results_dir / "deduplicated_expression_dataset.csv",
    metadata_file=TEST_DIR / "processed_data" / "sample_metadata.csv",
    output_dir=TEST_DIR / "cell_population_proxies",
    markers_csv=PROJECT_ROOT / "markers" / "pbmc_cell_markers.csv",
    metadata_sample_column="analysis_sample",
    metadata_phenotype_column="depression_status",
    control_label_keywords="control",
    case_label_keywords="case",
    control_label="Control",
    case_label="Case",
)

print(f"Expression output: {expression_csv}")
print(f"Metadata output:   {metadata_csv}")
print(f"Analysis outputs:  {TEST_DIR / 'cell_population_proxies'}")
