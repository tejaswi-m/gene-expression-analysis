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
GSM936718	S_PBMC_P1-0w_rep1
GSM936719	S_PBMC_P1-8w_rep1
GSM936720	S_PBMC_P2-0w_rep1
GSM936721	S_PBMC_P2-8w_rep1
GSM936722	S_PBMC_P3-0w_rep1
GSM936723	S_PBMC_P3-8w_rep1
GSM936724	S_PBMC_P4-0w_rep1
GSM936725	S_PBMC_P4-8w_rep1
GSM936726	S_PBMC_P5-0w_rep1
GSM936727	S_PBMC_P5-8w_rep1
GSM936728	S_PBMC_P6-0w_rep1
GSM936729	S_PBMC_P6-8w_rep1
GSM936730	S_PBMC_P7-0w_rep1
GSM936731	S_PBMC_P7-8w_rep1
GSM936732	S_PBMC_P8-0w_rep1
GSM936733	S_PBMC_P8-8w_rep1
GSM936734	S_PBMC_P9-0w_rep1
GSM936735	S_PBMC_P9-8w_rep1
GSM936736	S_PBMC_C1-0w_rep1
GSM936737	S_PBMC_C1-8w_rep1
GSM936738	S_PBMC_C2-0w_rep1
GSM936739	S_PBMC_C2-8w_rep1
GSM936740	S_PBMC_C3-0w_rep1
GSM936741	S_PBMC_C3-8w_rep1
GSM936742	S_PBMC_C4-0w_rep1
GSM936743	S_PBMC_C4-8w_rep1
GSM936744	S_PBMC_C5-0w_rep1
GSM936745	S_PBMC_C5-8w_rep1
GSM936746	S_PBMC_C6-0w_rep1
GSM936747	S_PBMC_C6-8w_rep1
GSM936748	S_PBMC_C7-0w_rep1
GSM936749	S_PBMC_C7-8w_rep1
GSM936750	S_PBMC_C8-0w_rep1
GSM936751	S_PBMC_C8-8w_rep1
GSM936752	S_PBMC_C9-0w_rep1
GSM936753	S_PBMC_C9-8w_rep1
"""


expression_csv, metadata_csv = run_preprocessor(
    mapping_text,
    TEST_DIR / "original_data/GSE38206_series_matrix.txt",
    TEST_DIR / "processed_data",
    annotation_path=TEST_DIR / "original_data" / "GPL13607_annotation.csv",
    mapping_strip_prefix=None,
    sample_title_prefix=None,
    phenotype_column_name="depression_status",
    phenotype_characteristic_name="status",
    phenotype_value_map={
        "control": "Control",
        "major depressive episode": "Case",
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
