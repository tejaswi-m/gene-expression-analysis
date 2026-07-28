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

# Study includes MDD, bipolar disorder (BD1/BD2), and healthy control (HC)
# subjects. Only MDD vs HC is scored as Case/Control here; BD samples are
# preprocessed and kept in the metadata but left unlabeled (phenotype blank)
# so they are excluded from the two-group differential-expression/cell-proxy
# steps rather than being folded into either group.
mapping_text = """
GSM973762	S_MDD-376
GSM973763	S_HC-104
GSM973764	S_HC-261
GSM973765	S_MDD-465
GSM973766	S_HC-513
GSM973767	S_BD2-231
GSM973768	S_HC-365
GSM973769	S_HC-3
GSM973770	S_MDD-456
GSM973771	S_HC-284
GSM973772	S_BD2-282
GSM973773	S_MDD-425
GSM973774	S_MDD-397
GSM973775	S_BD1-351
GSM973776	S_MDD-329
GSM973777	S_MDD-395
GSM973778	S_MDD-305
GSM973779	S_HC-45
GSM973780	S_HC-209
GSM973781	S_HC-592
GSM973782	S_BD1-588
GSM973783	S_HC-418
GSM973784	S_MDD-341
GSM973785	S_HC-84
GSM973786	S_HC-296
GSM973787	S_MDD-328
GSM973788	S_HC-172
GSM973789	S_HC-216
GSM973790	S_BD1-362
GSM973791	S_MDD-347
GSM973792	S_HC-321
GSM973793	S_HC-285
GSM973794	S_MDD-343
GSM973795	S_MDD-369
GSM973796	S_MDD-286
GSM973797	S_HC-350
GSM973798	S_MDD-398
GSM973799	S_HC-253
GSM973800	S_BD2-244
GSM973801	S_HC-262
GSM973802	S_MDD-580
GSM973803	S_MDD-407
GSM973804	S_MDD-573
GSM973805	S_MDD-402
GSM973806	S_BD2-338
GSM973807	S_HC-115
GSM973808	S_HC-204
GSM973809	S_BD1-568
GSM973810	S_HC-569
GSM973811	S_HC-279
GSM973812	S_MDD-566
GSM973813	S_HC-88
GSM973814	S_MDD-433
"""


expression_csv, metadata_csv = run_preprocessor(
    mapping_text,
    TEST_DIR / "original_data/GSE39653_series_matrix.txt",
    TEST_DIR / "processed_data",
    annotation_path=TEST_DIR / "original_data" / "GPL10558_annotation.csv",
    mapping_strip_prefix=None,
    sample_title_prefix=None,
    phenotype_column_name="depression_status",
    phenotype_characteristic_name="disease",
    phenotype_value_map={
        "healthy control": "Control",
        "major depressive disorder": "Case",
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
