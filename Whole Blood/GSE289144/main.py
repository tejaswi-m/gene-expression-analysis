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

# This series (Whole Blood) contains only unaffected controls (n=37) -- the
# study's affected cohort was deposited separately as GSE289146 [PBMC], a
# different tissue, so no case-vs-control comparison is run here.
from preprocessing import run_preprocessor
from run_qc_normalize_r import run_qc_normalization

mapping_text = """
GSM8784978	S_23L002656_S8
GSM8784979	S_23L002669_S21
GSM8784980	S_23L002731_S83
GSM8784981	S_23L002650_S2
GSM8784982	S_23L002716_S68
GSM8784983	S_23L002723_S75
GSM8784984	S_23L002689_S41
GSM8784985	S_23L002714_S66
GSM8784986	S_23L002691_S43
GSM8784987	S_23L002734_S86
GSM8784988	S_23L002683_S35
GSM8784989	S_23L002681_S33
GSM8784990	S_23L002702_S54
GSM8784991	S_23L002685_S37
GSM8784992	S_23L002738_S90
GSM8784993	S_23L002721_S73
GSM8784994	S_23L002660_S12
GSM8784995	S_23L002695_S47
GSM8784996	S_23L002730_S82
GSM8784997	S_23L002719_S71
GSM8784998	S_23L002677_S29
GSM8784999	S_23L002742_S94
GSM8785000	S_23L002662_S14
GSM8785001	S_23L002664_S16
GSM8785002	S_23L002727_S79
GSM8785003	S_23L002671_S23
GSM8785004	S_23L002709_S61
GSM8785005	S_23L002680_S32
GSM8785006	S_23L002688_S40
GSM8785007	S_23L002717_S69
GSM8785008	S_23L002735_S87
GSM8785009	S_23L002725_S77
GSM8785010	S_23L002653_S5
GSM8785011	S_23L002667_S19
GSM8785012	S_23L002694_S46
GSM8785013	S_23L002700_S52
GSM8785014	S_23L002676_S28
"""


expression_csv, metadata_csv = run_preprocessor(
    mapping_text,
    TEST_DIR / "original_data/GSE289144_series_matrix.txt",
    TEST_DIR / "processed_data",
    supplementary_tar_path=TEST_DIR / "original_data" / "GSE289144_RAW.tar",
    supplementary_value_format="rsem_gene_results",
    supplementary_value_column="TPM",
    supplementary_annotation_path=PROJECT_ROOT / "Ensembl" / "ensembl_gene_id_to_symbol.csv",
    mapping_strip_prefix=None,
    sample_title_prefix=None,
    phenotype_column_name="depression_status",
    phenotype_characteristic_name="status",
    phenotype_value_map={
        "unaffected": "Control",
        "affected": "Case",
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

print(f"Expression output:     {expression_csv}")
print(f"Metadata output:       {metadata_csv}")
print(f"QC/normalized output:  {qc_normalized_csv}")
print(
    "No differential-expression or cell-proxy analysis is run: this series "
    "has no case samples (see comment above)."
)
