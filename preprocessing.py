#!/usr/bin/env python3
"""Build processed expression and sample metadata tables from GEO series data.

This script consolidates the existing preprocessing steps into one entrypoint:
1. Load the original GEO series matrix table from the raw TXT file.
2. Merge probe-level gene annotation.
3. Rename GSM sample columns to control_*/case_* labels.
4. Convert mapped sample columns from log2 scale back to linear scale.
5. Extract sample-level metadata from the series matrix header.
6. Write the processed expression CSV and metadata CSV.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import ssl
import subprocess
import tarfile
import textwrap
from io import StringIO
from pathlib import Path
from typing import Callable
from urllib.error import URLError
from urllib.request import urlopen

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent
MAPPING_TEXT = ""


def parse_series_platform_id(series_matrix_path: Path) -> str | None:
    with series_matrix_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Series_platform_id"):
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    return None
                value = parts[1].strip().strip('"')
                return value or None
            if line.startswith("!series_matrix_table_begin"):
                break
    return None


def resolve_annotation_path(
    series_matrix_path: Path,
    explicit_annotation_path: str | Path | None,
) -> Path | None:
    if explicit_annotation_path:
        candidate = Path(explicit_annotation_path)
        return candidate if candidate.exists() else None

    platform_id = parse_series_platform_id(series_matrix_path)
    search_roots = [PROJECT_ROOT]

    if platform_id:
        platform_key = platform_id.lower().replace("-", "")
        direct_candidates = [
            PROJECT_ROOT / platform_id / f"{platform_id}_annotation.csv",
            PROJECT_ROOT / platform_id / "annotation.csv",
            PROJECT_ROOT / platform_id / "platform_annotation.csv",
            PROJECT_ROOT / f"{platform_id}_annotation.csv",
            PROJECT_ROOT / f"{platform_id}.csv",
        ]
        for candidate in direct_candidates:
            if candidate.exists():
                return candidate

        for root in search_roots:
            for file_path in root.rglob("*.csv"):
                path_lower = str(file_path).lower()
                name_lower = file_path.name.lower()
                if platform_key not in path_lower:
                    continue
                if "annotation" in name_lower:
                    return file_path

        return None
    return None


def _choose_column(columns: list[str], candidates: list[str]) -> str | None:
    lower_map = {col.lower().strip(): col for col in columns}
    for candidate in candidates:
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]
    return None


def _fetch_text_from_url(url: str, timeout: int = 90) -> str | None:
    try:
        with urlopen(url, timeout=timeout) as response:
            return response.read().decode("utf-8", errors="replace")
    except Exception:
        pass

    # Some local Python SSL setups fail validation for NCBI endpoints.
    try:
        insecure_ctx = ssl._create_unverified_context()
        with urlopen(url, timeout=timeout, context=insecure_ctx) as response:
            return response.read().decode("utf-8", errors="replace")
    except Exception:
        pass

    # Final fallback that works in environments where urllib TLS fails but curl succeeds.
    try:
        result = subprocess.run(
            ["curl", "-L", "--max-time", str(timeout), url],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode == 0 and result.stdout:
            return result.stdout
    except Exception:
        pass

    return None


def download_platform_annotation_from_geo(platform_id: str, output_path: Path) -> Path | None:
    # GEO text view contains a platform table between !platform_table_begin/end.
    url = (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        f"?acc={platform_id}&targ=self&form=text&view=full"
    )
    text = _fetch_text_from_url(url, timeout=90)
    if not text:
        return None

    begin_token = "!platform_table_begin"
    end_token = "!platform_table_end"
    begin = text.find(begin_token)
    end = text.find(end_token)
    if begin == -1 or end == -1 or end <= begin:
        return None

    table_text = text[begin + len(begin_token):end].strip()
    if not table_text:
        return None

    platform_df = pd.read_csv(
        StringIO(table_text),
        sep="\t",
        dtype=str,
        engine="python",
    )
    if platform_df.empty:
        return None

    cols = list(platform_df.columns)
    id_col = _choose_column(cols, ["ID", "ID_REF", "PROBEID", "Probe ID"])
    symbol_col = _choose_column(
        cols,
        [
            "Gene Symbol",
            "GENE_SYMBOL",
            "GENE SYMBOL",
            "Symbol",
            "Gene symbol",
            "Gene",
            "GENE",
            "GENE NAME",
        ],
    )
    title_col = _choose_column(
        cols,
        [
            "Gene Title",
            "GENE_TITLE",
            "GENE TITLE",
            "Gene title",
            "Description",
            "DESCRIPTION",
            "Gene Description",
            "Title",
            "GENE_NAME",
        ],
    )
    if id_col is None:
        return None

    out_df = pd.DataFrame({
        "ID": platform_df[id_col].astype(str),
        "Gene Symbol": platform_df[symbol_col].astype(str) if symbol_col else "",
        "Gene Title": platform_df[title_col].astype(str) if title_col else "",
    })
    out_df["ID"] = out_df["ID"].str.strip()
    out_df["Gene Symbol"] = out_df["Gene Symbol"].str.strip()
    out_df["Gene Title"] = out_df["Gene Title"].str.strip()
    out_df = out_df[out_df["ID"] != ""]
    out_df = out_df.drop_duplicates(subset=["ID"], keep="first")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(output_path, index=False)
    return output_path

def resolve_mapping_text(args: argparse.Namespace) -> str:
    if args.mapping_text and args.mapping_file:
        raise ValueError("Use only one of --mapping-text or --mapping-file")
    if args.mapping_text:
        return textwrap.dedent(args.mapping_text).strip()
    if args.mapping_file:
        mapping_path = Path(args.mapping_file)
        ensure_exists(mapping_path, "Mapping file")
        return mapping_path.read_text(encoding="utf-8").strip()
    return MAPPING_TEXT.strip()


def find_table_bounds(file_path: Path) -> tuple[int, int]:
    begin_line = None
    end_line = None

    with file_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, start=1):
            if line.startswith("!series_matrix_table_begin"):
                begin_line = line_no
            elif line.startswith("!series_matrix_table_end"):
                end_line = line_no
                break

    if begin_line is None or end_line is None:
        raise ValueError("Could not find series matrix table boundaries")
    if end_line <= begin_line + 1:
        raise ValueError("No table content found in the series matrix file")

    return begin_line, end_line


def load_series_matrix_as_dataframe(file_path: Path) -> pd.DataFrame:
    begin_line, end_line = find_table_bounds(file_path)
    df = pd.read_csv(
        file_path,
        sep="\t",
        skiprows=begin_line,
        nrows=end_line - begin_line - 1,
        dtype=str,
        engine="python",
    )
    return df.dropna(axis=1, how="all")


def drop_rows_with_only_reference(
    expression_df: pd.DataFrame,
    reference_column: str = "ID_REF",
) -> pd.DataFrame:
    if reference_column not in expression_df.columns:
        return expression_df

    non_ref_cols = [col for col in expression_df.columns if col != reference_column]
    if not non_ref_cols:
        return expression_df

    non_ref_values = expression_df[non_ref_cols].replace(r"^\s*$", np.nan, regex=True)
    has_any_signal = non_ref_values.notna().any(axis=1)
    return expression_df.loc[has_any_signal].reset_index(drop=True)


def load_annotation(
    annotation_path: Path,
    annotation_id_column: str,
    annotation_symbol_column: str,
    annotation_title_column: str,
) -> pd.DataFrame:
    annotation_df = pd.read_csv(
        annotation_path,
        dtype=str,
        low_memory=False,
        usecols=[
            annotation_id_column,
            annotation_symbol_column,
            annotation_title_column,
        ],
    )
    return annotation_df


def merge_annotation(
    expression_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    expression_id_column: str,
    annotation_id_column: str,
    annotation_symbol_column: str,
    annotation_title_column: str,
) -> pd.DataFrame:
    merged_df = expression_df.merge(
        annotation_df,
        left_on=expression_id_column,
        right_on=annotation_id_column,
        how="left",
    ).drop(columns=[annotation_id_column])

    columns = list(merged_df.columns)
    columns.remove(annotation_symbol_column)
    columns.remove(annotation_title_column)
    insert_position = columns.index(expression_id_column) + 1
    ordered_columns = (
        columns[:insert_position]
        + [annotation_symbol_column, annotation_title_column]
        + columns[insert_position:]
    )
    return merged_df[ordered_columns]


def build_column_mapping(mapping_text: str, strip_prefix: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for line in mapping_text.strip().splitlines():
        if not line.strip():
            continue
        parts = line.strip().split(maxsplit=1)
        if len(parts) != 2:
            raise ValueError(
                "Each mapping line must contain a GSM accession followed by a sample label"
            )
        gsm_accession, label = parts
        cleaned_label = label.removeprefix(strip_prefix) if strip_prefix else label
        mapping[gsm_accession] = cleaned_label
    return mapping


def rename_and_unlog(expression_df: pd.DataFrame, mapping: dict[str, str]) -> tuple[pd.DataFrame, list[str]]:
    present_mapping = {gsm: label for gsm, label in mapping.items() if gsm in expression_df.columns}
    if not present_mapping:
        raise ValueError("No mapped GSM columns were found in the expression matrix")

    renamed_df = expression_df.rename(columns=present_mapping)
    mapped_columns = list(present_mapping.values())
    numeric_values = renamed_df[mapped_columns].apply(pd.to_numeric, errors="coerce")
    renamed_df[mapped_columns] = np.power(2.0, numeric_values)
    return renamed_df, mapped_columns


def clean_gene_symbols(series: pd.Series) -> pd.Series:
    cleaned = series.fillna("").astype(str).str.strip()
    return cleaned.replace({"nan": "", "---": ""})


def deduplicate_probes_by_mean_signal(
    expression_df: pd.DataFrame,
    gene_symbol_column: str,
    sample_columns: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    working = expression_df.copy()
    if gene_symbol_column not in working.columns:
        raise ValueError(f"Gene identifier column not found: {gene_symbol_column}")

    working[gene_symbol_column] = clean_gene_symbols(working[gene_symbol_column])
    if (working[gene_symbol_column] != "").sum() == 0 and "ID_REF" in working.columns:
        working[gene_symbol_column] = working["ID_REF"].astype(str)

    working = working[working[gene_symbol_column] != ""].copy()

    numeric_samples = working[sample_columns].apply(pd.to_numeric, errors="coerce")
    working = pd.concat([working.drop(columns=sample_columns), numeric_samples], axis=1)
    working = working.assign(
        probe_mean_signal=numeric_samples.mean(axis=1, skipna=True)
    )
    working = working.sort_values(
        by=[gene_symbol_column, "probe_mean_signal"],
        ascending=[True, False],
    )

    duplicates = (
        working.groupby(gene_symbol_column, as_index=False)
        .agg(
            probes_for_gene=(gene_symbol_column, "size"),
            retained_probe_mean_signal=("probe_mean_signal", "first"),
        )
        .sort_values("probes_for_gene", ascending=False)
    )

    deduplicated = (
        working.drop_duplicates(subset=gene_symbol_column, keep="first")
        .drop(columns=["probe_mean_signal"])
        .sort_values(gene_symbol_column)
        .reset_index(drop=True)
    )
    return deduplicated, duplicates


def read_sample_header_rows(series_matrix_path: Path) -> dict[str, list[list[str]]]:
    sample_rows: dict[str, list[list[str]]] = {}
    with series_matrix_path.open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if row[0] == "!series_matrix_table_begin":
                break
            if row[0].startswith("!Sample_"):
                sample_rows.setdefault(row[0], []).append(row[1:])
    return sample_rows


def strip_quotes(value: str) -> str:
    return value.strip().strip('"')


def parse_characteristic_rows(rows: list[list[str]], sample_count: int) -> list[dict[str, str]]:
    metadata = [dict() for _ in range(sample_count)]
    for row in rows:
        values = [strip_quotes(value) for value in row]
        if len(values) != sample_count:
            raise ValueError("Characteristic row length does not match sample count")
        for index, value in enumerate(values):
            if ": " in value:
                key, parsed_value = value.split(": ", 1)
                metadata[index][key.strip()] = parsed_value.strip()
            else:
                metadata[index]["characteristic"] = value
    return metadata


def infer_group_label(
    sample_name: str,
    control_prefix: str | None,
    case_prefix: str | None,
    control_label: str | None,
    case_label: str | None,
) -> str:
    if control_prefix and control_label and sample_name.startswith(control_prefix):
        return control_label
    if case_prefix and case_label and sample_name.startswith(case_prefix):
        return case_label
    return ""


def build_metadata_frame(
    sample_rows: dict[str, list[list[str]]],
    mapping: dict[str, str],
    sample_title_prefix: str | None,
    phenotype_column_name: str | None,
    control_prefix: str | None,
    case_prefix: str | None,
    control_label: str | None,
    case_label: str | None,
    phenotype_characteristic_name: str | None = None,
    phenotype_value_map: dict[str, str] | None = None,
) -> pd.DataFrame:
    required_rows = [
        "!Sample_title",
        "!Sample_geo_accession",
        "!Sample_source_name_ch1",
    ]
    for key in required_rows:
        if key not in sample_rows or len(sample_rows[key]) != 1:
            raise ValueError(f"Expected exactly one row for {key}")

    titles = [strip_quotes(value) for value in sample_rows["!Sample_title"][0]]
    accessions = [strip_quotes(value) for value in sample_rows["!Sample_geo_accession"][0]]
    sources = [strip_quotes(value) for value in sample_rows["!Sample_source_name_ch1"][0]]
    sample_count = len(titles)

    if len(accessions) != sample_count or len(sources) != sample_count:
        raise ValueError("Sample metadata rows are not aligned")

    characteristics = parse_characteristic_rows(
        sample_rows.get("!Sample_characteristics_ch1", []),
        sample_count,
    )

    records = []
    for index, title in enumerate(titles):
        analysis_sample = mapping.get(accessions[index])
        if analysis_sample is None:
            analysis_sample = title.removeprefix(sample_title_prefix) if sample_title_prefix else title
        sample_characteristics = characteristics[index]

        if phenotype_characteristic_name and phenotype_value_map:
            phenotype_label = infer_group_label_from_characteristic(
                sample_characteristics.get(phenotype_characteristic_name, ""),
                phenotype_value_map,
            )
        else:
            phenotype_label = infer_group_label(
                analysis_sample,
                control_prefix,
                case_prefix,
                control_label,
                case_label,
            )

        record = {
            "analysis_sample": analysis_sample,
            "geo_sample_title": title,
            "gsm_accession": accessions[index],
            "source_name": sources[index],
            "subject_group": sample_characteristics.get("subject group", ""),
            (phenotype_column_name or "phenotype_group"): phenotype_label,
        }
        record.update(sample_characteristics)
        records.append(record)

    metadata_df = pd.DataFrame(records)
    if "age" in metadata_df.columns:
        metadata_df["age"] = pd.to_numeric(metadata_df["age"], errors="coerce")
    return metadata_df.sort_values("analysis_sample").reset_index(drop=True)


def validate_sample_alignment(
    metadata_df: pd.DataFrame,
    expression_df: pd.DataFrame,
    mapped_columns: list[str],
) -> None:
    expression_samples = [column for column in mapped_columns if column in expression_df.columns]
    metadata_samples = metadata_df["analysis_sample"].tolist()

    missing_in_metadata = sorted(set(expression_samples) - set(metadata_samples))
    missing_in_expression = sorted(set(metadata_samples) - set(expression_samples))
    if missing_in_metadata or missing_in_expression:
        raise ValueError(
            "Sample mapping mismatch. "
            f"Missing in metadata: {missing_in_metadata[:5]}; "
            f"Missing in expression: {missing_in_expression[:5]}"
        )


def ensure_exists(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{description} not found: {path}")


def generate_fallback_annotation(
    series_matrix_path: Path,
    output_path: Path,
    expression_id_column: str = "ID_REF",
    annotation_id_column: str = "ID",
    annotation_symbol_column: str = "Gene Symbol",
    annotation_title_column: str = "Gene Title",
) -> Path:
    """Create a minimal probe annotation table from expression IDs.

    This allows pipelines to run when platform annotation files are unavailable.
    Gene symbol defaults to the probe ID and title is left blank.
    """

    expression_df = load_series_matrix_as_dataframe(series_matrix_path)
    if expression_id_column not in expression_df.columns:
        raise ValueError(
            f"Cannot generate fallback annotation: missing {expression_id_column} column"
        )

    probe_ids = (
        expression_df[expression_id_column]
        .dropna()
        .astype(str)
        .str.strip()
    )
    probe_ids = probe_ids[probe_ids != ""]
    probe_ids = probe_ids.drop_duplicates().reset_index(drop=True)

    fallback_annotation = pd.DataFrame(
        {
            annotation_id_column: probe_ids,
            annotation_symbol_column: probe_ids,
            annotation_title_column: "",
        }
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fallback_annotation.to_csv(output_path, index=False)
    return output_path


def parse_supplementary_file_map(series_matrix_path: Path) -> dict[str, str]:
    """Map each GSM accession to the basename of its per-sample supplementary file.

    Used for series (e.g. RNA-seq) where per-gene values are not embedded in the
    series matrix table itself but shipped as one file per sample in a RAW.tar.
    """
    sample_rows = read_sample_header_rows(series_matrix_path)
    accessions = [strip_quotes(value) for value in sample_rows["!Sample_geo_accession"][0]]
    supplementary_key = next(
        (key for key in sample_rows if key.startswith("!Sample_supplementary_file")),
        None,
    )
    if supplementary_key is None:
        raise ValueError("No !Sample_supplementary_file rows found in series matrix")

    filenames = [strip_quotes(value) for value in sample_rows[supplementary_key][0]]
    if len(filenames) != len(accessions):
        raise ValueError("Supplementary file row is not aligned with sample accessions")

    return {
        accession: filename.rsplit("/", 1)[-1]
        for accession, filename in zip(accessions, filenames)
    }


def extract_supplementary_archive(tar_path: Path, extract_dir: Path) -> Path:
    extract_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r") as archive:
        archive.extractall(extract_dir)
    return extract_dir


def load_gene_value_file(file_path: Path) -> pd.Series:
    """Read a headerless two-column (gene_symbol, value) per-sample file."""
    opener = gzip.open if file_path.suffix == ".gz" else open
    with opener(file_path, "rt") as handle:
        series = pd.read_csv(
            handle,
            sep="\t",
            header=None,
            names=["Gene Symbol", "value"],
            dtype={"Gene Symbol": str},
        ).set_index("Gene Symbol")["value"]
    return series.groupby(level=0).mean()


def load_rsem_gene_results_file(file_path: Path, value_column: str = "TPM") -> pd.Series:
    """Read an RSEM *.genes.results file and return log2(value + 1) per gene.

    The log2 transform mirrors the scale rename_and_unlog() expects (it un-logs
    with 2**x), keeping this format interchangeable with the plain two-column
    gene-value files handled by load_gene_value_file().
    """
    opener = gzip.open if file_path.suffix == ".gz" else open
    with opener(file_path, "rt") as handle:
        results_df = pd.read_csv(handle, sep="\t", dtype=str)
    values = pd.to_numeric(results_df[value_column], errors="coerce")
    series = pd.Series(values.to_numpy(), index=results_df["gene_id"])
    series = series.groupby(level=0).sum()
    return np.log2(series + 1.0)


def load_star_gene_counts_file(file_path: Path) -> pd.Series:
    """Read a STAR ReadsPerGene.out.tab-style per-sample file (header row,
    gene_id + raw count, plus leading N_unmapped/N_multimapping/N_noFeature/
    N_ambiguous summary rows) and return log2(CPM + 1) per gene.
    """
    opener = gzip.open if file_path.suffix == ".gz" else open
    with opener(file_path, "rt") as handle:
        counts_df = pd.read_csv(handle, sep="\t", header=0)
    gene_col, value_col = counts_df.columns[0], counts_df.columns[1]
    counts_df = counts_df[~counts_df[gene_col].astype(str).str.startswith("N_")]
    values = pd.to_numeric(counts_df[value_col], errors="coerce")
    series = pd.Series(values.to_numpy(), index=counts_df[gene_col].astype(str))
    series = series.groupby(level=0).sum()
    cpm = series / series.sum() * 1.0e6
    return np.log2(cpm + 1.0)


def build_expression_from_supplementary_files(
    series_matrix_path: Path,
    tar_path: Path,
    extract_dir: Path,
    file_loader: Callable[[Path], pd.Series] = load_gene_value_file,
    annotation_path: Path | None = None,
) -> pd.DataFrame:
    """Merge per-sample gene-level supplementary files into one wide expression
    table shaped like the annotated series-matrix table (ID_REF, Gene Symbol,
    Gene Title, one column per GSM accession).

    file_loader controls how each per-sample file is parsed (e.g. plain
    gene/value pairs vs. RSEM genes.results). If annotation_path is given, its
    ID column is matched against the per-sample files' identifiers (e.g.
    Ensembl gene IDs) to look up Gene Symbol/Gene Title; otherwise the raw
    identifier is used as the gene symbol.
    """
    file_map = parse_supplementary_file_map(series_matrix_path)
    extract_supplementary_archive(tar_path, extract_dir)

    sample_series = {}
    for accession, filename in file_map.items():
        file_path = extract_dir / filename
        ensure_exists(file_path, f"Supplementary file for {accession}")
        sample_series[accession] = file_loader(file_path)

    merged_df = pd.DataFrame(sample_series)
    merged_df.index.name = "ID_REF"
    merged_df = merged_df.reset_index()

    if annotation_path is not None:
        annotation_df = load_annotation(
            Path(annotation_path), "ID", "Gene Symbol", "Gene Title"
        )
        merged_df = merge_annotation(
            merged_df, annotation_df, "ID_REF", "ID", "Gene Symbol", "Gene Title"
        )
    else:
        merged_df.insert(1, "Gene Symbol", merged_df["ID_REF"])
        merged_df.insert(2, "Gene Title", "")
    return merged_df


def infer_group_label_from_characteristic(
    characteristic_value: str,
    value_map: dict[str, str],
) -> str:
    lowered = characteristic_value.lower()
    for substring, label in value_map.items():
        if substring.lower() in lowered:
            return label
    return ""


def run_preprocessor(
    mapping_text: str,
    input_file: str | Path,
    output_directory: str | Path,
    annotation_path: str | Path | None = None,
    supplementary_tar_path: str | Path | None = None,
    supplementary_value_format: str = "two_column",
    supplementary_value_column: str = "TPM",
    supplementary_annotation_path: str | Path | None = None,
    mapping_strip_prefix: str | None = None,
    sample_title_prefix: str | None = None,
    phenotype_column_name: str | None = None,
    control_prefix: str | None = None,
    case_prefix: str | None = None,
    control_label: str | None = None,
    case_label: str | None = None,
    phenotype_characteristic_name: str | None = None,
    phenotype_value_map: dict[str, str] | None = None,
) -> tuple[Path, Path]:
    series_matrix_path = Path(input_file)
    output_directory_path = Path(output_directory)

    ensure_exists(series_matrix_path, "Series matrix file")

    mapping = build_column_mapping(mapping_text, mapping_strip_prefix or "")
    original_output_path = output_directory_path / "original_matrix.csv"
    annotated_output_path = output_directory_path / "annotated.csv"
    expression_output_path = output_directory_path / "annotated_renamed_unlogged.csv"
    deduplicated_output_path = output_directory_path / "deduplicated_expression_dataset.csv"
    duplicate_summary_output_path = output_directory_path / "duplicate_probe_summary.csv"
    metadata_output_path = output_directory_path / "sample_metadata.csv"

    output_directory_path.mkdir(parents=True, exist_ok=True)

    if supplementary_tar_path is not None:
        # Series (e.g. RNA-seq) whose values aren't embedded in the series
        # matrix table; per-sample gene-level files are merged directly.
        if supplementary_value_format == "rsem_gene_results":
            file_loader = lambda path: load_rsem_gene_results_file(
                path, value_column=supplementary_value_column
            )
        elif supplementary_value_format == "two_column":
            file_loader = load_gene_value_file
        elif supplementary_value_format == "star_gene_counts":
            file_loader = load_star_gene_counts_file
        else:
            raise ValueError(
                f"Unknown supplementary_value_format: {supplementary_value_format!r}"
            )

        annotated_df = build_expression_from_supplementary_files(
            series_matrix_path,
            Path(supplementary_tar_path),
            output_directory_path / "supplementary_extracted",
            file_loader=file_loader,
            annotation_path=Path(supplementary_annotation_path) if supplementary_annotation_path else None,
        )
        annotated_df.to_csv(original_output_path, index=False)
        annotated_df.to_csv(annotated_output_path, index=False)
    else:
        resolved_annotation_path = resolve_annotation_path(series_matrix_path, annotation_path)
        if resolved_annotation_path is None:
            platform_id = parse_series_platform_id(series_matrix_path)
            if platform_id:
                auto_annotation_path = (
                    output_directory_path / f"{platform_id}_annotation_auto.csv"
                )
                downloaded = download_platform_annotation_from_geo(
                    platform_id,
                    auto_annotation_path,
                )
                if downloaded is not None:
                    print(f"Downloaded platform annotation from GEO: {downloaded}")
                    resolved_annotation_path = downloaded

        if resolved_annotation_path is None:
            platform_id = parse_series_platform_id(series_matrix_path)
            print(
                "WARNING: no matching annotation file was found"
                + (f" for platform {platform_id}" if platform_id else "")
                + "; "
                "falling back to probe-ID-as-symbol annotation. Cell-type marker "
                "matching and gene-level interpretation will not work correctly."
            )
            resolved_annotation_path = generate_fallback_annotation(
                series_matrix_path,
                output_directory_path / "generated_annotation.csv",
            )
        else:
            print(f"Using annotation file: {resolved_annotation_path}")

        expression_df = load_series_matrix_as_dataframe(series_matrix_path)
        expression_df = drop_rows_with_only_reference(expression_df, reference_column="ID_REF")
        expression_df.to_csv(original_output_path, index=False)

        annotation_df = load_annotation(
            resolved_annotation_path,
            "ID",
            "Gene Symbol",
            "Gene Title",
        )
        annotated_df = merge_annotation(
            expression_df,
            annotation_df,
            "ID_REF",
            "ID",
            "Gene Symbol",
            "Gene Title",
        )

        annotated_df.to_csv(annotated_output_path, index=False)

    processed_expression_df, mapped_columns = rename_and_unlog(
        pd.read_csv(annotated_output_path, dtype=str, low_memory=False),
        mapping,
    )

    sample_rows = read_sample_header_rows(series_matrix_path)
    metadata_df = build_metadata_frame(
        sample_rows,
        mapping,
        sample_title_prefix,
        phenotype_column_name,
        control_prefix,
        case_prefix,
        control_label,
        case_label,
        phenotype_characteristic_name,
        phenotype_value_map,
    )
    validate_sample_alignment(metadata_df, processed_expression_df, mapped_columns)

    deduplicated_df, duplicate_summary_df = deduplicate_probes_by_mean_signal(
        processed_expression_df,
        "Gene Symbol",
        mapped_columns,
    )

    processed_expression_df.to_csv(expression_output_path, index=False)
    deduplicated_df.to_csv(deduplicated_output_path, index=False)
    duplicate_summary_df.to_csv(duplicate_summary_output_path, index=False)
    metadata_df.to_csv(metadata_output_path, index=False)

    return expression_output_path, metadata_output_path


# def main() -> None:
#     args = parse_args()
#     series_matrix_path = Path(args.series_matrix)
#     annotation_path = Path(args.annotation)
#     expression_output_path = Path(args.expression_output)
#     metadata_output_path = Path(args.metadata_output)
#     mapping = build_column_mapping(
#         resolve_mapping_text(args),
#         args.mapping_strip_prefix,
#     )

#     ensure_exists(series_matrix_path, "Series matrix file")
#     ensure_exists(annotation_path, "Annotation file")

#     expression_df = load_series_matrix_as_dataframe(series_matrix_path)
#     annotation_df = load_annotation(
#         annotation_path,
#         args.annotation_id_column,
#         args.annotation_symbol_column,
#         args.annotation_title_column,
#     )
#     processed_expression_df, mapped_columns = rename_and_unlog(
#         merge_annotation(
#             expression_df,
#             annotation_df,
#             args.expression_id_column,
#             args.annotation_id_column,
#             args.annotation_symbol_column,
#             args.annotation_title_column,
#         ),
#         mapping,
#     )

#     sample_rows = read_sample_header_rows(series_matrix_path)
#     metadata_df = build_metadata_frame(
#         sample_rows,
#         mapping,
#         args.sample_title_prefix,
#         args.control_prefix,
#         args.case_prefix,
#         args.control_label,
#         args.case_label,
#     )
#     validate_sample_alignment(metadata_df, processed_expression_df, mapped_columns)

#     expression_output_path.parent.mkdir(parents=True, exist_ok=True)
#     metadata_output_path.parent.mkdir(parents=True, exist_ok=True)
#     processed_expression_df.to_csv(expression_output_path, index=False)
#     metadata_df.to_csv(metadata_output_path, index=False)

#     print(f"Expression rows:      {len(processed_expression_df)}")
#     print(f"Mapped sample columns:{len(mapped_columns):>9}")
#     print(f"Expression output:    {expression_output_path}")
#     print(f"Metadata rows:        {len(metadata_df)}")
#     print(f"Metadata output:      {metadata_output_path}")


# if __name__ == "__main__":
#     main()