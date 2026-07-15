from __future__ import annotations

import json
import tempfile
from pathlib import Path
from time import perf_counter

import pandas as pd
import streamlit as st
import yaml

from src.core import run

CONFIG_PATH = Path("configs/default.yaml")
MAX_TOTAL_UPLOAD_MB = 25
MAX_SAMPLES = 25


def load_config() -> dict:
    return yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8")) or {}


@st.cache_data(show_spinner=False, max_entries=20)
def analyze_cached(
    reference_name: str,
    reference_bytes: bytes,
    sample_payload: tuple[tuple[str, bytes], ...],
    motifs_json: str,
) -> dict:
    motifs = json.loads(motifs_json)
    with tempfile.TemporaryDirectory() as directory:
        base = Path(directory)
        reference_path = base / Path(reference_name).name
        reference_path.write_bytes(reference_bytes)

        sample_paths: list[str] = []
        for index, (name, content) in enumerate(sample_payload, start=1):
            safe_name = f"{index:03d}_{Path(name).name}"
            sample_path = base / safe_name
            sample_path.write_bytes(content)
            sample_paths.append(str(sample_path))

        return run(str(reference_path), sample_paths, motifs)


cfg = load_config()
st.set_page_config(page_title=cfg.get("app_title", "Sequence & Construct QC Workbench"), layout="wide")
st.title(cfg.get("app_title", "Sequence & Construct QC Workbench"))
st.caption(cfg.get("app_subtitle", "Annotation-aware comparison of engineered constructs"))
st.info("Use synthetic, public, or authorized non-confidential files only.")

with st.form("construct_qc_form"):
    reference_file = st.file_uploader(
        "Reference GenBank or FASTA",
        type=["gb", "gbk", "genbank", "fa", "fasta"],
    )
    sample_files = st.file_uploader(
        "Sample constructs",
        type=["gb", "gbk", "genbank", "fa", "fasta"],
        accept_multiple_files=True,
    )
    submitted = st.form_submit_button("Run construct QC", type="primary")

if submitted:
    if reference_file is None or not sample_files:
        st.error("Upload one reference and at least one sample construct.")
    elif len(sample_files) > MAX_SAMPLES:
        st.error(f"Upload no more than {MAX_SAMPLES} samples per run.")
    else:
        reference_bytes = reference_file.getvalue()
        sample_payload = tuple((file.name, file.getvalue()) for file in sample_files)
        total_bytes = len(reference_bytes) + sum(len(content) for _, content in sample_payload)

        if total_bytes > MAX_TOTAL_UPLOAD_MB * 1024 * 1024:
            st.error(f"The combined upload exceeds {MAX_TOTAL_UPLOAD_MB} MB.")
        else:
            started = perf_counter()
            try:
                with st.spinner("Comparing constructs and annotating variants..."):
                    result = analyze_cached(
                        reference_file.name,
                        reference_bytes,
                        sample_payload,
                        json.dumps(cfg.get("motifs", []), sort_keys=True),
                    )
                st.session_state["construct_qc_result"] = result
                st.session_state["construct_qc_elapsed"] = perf_counter() - started
            except Exception as exc:
                st.exception(exc)

result = st.session_state.get("construct_qc_result")
if result:
    elapsed = st.session_state.get("construct_qc_elapsed")
    if elapsed is not None:
        st.success(f"Analysis completed in {elapsed:.2f} seconds.")

    sample_df = pd.DataFrame(result.get("samples", []))
    variant_df = pd.DataFrame(result.get("variants", []))

    st.subheader("Sample summary")
    st.dataframe(sample_df, use_container_width=True, hide_index=True)
    st.subheader("Variants")
    st.dataframe(variant_df, use_container_width=True, hide_index=True)

    st.download_button(
        "Download variants CSV",
        variant_df.to_csv(index=False).encode("utf-8"),
        "variants.csv",
        "text/csv",
    )
