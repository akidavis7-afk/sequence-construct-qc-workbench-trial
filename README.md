# Sequence & Construct QC Workbench

A reusable CLI + Streamlit workbench for plasmids, viral constructs, CRISPR vectors and transgenes.

[Click here to view the Live Interactive Web App Demo](https://sequence-construct-qc-workbench-trial-7g4xmbmzfieikkmxhbtovf.streamlit.app/)

## Capabilities

- Read FASTA and annotated GenBank references/samples
- Pairwise global alignment
- Detect substitutions, insertions and deletions
- Map variants to GenBank features
- First-pass CDS consequence labels
- IUPAC motif gain/loss scanning
- Batch CSV, JSON and self-contained HTML reports

## Repository structure

```text
app.py
cli.py
src/
configs/
examples/
tests/
Dockerfile
README.md
```

## Local setup (Windows CMD)

```cmd
py -m venv .venv
.venv\Scripts\activate
python -m pip install --upgrade pip
pip install -r requirements.txt
pytest -q
python cli.py --reference examples/reference.gb --samples examples/clone_a.fasta examples/clone_b.fasta --config configs/default.yaml --output outputs
streamlit run app.py
```

Open `http://localhost:8501`.

## Docker

```cmd
docker build -t sequence-construct-qc .
docker run --rm -p 8501:8501 sequence-construct-qc
```

## Personalize for a lab

Change only:

- `configs/lab_profile.yaml`
- files under `examples/<lab_name>/`
- app title and introduction values in YAML

Do not overwrite the default examples. Put personalized files in a lab-specific subfolder and point the YAML configuration to that folder.

## Scientific limitations

This project provides first-pass QC. Complex rearrangements, biological function and clinical significance require manual review and specialist tools. Example files are synthetic.

## Streamlit performance

The app uses a submit form, `st.cache_data`, immutable upload bytes, and session-state result persistence. This prevents the analysis from running again when a widget changes. The sequence-comparison apps also use a direct linear comparison for closely related equal-length sequences and reserve global alignment for likely indels or larger differences.
