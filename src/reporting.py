from pathlib import Path
import json
import pandas as pd
from jinja2 import Template

def write_outputs(result: dict, output: str):
    out = Path(output); out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(result['samples']).to_csv(out/'sample_summary.csv', index=False)
    pd.DataFrame(result['variants']).to_csv(out/'variants.csv', index=False)
    (out/'result.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    html = Template("""<html><body><h1>Sequence & Construct QC Report</h1>
    <h2>Samples</h2>{{ samples|safe }}<h2>Variants</h2>{{ variants|safe }}
    <p>First-pass computational report; review complex changes manually.</p></body></html>""").render(
        samples=pd.DataFrame(result['samples']).to_html(index=False),
        variants=pd.DataFrame(result['variants']).to_html(index=False))
    (out/'report.html').write_text(html, encoding='utf-8')
