import argparse, yaml
from src.core import run
from src.reporting import write_outputs

p=argparse.ArgumentParser()
p.add_argument('--reference',required=True)
p.add_argument('--samples',nargs='+',required=True)
p.add_argument('--config',default='configs/default.yaml')
p.add_argument('--output',default='outputs')
a=p.parse_args()
config=yaml.safe_load(open(a.config,encoding='utf-8'))
write_outputs(run(a.reference,a.samples,config.get('motifs',[])),a.output)
print(f'Report written to {a.output}')
