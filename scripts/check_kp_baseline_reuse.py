"""Verify reused baseline evaluation and every signed endpoint stay identical."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.io import prepare_output_directory,write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--old',type=Path,required=True);p.add_argument('--new',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    checks={}
    for day in [30,60]:
        old=a.old/f'kd_1_{day}d';new=a.new/f'kp_1_{day}d'
        before=json.loads((old/'residual_metrics.json').read_text(encoding='utf-8'));after=json.loads((new/'residual_metrics.json').read_text(encoding='utf-8'))
        rows={key:before[key]==after[key] for key in ['local_cycles','group_cycles','quiet_cycles']}
        summary={key:value==after['summary'][key] for key,value in before['summary'].items()}
        with np.load(old/'cycle_end_states.npz') as x,np.load(new/'cycle_end_states.npz') as y:
            arrays={key:bool(np.array_equal(x[key],y[key])) for key in x.files}
        checks[str(day)]=dict(rows=rows,existing_summary_fields=summary,endpoint_arrays=arrays,
            passed=all(rows.values()) and all(summary.values()) and all(arrays.values()))
    passed=all(x['passed'] for x in checks.values())
    out=prepare_output_directory(a.output);write_json(out/'baseline_reuse.json',dict(passed=passed,cases=checks))
    print(dict(passed=passed),flush=True);return 0 if passed else 1


if __name__=='__main__':raise SystemExit(main())
