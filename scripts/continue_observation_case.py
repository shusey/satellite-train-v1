#!/usr/bin/env python3
"""Extend a saved fixed-input run without replaying its existing prefix."""
import argparse
from pathlib import Path
import sys
from time import perf_counter
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sat_string.config import load_config,config_from_dict
from sat_string.continuation import continue_case
from sat_string.io import load_case_arrays,load_metadata,load_switch_events,prepare_output_directory,save_case,write_json
from sat_string.metrics import compute_case_metrics,save_metrics
from sat_string.simulator import CaseDefinition


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--days",type=float,required=True)
    args=p.parse_args()
    original=load_config(args.source/"normalized_config.yaml")
    data=original.to_dict();data["simulation"]["duration_s"]=args.days*86400
    config=config_from_dict(data)
    case=CaseDefinition(**load_metadata(args.source)["case"])
    arrays=load_case_arrays(args.source);events=load_switch_events(args.source)
    output=prepare_output_directory(args.output)
    start=perf_counter()
    result=continue_case(arrays,original,config,case,events,
                        progress_callback=lambda t,end:print(f"{case.name}: day {t/86400:.2f}/{end/86400:.0f} ({perf_counter()-start:.1f} s)",flush=True))
    save_case(result,config,output)
    save_metrics(compute_case_metrics(result.arrays(),config,result.switch_events),output)
    write_json(output/"continuation.json",dict(source_directory=str(args.source.resolve()),
        resume_time_s=float(arrays["time_s"][-1]),configuration_change="simulation.duration_s only",elapsed_s=perf_counter()-start))
    print(f"Saved {output} ({perf_counter()-start:.1f} s)",flush=True)


if __name__=="__main__":
    main()
