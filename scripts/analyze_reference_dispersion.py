"""Reuse transfer/residual analysis for six frozen reference batch members."""
import yaml
import analyze_dispersion_control as analysis
from run_reference_dispersion import ROOT, SPEC, DOC, DATA
from run_dispersion_control import configuration as baseline_configuration
from sat_string.reference_condition import load_frozen
from sat_string.reference_metrics import evaluate_reference

spec=yaml.safe_load(SPEC.read_text())
# Legacy analyzer uses this path only to label reuse provenance. All reference
# OFF trajectories are new because their exact normalized profiles are new.
spec['screening_output']='outputs/initial_dispersion_20261011'
analysis.OUT=ROOT/spec['output']; analysis.DATA=DATA
analysis.configuration=lambda days:(spec,baseline_configuration(days)[1])
analysis.frozen_profiles=lambda _:load_frozen(ROOT/spec['frozen_profiles'])['profiles']


def sources(days,case,p,spec,on=True):
    paths=[analysis.OUT/'batch_30d']+([analysis.OUT/'batch_60d'] if days==60 else [])
    return paths,p['id'],case if on else case+'_OFF'


analysis.sources=sources
analysis.evaluate_transfer=evaluate_reference

if __name__=='__main__': analysis.main()
