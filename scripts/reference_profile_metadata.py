"""Publish exact offsets and quantify the old-versus-reference normalization."""
import json
import numpy as np
from run_reference_dispersion import ROOT,DATA
from run_dispersion_control import configuration
from analyze_dispersion_control import csv_rows
from sat_string.reference_condition import load_frozen
from sat_string.passive_dispersion import initial_profile_state
from sat_string.controller import gap_kinematics
from sat_string.io import write_json


def main():
    obj=load_frozen(DATA/'profiles.json');ps=obj['profiles'];_,cfg=configuration(30)
    old=json.loads((ROOT/'docs/dispersion_control_transfer/data/profiles.json').read_text())['profiles'][:3]
    rows=[]
    for p in ps[:3]:
        state,_=initial_profile_state(cfg,p['delta_a_m'],p['B_multiplier'])
        e,v,_=gap_kinematics(*state,reference_radius_m=cfg.derived.reference_radius_m,
            target_spacing_rad=cfg.derived.target_spacing_rad,mu_m3_s2=cfg.constants.mu_m3_s2)
        # Actual initialized radii also expose floating-point offset rounding.
        realized=state[0]-cfg.derived.reference_radius_m
        n0=np.sqrt(cfg.constants.mu_m3_s2/cfg.derived.reference_radius_m**3)
        theory=-1.5*n0*np.diff(p['delta_a_m'])
        rows.append(dict(case=p['case'],specified_a_sample_std_m=float(np.std(p['delta_a_m'],ddof=1)),
            initialized_a_sample_std_m=float(np.std(realized,ddof=1)),
            a_mean_m=float(np.mean(p['delta_a_m'])),B_mean=float(np.mean(p['epsilon'])),
            B_sample_std=float(np.std(p['epsilon'],ddof=1)),
            initial_max_error_m=float(np.abs(e).max()),initial_max_rate_m_s=float(np.abs(v).max()),
            initial_drift_first_order_max_absolute_difference_m_s=float(np.abs(v-theory).max())))
    csv_rows(DATA/'reference_profiles.csv',[dict(satellite_id=j+1,
        reference_delta_a_m=obj['reference_delta_a_profile'][j],reference_epsilon=obj['reference_B_profile'][j],
        B_LOW_m2_kg=.0055*(1+obj['reference_B_profile'][j]),B_HIGH_m2_kg=.0165*(1+obj['reference_B_profile'][j])) for j in range(21)])
    low=[]
    for p in old:
        a=np.asarray(p['delta_a_m']);b=np.asarray(p['epsilon'])
        low.append(dict(case=p['case'],id=p['id'],seed=p['seed'],
            a_population_std_m=float(a.std()),a_sample_std_m=float(a.std(ddof=1)),
            B_population_std=float(b.std()),B_sample_std=float(b.std(ddof=1)),
            profile_shape_comparable_to_reference=p['case']!='AB',
            nominal_a_m=.1 if p['case'] in ('A','AB') else 0,
            nominal_B_fraction=.003 if p['case'] in ('B','AB') else 0))
    write_json(DATA/'profile_metadata.json',dict(reference=rows,old_low=low,
        A_reference_to_old_profile_scale=.5/obj['source_A_sample_std'],
        B_reference_to_old_profile_scale=.005/obj['source_B_sample_std'],
        B_model='fixed multiplier of Cd*A/m through LOW/HIGH/slew; not a detailed attitude-error model',
        initial_position='fixed-reference 100km, unchanged mean longitudes; error zero to floating-point precision',
        ddof=1,statement=obj['statement']))


if __name__=='__main__':main()
