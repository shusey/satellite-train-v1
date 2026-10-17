from pathlib import Path
import runpy
import pytest
from sat_string.config import load_config
from sat_string.dwell_intervention_metrics import distribution,followup_statistics


def test_intervention_changes_only_dwell_and_duration():
    root=Path(__file__).parents[1]
    fn=runpy.run_path(str(root/'scripts/run_reference_a_dwell.py'))['intervention_config']
    base=load_config(root/'tests/fixtures/regression_24h/config.yaml')
    result=fn(base)
    expected=base.to_dict();expected['controller']['dwell_time_s']=150.;expected['simulation']['duration_s']=604800.
    assert result.to_dict()==expected
    assert result.derived.controller_kp_s2_inv==base.derived.controller_kp_s2_inv
    assert result.derived.controller_kd_s_inv==base.derived.controller_kd_s_inv
    for value in [60,100,200,300,600]:
        with pytest.raises(ValueError):fn(base,dwell=value)
    with pytest.raises(ValueError):fn(base,days=60)


def test_residual_distribution_percentiles():
    stats=distribution([0.,1.,2.,3.,4.])
    assert stats['median']==2 and stats['p25']==1 and stats['p75']==3
    assert distribution([])['median'] is None


def test_followup_does_not_turn_censoring_into_no_retrigger():
    rows=[dict(end_s=0,next_wait_s=300),dict(end_s=100,next_wait_s=None),
          dict(end_s=9000,next_wait_s=None)]
    stats=followup_statistics(rows,10000)
    assert stats['right_censored']==2
    assert stats['within_1800s']['known_outcomes']==2
    assert stats['within_1800s']['fraction']==.5
    assert stats['within_1800s']['observed_waits_only_fraction']==1
