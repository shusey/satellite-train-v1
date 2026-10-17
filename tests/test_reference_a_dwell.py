from pathlib import Path
import runpy
import pytest
import numpy as np
from sat_string.config import load_config
from sat_string.dwell_intervention_metrics import distribution,followup_statistics
from sat_string.dwell_horizon import joined_saved_parts


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
    for days in [30,60]:
        longer=fn(base,days=days)
        check=result.to_dict();check['simulation']['duration_s']=days*86400.
        assert longer.to_dict()==check
    with pytest.raises(ValueError):fn(base,days=14)


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


def test_checkpoint_join_drops_duplicate_and_preserves_state(tmp_path):
    parts=[tmp_path/'first',tmp_path/'next']
    for folder,times in zip(parts,[[0.,10.,20.],[20.,30.,40.]]):
        folder.mkdir();np.save(folder/'time_s.npy',times)
        np.save(folder/'a_m.npy',np.asarray(times)[:,None,None]+100)
    t,x,checks=joined_saved_parts(parts,0,['a_m'])
    np.testing.assert_array_equal(t,[0,10,20,30,40])
    np.testing.assert_array_equal(x['a_m'][:,0],t+100)
    assert checks[0]['max_abs_difference']==0


def test_checkpoint_join_rejects_state_change(tmp_path):
    parts=[tmp_path/'first',tmp_path/'next']
    for folder,times,value in zip(parts,[[0.,10.],[10.,20.]],[1.,2.]):
        folder.mkdir();np.save(folder/'time_s.npy',times)
        np.save(folder/'a_m.npy',np.full((2,1,1),value))
    with pytest.raises(ValueError,match='Checkpoint state changed'):
        joined_saved_parts(parts,0,['a_m'])
