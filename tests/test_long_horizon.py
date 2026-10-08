from pathlib import Path
import json
import numpy as np
import pytest
from sat_string.config import load_config
from sat_string.controller import switching_function
from sat_string.io import load_case_arrays,load_switch_events
from sat_string.long_horizon import decompose_q,group_episodes,infer_episode_gap,analyze_long_horizon


def event(i,start,rising=True):
    return dict(satellite_id=i,command_time=start,start_time=start,complete_time=start+2,
                from_mode="L" if rising else "H",to_mode="H" if rising else "L")


def test_q_decomposition_matches_implemented_controller():
    e=np.array([[1.,-2.,3.],[4.,5.,-6.]])
    v=np.array([[.001,.002,-.003],[.004,-.005,.006]])
    qp,qd=decompose_q(e,v,2e-9,9e-5)
    expected=np.array([switching_function(x,y,2e-9,9e-5) for x,y in zip(e,v)])
    np.testing.assert_allclose(qp+qd,expected,rtol=1e-14,atol=1e-22)


def test_global_wave_includes_central_first_cycle_in_second_wave():
    events=[event(1,0),event(1,5,False),event(3,10),event(3,15,False),
            event(1,1000),event(1,1005,False),event(3,1010),event(3,1015,False),
            event(2,1020),event(2,1025,False)]
    gap=infer_episode_gap(events)
    groups=group_episodes(events,1100,gap["gap_s"])
    assert len(groups)==2 and len(groups[0])==4 and len(groups[1])==6
    assert {e["satellite_id"] for e in groups[1]}=={1,2,3}


def test_group_preserves_uncompleted_cycle():
    groups=group_episodes([event(1,10),event(1,20,False)],21,100)
    assert groups[0][-1]["complete_time"]==22


def test_inference_rejects_ambiguous_gaps():
    with pytest.raises(ValueError,match="clear"):
        infer_episode_gap([event(1,t) for t in [0,10,20,30]])


def test_24h_episode_matches_frozen_events_and_cost():
    root=Path(__file__).parent/"fixtures/regression_24h"
    config=load_config(root/"config.yaml")
    arrays=load_case_arrays(root/"controlled_disturbed")
    low=load_case_arrays(root/"all_low_drag_baseline")
    qp,qd=decompose_q(arrays["gap_error_m"],arrays["gap_rate_m_s"],
                       config.derived.controller_kp_s2_inv,config.derived.controller_kd_s_inv)
    arrays["q_m_s2"]=qp+qd
    report,derived=analyze_long_horizon(arrays,config,load_switch_events(root/"controlled_disturbed"),low,86400.)
    wave=report["waves"][0]
    assert wave["transition_starts"]==40 and wave["complete"]
    assert wave["l_to_h_satellite_count"]==20
    assert 11 not in wave["l_to_h_satellite_ids"]
    assert wave["first_satellite_ids"]==[12]
    assert wave["origin_candidates"]==[10,12]
    assert wave["outer_branch_fits"]["left"]["speed_m_s"]==pytest.approx(20.263407114)
    assert report["q_decomposition_max_residual_m_s2"]==0
    json.dumps(report,allow_nan=False)
    assert derived["extra_mean_loss_m"][-1]==pytest.approx(.8211498873)
