from dataclasses import replace
import numpy as np
import pytest
from sat_string.attitude import H,L,L_TO_H
from sat_string.continuation import continue_case,restore_attitude
from sat_string.simulator import STANDARD_CASES,simulate_case


@pytest.mark.parametrize("case",STANDARD_CASES)
def test_continuation_matches_uninterrupted_every_array_and_event(short_config,case):
    first=replace(short_config,simulation=replace(short_config.simulation,duration_s=70.))
    final=replace(short_config,simulation=replace(short_config.simulation,duration_s=140.))
    prefix=simulate_case(first,case)
    saved={key:value.copy() for key,value in prefix.arrays().items()}
    continued=continue_case(prefix.arrays(),first,final,case,prefix.switch_events)
    full=simulate_case(final,case)
    for key,value in full.arrays().items():
        np.testing.assert_allclose(continued.arrays()[key],value,atol=1e-12,rtol=0,err_msg=key)
    assert continued.switch_events==full.switch_events
    assert all(np.array_equal(saved[key],value) for key,value in prefix.arrays().items())


def test_restore_active_pending_and_dwell(short_config):
    events=[dict(satellite_id=1,command_time=1.,start_time=2.3,complete_time=8.8,from_mode="L",to_mode="H")]
    pending=restore_attitude(short_config,events,2.)
    assert pending.mode[0]==L and pending.pending_target[0]==H
    active=restore_attitude(short_config,events,5.)
    assert active.mode[0]==L_TO_H
    complete=restore_attitude(short_config,events,10.)
    assert complete.mode[0]==H
    assert complete.next_allowed_time[0]==pytest.approx(18.3)


def test_continuation_rejects_parameter_changes(short_config):
    first=replace(short_config,simulation=replace(short_config.simulation,duration_s=70.))
    result=simulate_case(first,STANDARD_CASES[0])
    changed=replace(first,controller=replace(first.controller,dwell_time_s=10.),simulation=replace(first.simulation,duration_s=140.))
    with pytest.raises(ValueError,match="only"):
        continue_case(result.arrays(),first,changed,STANDARD_CASES[0],result.switch_events)
