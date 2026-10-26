import copy
import json
from dataclasses import replace
import numpy as np
import pytest
from sat_string.config import load_config
from sat_string.n_scaling import CentralNodePulse,scaling_configuration,thirty_day_extension
from sat_string.dispersion_control import integrate_control_profiles


def checkpoint(n=21):
    cfg=scaling_configuration(load_config('configs/baseline.yaml'),n,14)
    profiles=[dict(control_enabled=True,delta_a_m=[0.]*n,B_multiplier=[1.]*n)]*2
    return dict(configuration=cfg.to_dict(),time_s=14*86400.,profiles=profiles,disturbance_enabled=False,
        node_input=CentralNodePulse(n).metadata(),state=[np.full(2*n,cfg.derived.reference_radius_m).tolist(),
        np.tile(np.arange(n)*cfg.derived.target_spacing_rad,2).tolist()])


@pytest.mark.parametrize('n',[21,301])
def test_only_duration_changes(n):
    cp=checkpoint(n);original=copy.deepcopy(cp)
    cfg,pulse=thirty_day_extension(cp)
    expected=copy.deepcopy(cp['configuration']);expected['simulation']['duration_s']=30*86400.
    assert cfg.to_dict()==expected and cp==original
    assert np.count_nonzero(pulse(14*86400+1))==0


@pytest.mark.parametrize('bad',['time','input','dispersion','density','order','state','N'])
def test_reject_changed_inputs(bad):
    cp=checkpoint()
    if bad=='time':cp['time_s']-=10
    if bad=='input':cp['node_input']['amplitude_m_s2']*=2
    if bad=='dispersion':cp['profiles'][0]['delta_a_m'][0]=.1
    if bad=='density':cp['disturbance_enabled']=True
    if bad=='order':cp['state'][1][1]=-1
    if bad=='state':cp['state'][0][0]=float('nan')
    if bad=='N':cp=checkpoint(11)
    with pytest.raises(ValueError):thirty_day_extension(cp)


def test_split_during_slew_exactly_matches_uninterrupted(tmp_path):
    cfg=scaling_configuration(load_config('configs/baseline.yaml'),11,14)
    cfg=replace(cfg,simulation=replace(cfg.simulation,duration_s=1000.))
    prefix=replace(cfg,simulation=replace(cfg.simulation,duration_s=100.))
    profiles=[dict(control_enabled=True,delta_a_m=[0.]*11,B_multiplier=[1.]*11)]*2
    pulse=CentralNodePulse(11,1e-3,0,20)
    kw=dict(node_acceleration=pulse,node_input_metadata=pulse.metadata())
    integrate_control_profiles(cfg,profiles,None,tmp_path/'whole',**kw)
    integrate_control_profiles(prefix,profiles,None,tmp_path/'prefix',**kw)
    cp=json.loads((tmp_path/'prefix/checkpoint.json').read_text())
    assert any(mode in [1,3] for mode in cp['modes'][1])
    integrate_control_profiles(cfg,profiles,None,tmp_path/'suffix',checkpoint=cp,**kw)
    for key in ['a_m','lambda_rad','mode','area_m2','density_kg_m3']:
        joined=np.concatenate([np.load(tmp_path/'prefix'/f'{key}.npy'),np.load(tmp_path/'suffix'/f'{key}.npy')[1:]])
        np.testing.assert_array_equal(joined,np.load(tmp_path/'whole'/f'{key}.npy'))
    full=json.loads((tmp_path/'whole/checkpoint.json').read_text())
    suffix=json.loads((tmp_path/'suffix/checkpoint.json').read_text())
    assert full['events']==suffix['events']
