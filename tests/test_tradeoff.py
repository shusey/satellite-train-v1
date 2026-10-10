import pytest
import runpy
from pathlib import Path
from sat_string.tradeoff import relative_to_reference, adjacent_cost, dominated_by


def row(factor=1., error=100., rms=30., switches=100, loss=10.):
    return dict(kp_factor=factor, e_max_m=error, rms_m=rms, transitions=switches,
                extra_altitude_loss_m=loss, high_equivalent_total_s=switches*210.)


def test_reference_percentages_and_signs():
    r=relative_to_reference(row(2.,80.,24.,150,15.),row())
    assert r['e_max_m_reduction_pct']==pytest.approx(20)
    assert r['rms_m_reduction_pct']==pytest.approx(20)
    assert r['transitions_increase_pct']==pytest.approx(50)
    assert r['extra_altitude_loss_m_increase_pct']==pytest.approx(50)


def test_adjacent_deterioration_is_not_an_improvement_cost():
    r=adjacent_cost(row(),row(2.,110.,30.,150,15.))
    assert r['e_max_m_improvement_m']==-10
    assert r['e_max_m_extra_transitions_per_m'] is None
    assert r['rms_m_extra_altitude_loss_per_m'] is None


def test_adjacent_cost_reduction_keeps_signed_value():
    r=adjacent_cost(row(),row(2.,80.,25.,80,8.))
    assert r['e_max_m_extra_transitions_per_m']==-1
    assert r['rms_m_extra_altitude_loss_per_m']==pytest.approx(-.4)


def test_pareto_distinguishes_metric_and_requires_strict_improvement():
    a=row(); b=row(2.,80.,35.,100,10.); c=row(3.,70.,25.,200,20.)
    assert dominated_by(a,[a,b,c],'e_max_m')==[2.]
    assert dominated_by(a,[a,b,c],'rms_m')==[]
    assert dominated_by(a,[a,row(2.)],'e_max_m')==[]


def test_zero_reference_is_undefined_not_infinite():
    r=relative_to_reference(row(),row(error=0,rms=0,switches=0,loss=0))
    assert all(value is None for value in r.values())


def test_compact_keeps_metrics_and_separates_LOW_wait_from_LH_period(tmp_path):
    fn=runpy.run_path(str(Path(__file__).parents[1]/'scripts/compare_kp_tradeoff.py'))['compact']
    (tmp_path/'local_cycles.csv').write_text('actual_next_wait_s,next_same_satellite_start_interval_s\n300,840\n500,1040\n,\n',encoding='utf-8')
    old=dict(analysis_directory=str(tmp_path),first_response=dict(first_local_cycles=[dict(dummy=1)],first_command_time_s=10.))
    new=fn(old)
    assert new['LOW_completion_to_next_HIGH_s']==dict(n=2,min=300.,median=400.,max=500.)
    assert 'first_local_cycles' not in new['first_response']
    assert new['first_response']['first_command_time_s']==10.
    assert 'first_local_cycles' in old['first_response']
