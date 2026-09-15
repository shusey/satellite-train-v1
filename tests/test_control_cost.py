import json
from dataclasses import replace
import numpy as np
import pytest
from sat_string.attitude import SwitchEvent
from sat_string.control_cost import event_activity, evaluate_cost, propagation_metrics, linear_modes, clean_json
from sat_string.controller import switching_function
from sat_string.simulator import simulate_case, CaseDefinition


def event(start,complete,source="L",target="H",sat=1):
    return SwitchEvent(sat,start,start,complete,source,target)


def test_exact_slew_integral_and_high_dwell():
    records=[event(10,130),event(430,550,"H","L")]
    m=event_activity(records,2,600)
    assert m["switch_count_by_satellite"]==[2,0]
    assert m["high_stable_s_by_satellite"]==[300,0]
    assert m["high_equivalent_s_by_satellite"]==[420,0]


def test_partial_slews_delayed_unstarted_and_end_boundary():
    records=[event(10,130),event(430,550,"H","L"),event(900,1020)]
    assert event_activity(records,1,70)["high_equivalent_total_s"]==pytest.approx(11.25)
    assert event_activity(records,1,490)["high_equivalent_total_s"]==pytest.approx(408.75)
    m=event_activity(records,1,430)
    assert m["switch_count"]==2
    assert m["high_stable_total_s"]==300
    assert event_activity(records,1,600)["switch_count"]==2


def test_invalid_event_sequence_rejected():
    with pytest.raises(ValueError):
        event_activity([event(10,130),event(20,140,"H","L")],1,500)


def test_missing_nonmonotone_front_not_reported_as_speed(short_config):
    c=replace(short_config,formation=replace(short_config.formation,satellite_count=21),
              simulation=replace(short_config.simulation,duration_s=1000))
    records=[event(100+20*i,110+20*i,sat=11+i) for i in range(1,11) if i!=3]
    m=propagation_metrics(records,c)
    assert m["range_hops"]==10
    assert m["right"]["first_unreached_satellite_id"]==14
    assert m["right"]["speed_m_s"] is None
    assert "gap_in_front" in m["right"]["fit_status"]
    json.dumps(m,allow_nan=False)


def test_orbital_energy_and_identical_reference(short_config):
    result=simulate_case(short_config,CaseDefinition("low",False,True))
    m=evaluate_cost(result.arrays(),short_config,[],result.arrays())
    assert m["extra_mean_a_loss_vs_low_m"]==0
    assert m["extra_energy_loss_vs_low_J"]==0
    a=result.a_m
    energy_difference=(-short_config.constants.mu_m3_s2/(2*a[-1])
                       +short_config.constants.mu_m3_s2/(2*a[0]))
    assert np.allclose(m["delta_specific_energy_J_kg_by_satellite"],energy_difference,atol=1e-8,rtol=1e-8)
    assert m["mean_a_loss_m"]>0
    assert m["total_energy_loss_J"]>0


def test_laplacian_matches_actual_controller_and_modes(short_config):
    n=short_config.formation.satellite_count
    D=np.zeros((n-1,n))
    for i in range(n-1):
        D[i,i]=-1
        D[i,i+1]=1
    L=D.T@D
    x=np.arange(n,dtype=float)**2
    v=np.cos(x)
    kp=short_config.derived.controller_kp_s2_inv
    kd=short_config.derived.controller_kd_s_inv
    assert np.allclose(switching_function(D@x,D@v,kp,kd),-L@(kp*x+kd*v))
    out=linear_modes(short_config)
    assert np.allclose(np.linalg.eigvalsh(L),[r["laplacian_eigenvalue"] for r in out["modes"]])
    b=out["b_s_inv"]
    full=np.block([[np.zeros((n,n)),np.eye(n)],[-kp*L,b*np.eye(n)-kd*L]])
    expected=np.array([complex(r[f"real_{j}_s_inv"],r[f"imag_{j}_s_inv"]) for r in out["modes"] for j in (1,2)])
    actual=np.linalg.eigvals(full)
    # Compare sets by nearest distance; complex sorting is numerically fragile.
    assert max(min(abs(z-expected)) for z in actual)<1e-10


def test_nonlinear_orbit_acceleration_linearization(short_config):
    c=short_config
    a=c.derived.reference_radius_m
    low=c.spacecraft.area_low_m2
    rho=c.atmosphere.reference_density_kg_m3
    mu=c.constants.mu_m3_s2
    H=c.atmosphere.effective_scale_height_m
    factor=c.spacecraft.drag_coefficient/c.spacecraft.mass_kg
    def acceleration(r,area):
        return 1.5*a*mu/r**2 * rho*np.exp(-(r-a)/H)*factor*area
    da=1.0
    derivative=(acceleration(a+da,low)-acceleration(a-da,low))/(2*da)
    velocity_derivative=-1.5*np.sqrt(mu/a**3)
    model=linear_modes(c)
    assert derivative/velocity_derivative==pytest.approx(model["b_s_inv"],rel=1e-8)
    assert acceleration(a,1.0)==pytest.approx(model["g_m_inv_s2"],rel=1e-12)



def test_analysis_cli_roundtrip_and_output_protection(tmp_path, short_config):
    import subprocess
    import sys
    from pathlib import Path
    import yaml
    root=Path(__file__).parents[1]
    config=tmp_path/"config.yaml"
    config.write_text(yaml.safe_dump(short_config.to_dict()),encoding="utf-8")
    out=tmp_path/"cost"
    command=[sys.executable,str(root/"scripts"/"run_control_cost_analysis.py"),
             "--config",str(config),"--output",str(out)]
    completed=subprocess.run(command+["--stage","baseline"],cwd=root,capture_output=True,text=True,timeout=90)
    assert completed.returncode==0,completed.stderr
    saved=out/"baseline"/"controlled_disturbed"/"results.npz"
    before=saved.stat().st_mtime_ns
    second=subprocess.run(command+["--stage","postprocess"],cwd=root,capture_output=True,text=True,timeout=90)
    assert second.returncode==0,second.stderr
    assert saved.stat().st_mtime_ns==before
    assert (out/"plots"/"01_switching_states.png").is_file()
    assert (out/"baseline_comparison.csv").is_file()
    assert (out/"linear_modes.csv").is_file()
    third=subprocess.run(command+["--stage","baseline"],cwd=root,capture_output=True,text=True,timeout=90)
    assert third.returncode!=0
    assert saved.stat().st_mtime_ns==before

