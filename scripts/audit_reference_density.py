"""Use existing direct NRLMSIS spot audit without repeating any trajectory."""
import yaml
import audit_dispersion_control_density as audit
from run_reference_dispersion import ROOT,SPEC,DATA
from run_dispersion_control import configuration
spec=yaml.safe_load(SPEC.read_text())
audit.OUT=ROOT/spec['output'];audit.DATA=DATA
audit.configuration=lambda days:(spec,configuration(days)[1])
if __name__=='__main__':raise SystemExit(audit.main())
