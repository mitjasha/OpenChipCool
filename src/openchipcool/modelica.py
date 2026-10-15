"""Export physical Python RC coefficients to a compilable Modelica ODE model.

Not a Modelica.Fluid distributed network. Modelica includes lumped reservoir,
pump-dependent coldplate UA and a radiator with fixed thermal conductance.
"""
from pathlib import Path
from .rc import RCParameters
from .config import Spec
from .literature import ChannelDesign, reference_nusselt, reynolds


def render_modelica(rc:RCParameters,s=Spec()):
    plate=ChannelDesign(plate_width_m=s.side,plate_length_m=s.side,length_m=s.side)
    nb,ns=reference_nusselt()
    # h(flow) is an idealized, paper-anchored Nu integration (no own CFD).
    bottom_weight=plate.bottom_area/plate.planform_area
    side_weight=plate.side_area/plate.planform_area
    re_nominal=reynolds(s.flow_ref_lpm,plate)
    return f'''within ;
model OpenChipCoolLoop
  "Four die nodes, four plate nodes, and one well-mixed coolant loop"
  parameter Real Cdie(unit="J/K") = {rc.C_die:.12g};
  parameter Real Cplate(unit="J/K") = {rc.C_plate:.12g};
  parameter Real Ccool(unit="J/K") = {s.coolant_capacity:.12g};
  parameter Real Gvertical(unit="W/K") = {rc.G_vertical:.12g};
  parameter Real Glateral(unit="W/K") = {rc.G_lateral:.12g};
  parameter Real GconvRef(unit="W/K") = {rc.G_conv_ref:.12g};
  parameter Real Grad(unit="W/K") = {s.radiator_UA:.12g};
  parameter Real ambient(unit="K") = {s.ambient+273.15:.12g};
  parameter Real Tinitial(unit="K") = {s.coolant_initial+273.15:.12g};
  parameter Real nominalPower(unit="W") = 200;
  parameter Real nominalFlow(unit="L/min") = {s.flow_ref_lpm:.12g};
  parameter Real workloadStepTime(unit="s") = 5;
  parameter Real pumpFailureTime(unit="s") = 25;
  parameter Real failedFlowFraction = 0.05;
  parameter Real hotspotMultiplier = 1.75;
  parameter Real ReNominal = {re_nominal:.12g};
  parameter Real ReAnchor = {plate.reference_re:.12g};
  parameter Real ReMinimum = {plate.min_measured_re:.12g};
  parameter Real NuBottom = {nb:.12g};
  parameter Real NuSide = {ns:.12g};
  parameter Real NuBottomLimit = 3.9;
  parameter Real NuSideLimit = 6.8;
  parameter Real bottomWeight = {bottom_weight:.12g};
  parameter Real sideWeight = {side_weight:.12g};
  parameter Real adjacent[4,4] = [0,1,1,0; 1,0,0,1; 1,0,0,1; 0,1,1,0];
  Real Tdie[4](each unit="K", each start=Tinitial, each fixed=true);
  Real Tplate[4](each unit="K", each start=Tinitial, each fixed=true);
  Real Tcool(unit="K", start=Tinitial, fixed=true);
  Real power[4](each unit="W");
  Real flowFraction;
  Real Gconv(unit="W/K");
  Real ReCurrent;
  Real nusseltCurrent;
  parameter Real nusseltNominal = bottomWeight*(NuBottomLimit+(NuBottom-NuBottomLimit)*(ReNominal/ReAnchor)^(1.0/3.0)) + sideWeight*(NuSideLimit+(NuSide-NuSideLimit)*(ReNominal/ReAnchor)^(1.0/3.0));
  Real heatToCoolant(unit="W");
equation
  flowFraction = if time < pumpFailureTime then 1 else failedFlowFraction;
  ReCurrent = ReNominal * max(0, flowFraction);
  nusseltCurrent = bottomWeight*(NuBottomLimit+(NuBottom-NuBottomLimit)*(max(0.000001, ReCurrent)/ReAnchor)^(1.0/3.0)) + sideWeight*(NuSideLimit+(NuSide-NuSideLimit)*(max(0.000001, ReCurrent)/ReAnchor)^(1.0/3.0));
  // NOTE: Re < 139 scaling is an explicitly unvalidated pump-stop extrapolation.
  Gconv = GconvRef * (nusseltCurrent/nusseltNominal) * min(1, ReCurrent/ReMinimum);
  for i in 1:4 loop
    power[i] = nominalPower * (if time < workloadStepTime then 0.5 else
      if (time >= 15 and time < 35 and i == 1) then hotspotMultiplier else 1);
    Cdie * der(Tdie[i]) = power[i] - Gvertical*(Tdie[i]-Tplate[i]);
    Cplate * der(Tplate[i]) = Gvertical*(Tdie[i]-Tplate[i])
      + Glateral * sum(adjacent[i,j]*(Tplate[j]-Tplate[i]) for j in 1:4)
      - Gconv*(Tplate[i]-Tcool);
  end for;
  heatToCoolant = Gconv * sum(Tplate[i]-Tcool for i in 1:4);
  Ccool * der(Tcool) = heatToCoolant - Grad*(Tcool-ambient);
  annotation(experiment(StartTime=0,StopTime=60,Interval=0.1,Tolerance=1e-6));
end OpenChipCoolLoop;
'''


def write_modelica(path,rc,s=Spec()):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(render_modelica(rc,s),encoding="utf-8")
