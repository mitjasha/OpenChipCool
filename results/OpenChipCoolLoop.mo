within ;
model OpenChipCoolLoop
  "Four die nodes, four plate nodes, and one well-mixed coolant loop"
  parameter Real Cdie(unit="J/K") = 0.127037888605;
  parameter Real Cplate(unit="J/K") = 1.73218888704;
  parameter Real Ccool(unit="J/K") = 42000;
  parameter Real Gvertical(unit="W/K") = 6.41328380251;
  parameter Real Glateral(unit="W/K") = 2.32928777482;
  parameter Real GconvRef(unit="W/K") = 12.2880164934;
  parameter Real Grad(unit="W/K") = 300;
  parameter Real ambient(unit="K") = 298.15;
  parameter Real Tinitial(unit="K") = 303.15;
  parameter Real nominalPower(unit="W") = 200;
  parameter Real nominalFlow(unit="L/min") = 4;
  parameter Real workloadStepTime(unit="s") = 5;
  parameter Real pumpFailureTime(unit="s") = 25;
  parameter Real failedFlowFraction = 0.05;
  parameter Real hotspotMultiplier = 1.75;
  parameter Real ReNominal = 1631.80628911;
  parameter Real ReAnchor = 890;
  parameter Real ReMinimum = 139;
  parameter Real NuBottom = 5.97120535714;
  parameter Real NuSide = 9.1515625;
  parameter Real NuBottomLimit = 3.9;
  parameter Real NuSideLimit = 6.8;
  parameter Real bottomWeight = 0.693;
  parameter Real sideWeight = 4.278;
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
