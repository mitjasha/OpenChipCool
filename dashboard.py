"""Streamlit interface for the calibrated reduced model (NOT an FEM/CFD solver)."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from openchipcool.config import Spec
from openchipcool.rc import from_geometry,RCParameters,simulate_loop
from openchipcool.scenarios import dashboard_power

ROOT=Path(__file__).resolve().parent
s=Spec()
st.set_page_config(page_title="OpenChipCool | Thermal Digital Twin",layout="wide")
st.title("OpenChipCool — Thermal Digital Twin")
st.caption("3D FEM-derived low-order thermal model + simplified 1D lumped coolant loop. Synthetic package geometry, Qu/Mudawar (2002) CFD-Nu-anchored equivalent HTC; NOT experimentally validated.")
cal=ROOT/"results"/"rc_calibration.json"
if cal.exists():
    rc=RCParameters.from_dict(json.loads(cal.read_text()))
    st.success("FEM-calibrated RC coefficients loaded. Model validated only against held-out synthetic FEM traces at nominal flow.")
else:
    rc=from_geometry(s)
    st.warning("Geometry-estimated, uncalibrated RC parameters. Run `openchipcool benchmark` to calibrate against FEM.")
left,right=st.columns([1,2])
with left:
    tdp=st.slider("Nominal total TDP [W]",200,1200,800,50)
    nominal_flow=st.slider("Nominal coolant flow [L/min]",1.,10.,4.,.5)
    event=st.selectbox("Load profile",["step","hotspot","surge"],format_func=lambda v:{"step":"Base load step","hotspot":"One-chip hotspot event","surge":"Total power surge"}[v])
    failure=st.checkbox("Pump flow reduction at 25 s",value=True)
    reduction=st.slider("Remaining flow after fault [%]",0,100,5,5)
    duration=st.slider("Simulated duration [s]",30,180,60,10)
times=np.linspace(0,duration,int(duration/.2)+1)
powers=dashboard_power(times,tdp,event)
flows=np.full(len(times),nominal_flow)
if failure: flows[times>=25]*=reduction/100
res=simulate_loop(times,powers,flows,rc,s)
with right:
    stat1,stat2,stat3=st.columns(3)
    stat1.metric("Max die mean",f"{res['die'].max():.1f} °C")
    stat2.metric("Coolant at end",f"{res['coolant'][-1]:.1f} °C")
    stat3.metric("Die mean threshold",">95 °C" if res['die'].max()>95 else "≤95 °C")
    fig=go.Figure()
    for i in range(4):
        fig.add_trace(go.Scatter(x=times,y=res['die'][:,i],name=f"Die {i+1} avg"))
    fig.add_trace(go.Scatter(x=times,y=res['coolant'],name="Well-mixed coolant",line=dict(dash="dot")))
    fig.add_hline(y=s.junction_limit,line_dash="dash",annotation_text="95 °C (illustrative)")
    if failure:fig.add_vline(x=25,line_dash="dot",annotation_text="Pump fault")
    fig.update_layout(title="Transient average die temperatures (RC only)",xaxis_title="Time [s]",yaxis_title="Temperature [°C]",height=410)
    st.plotly_chart(fig,use_container_width=True)
    fig2=go.Figure()
    fig2.add_trace(go.Scatter(x=times,y=np.sum(powers,axis=1),name="Die power",yaxis="y"))
    fig2.add_trace(go.Scatter(x=times,y=flows,name="Flow",yaxis="y2"))
    fig2.update_layout(title="Workload and pump fault inputs",xaxis_title="Time [s]",yaxis=dict(title="Power [W]"),yaxis2=dict(title="Flow [L/min]",overlaying="y",side="right"),height=280)
    st.plotly_chart(fig2,use_container_width=True)
    df=pd.DataFrame({"time_s":times,"total_power_W":np.sum(powers,axis=1),"flow_L_min":flows,
                     "coolant_C":res['coolant'],**{f"die_{i+1}_mean_C":res['die'][:,i] for i in range(4)}})
    st.download_button("Download RC trace (.csv)",df.to_csv(index=False),"openchipcool_rc_trace.csv","text/csv")
st.info("The displayed die states approximate FEM top-face AVERAGES, not local hotspot maxima. This dashboard does not run 3D FEM or OpenModelica. Near-stall HTC is unvalidated extrapolation from CFD-Nu-based operating conditions; no pressure-drop or two-phase physics is included. The 95 °C limit is a synthetic design assumption.")
bench=ROOT/"results"/"benchmark_summary.json"
if bench.exists():
    summary=json.loads(bench.read_text())
    st.caption(f"Offline RC-vs-FEM hold-out RMSE (nominal flow only): {summary['holdout']['rmse_die_mean_K']:.3f} K | FEM steady maximum hotspot: {max(summary['steady']['die_top_peak_C']):.1f} °C")
    if (ROOT/"results"/"benchmark.png").exists():
        st.image(str(ROOT/"results"/"benchmark.png"),caption="Offline FEM/RC benchmark; distinct from above pump-fault scenario")
