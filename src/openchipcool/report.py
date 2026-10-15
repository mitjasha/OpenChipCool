"""Optional offline-interactive HTML plots; no Streamlit server required."""
from pathlib import Path
import numpy as np
from .config import Spec
from .rc import simulate_loop
from .scenarios import dashboard_power


def write_interactive_html(path,rc,s=Spec()):
    try:
        import plotly.graph_objects as go
    except ImportError:
        return False  # optional package; benchmark itself needs NumPy/SciPy only
    t=np.linspace(0,60,301)
    flows=np.full(len(t),s.flow_ref_lpm)
    cases=[
        ("Base load step",dashboard_power(t,800,"step"),flows),
        ("One-chip power pulse",dashboard_power(t,800,"hotspot"),flows),
        ("Global workload surge",dashboard_power(t,800,"surge"),flows),
        ("Pump flow collapse (25 s)",dashboard_power(t,800,"step"),np.where(t<25,flows,flows*.05)),
    ]
    fig=go.Figure()
    for k,(label,power,flow) in enumerate(cases):
        result=simulate_loop(t,power,flow,rc,s)
        for j in range(4):
            fig.add_trace(go.Scatter(x=t,y=result["die"][:,j],mode="lines",name=f"Die {j+1} mean",
                                     hovertemplate="%{x:.1f} s: %{y:.2f} °C<extra></extra>",visible=k==0))
        fig.add_trace(go.Scatter(x=t,y=result["coolant"],name="Coolant reservoir",line=dict(dash="dot"),visible=k==0))
    buttons=[]
    for k,(label,_,__) in enumerate(cases):
        buttons.append(dict(label=label,method="update",args=[{"visible":[i//5==k for i in range(5*len(cases))]},
                      {"title":f"OpenChipCool RC + 1D cooling loop: {label}"}]))
    fig.add_hline(y=s.junction_limit,line_dash="dash",annotation_text="Illustrative 95 °C die-mean line, not a hotspot prediction")
    fig.update_layout(title="OpenChipCool RC + 1D cooling loop: Base load step",
                      xaxis_title="Time [s]",yaxis_title="Temperature [°C]",
                      template="plotly_white",height=630,
                      updatemenus=[dict(type="dropdown",x=.02,y=1.14,buttons=buttons,active=0)],
                      margin=dict(t=145,b=90),legend=dict(orientation="h",y=-.2))
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    fig.write_html(str(p),include_plotlyjs=True,full_html=True,auto_open=False)
    return True
