"""Reproducible OpenChipCool command-line benchmarks."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .config import Spec
from .fem import build_fem,build_aligned_fem,run_fem
from .rc import fit_to_fem,run_rc_fixed_inlet
from .scenarios import make_times,training_power,holdout_power
from .modelica import write_modelica
from .literature import summary as literature_summary
from .visualization import export_legacy_vtk,plot_steady_field
from .report import write_interactive_html


def save_csv(path,headers,columns):
    with open(path,"w",newline="",encoding="utf-8") as fp:
        w=csv.writer(fp); w.writerow(headers); w.writerows(np.column_stack(columns))


def plot_results(path,train,holdout,rc_result,steady):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,2,figsize=(12,4.6))
    t=holdout["times"]
    for i in range(4):
        ax[0].plot(t,holdout["means"][:,i],label=f"Die {i+1} FEM")
        ax[0].plot(t,rc_result[:,i],linestyle="--",alpha=.7,label=f"Die {i+1} RC")
    ax[0].set(title="Unseen transient: FEM vs RC",xlabel="Time [s]",ylabel="Die top mean [°C]")
    ax[0].legend(fontsize=7,ncol=2)
    ax[0].grid(alpha=.2)
    ax[1].plot(t,np.max(holdout["peaks"],axis=1),label="Peak hotspot (FEM)")
    ax[1].plot(t,np.max(holdout["means"],axis=1),label="Highest die mean (FEM)")
    ax[1].plot(t,np.max(rc_result[:,:4],axis=1),linestyle="--",label="Highest die mean (RC)")
    ax[1].axhline(95,color="tab:red",linestyle=":",label="Illustrative 95 °C junction limit")
    ax[1].set(title="Averages are not hotspots",xlabel="Time [s]",ylabel="Temperature [°C]")
    ax[1].legend(fontsize=8);ax[1].grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(path,dpi=170)
    plt.close(fig)


def benchmark(args):
    s=Spec()
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    f=(build_aligned_fem(s,refinement=args.refinement) if args.nx is None else
       build_fem(s,nx=args.nx,n_copper=args.copper_layers,n_tim=args.tim_layers,n_die=args.die_layers))
    powers0=np.full(4,200.)
    field=f.steady(powers0,s.coolant_initial,s.h_ref)
    q_out=f.boundary_heat(field,s.coolant_initial,s.h_ref)
    q_in=float(np.sum(powers0))
    times=make_times(args.duration,args.dt)
    train_p=training_power(times)
    valid_p=holdout_power(times)
    train=run_fem(f,times,train_p,s.coolant_initial,s.h_ref)
    holdout=run_fem(f,times,valid_p,s.coolant_initial,s.h_ref)
    params=fit_to_fem(times,train_p,train["means"],s)
    rc_trace=run_rc_fixed_inlet(times,valid_p,params,s.coolant_initial)
    err=rc_trace[:,:4]-holdout["means"]
    params.validation_rmse_K=float(np.sqrt(np.mean(err**2)))
    summary={
        "project":"OpenChipCool", "status":"Numerical verification only; no experimental validation",
        "mesh":{"style":"interface aligned" if args.nx is None else "legacy unaligned",
                "refinement":args.refinement if args.nx is None else None,
                "nx":len(np.unique(f.coords[:,0]))-1,"ny":len(np.unique(f.coords[:,1]))-1,
                "nodes":len(f.coords),"tetrahedra":len(f.tets),
                "z_copper_layers":2*args.refinement if args.nx is None else args.copper_layers,
                "z_TIM_layers":args.refinement if args.nx is None else args.tim_layers,
                "z_die_layers":2*args.refinement if args.nx is None else args.die_layers},
        "literature_coldplate": literature_summary(s.flow_ref_lpm),
        "assumptions":{"chiplets":4,"rated_power_per_chiplet_W":200,"coolant_inlet_C":s.coolant_initial,
                       "h_ref_W_m2K":s.h_ref,"junction_limit_C":s.junction_limit,
                       "fem_top_averages_not_hotspot_peaks":True},
        "steady":{"thermal_power_W":q_in,"removed_heat_W":q_out,
                  "relative_energy_imbalance":abs(q_out-q_in)/q_in,
                  "die_top_mean_C":f.average_die_temps(field).tolist(),
                  "die_top_peak_C":f.peak_die_temps(field).tolist(),
                  "fails_assumed_junction_limit":bool(np.max(f.peak_die_temps(field))>s.junction_limit)},
        "train":{"power_case":"training_power()","steps":len(times)-1,
                 "rmse_die_mean_K":params.training_rmse_K},
        "holdout":{"power_case":"holdout_power(), disjoint from calibration",
                   "rmse_die_mean_K":params.validation_rmse_K,
                   "max_abs_die_mean_error_K":float(np.max(np.abs(err))),
                   "mean_abs_die_mean_error_K":float(np.mean(np.abs(err)))},
        "important_limitations":["No OWN CFD or experimental validation; h anchored to manually digitized published CFD Nu at one Re", "No fluid pressure-drop / pump curve",
                                 "3D mesh still requires stronger independent verification", "die-top means computed by surface quadrature",
                                 "RC calibrated and tested at the same nominal convection coefficient",
                                 "RC cannot predict local peak hotspot temperature", "No experimental validation"]
    }
    (out/"rc_calibration.json").write_text(json.dumps(params.as_dict(),indent=2),encoding="utf-8")
    (out/"benchmark_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    save_csv(out/"holdout_comparison.csv",["time_s"]+[f"fem_mean_die{i+1}_C" for i in range(4)]+
             [f"fem_peak_die{i+1}_C" for i in range(4)]+[f"rc_mean_die{i+1}_C" for i in range(4)],
             [times]+[holdout["means"][:,i] for i in range(4)]+[holdout["peaks"][:,i] for i in range(4)]+[rc_trace[:,i] for i in range(4)])
    save_csv(out/"steady_temperature_points.csv",["x_m","y_m","z_m","T_C"],
             [f.coords[:,0],f.coords[:,1],f.coords[:,2],field])
    export_legacy_vtk(f,field,out/"steady_temperature.vtk")
    plot_steady_field(f,field,out/"temperature_field.png",s)
    write_modelica(out/"OpenChipCoolLoop.mo",params,s)
    write_interactive_html(out/"interactive_report.html",params,s)
    plot_results(out/"benchmark.png",train,holdout,rc_trace,summary)
    print(json.dumps(summary,indent=2))
    print("Saved:",out.resolve())


def richardson_gci(values,refinements=(1,2,3)):
    """Estimate order and fine-grid GCI for systematic three-grid refinement.

    Values ordered coarse, medium, fine; length scale h~1/refinement.
    If non-monotonic, diverging or noisy, report inconclusive rather than
    inventing a positive uncertainty. Safety factor 1.25 is a convention,
    not a confidence interval. Refs: Celik et al. (2008); Roache (1994).
    """
    from scipy.optimize import brentq
    v=np.asarray(values,dtype=float)
    r=np.asarray(refinements,dtype=float)
    if v.shape!=(3,) or r.shape!=(3,) or np.any(np.diff(r)<=0):
        raise ValueError("Expected three increasingly refined grids")
    d21=v[0]-v[1];d32=v[1]-v[2]
    status={"status":"inconclusive", "observed_order":None,
            "gci_fine_percent":None,"extrapolated_value":None,
            "asymptotic_ratio":None}
    if abs(d21)<1e-12 or abs(d32)<1e-12 or d21*d32<=0:
        status["reason"]="zero or oscillatory differences; cannot use classical GCI"
        return status
    h=1/r
    ratio=d21/d32
    def fn(order):
        return (h[0]**order-h[1]**order)/(h[1]**order-h[2]**order)-ratio
    # Implicit minimum observed-order ratio for this unequal refinement.
    if fn(0.00001)>0 or fn(12)<0:
        status["reason"]="observed convergence outside p=(0,12) or asymptotic region"
        return status
    order=float(brentq(fn,0.00001,12))
    fine_ratio=r[2]/r[1]
    err=(v[2]-v[1])/(fine_ratio**order-1)
    extrap=float(v[2]+err)
    gci=125*abs(err)/max(abs(v[2]),1e-12)
    # ASME-Celik asymptotic check: GCI_medium / (r^p GCI_fine) ~ 1.
    coarse_ratio=r[1]/r[0]
    gci_middle=1.25*abs((v[1]-v[0])/v[1])/(coarse_ratio**order-1)
    gci_fine=1.25*abs((v[2]-v[1])/v[2])/(fine_ratio**order-1)
    asym=float(gci_middle/(fine_ratio**order*gci_fine)) if gci_fine else None
    return {"status":"estimated; asymptotic assumption must be checked",
            "observed_order":order,"gci_fine_percent":float(gci),
            "extrapolated_value":extrap,"asymptotic_ratio":asym}


def mesh_study(args):
    """Interface-aligned refinement of ALL three dimensions with GCI when valid."""
    s=Spec();out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    refinements=[int(v) for v in args.refinements.split(",")]
    if len(refinements)!=3 or any(refinements[i]>=refinements[i+1] for i in range(2)):
        raise ValueError("Specify exactly three increasing factors, e.g., 1,2,3")
    records=[]
    for scale in refinements:
        f=build_aligned_fem(s,refinement=scale)
        u=f.steady(h=s.h_ref,inlet=s.coolant_initial)
        m=f.average_die_temps(u); peak=f.peak_die_temps(u)
        q=f.boundary_heat(u,s.coolant_initial,s.h_ref)
        records.append([scale,len(f.coords),len(f.tets),float(np.mean(m)),
                        float(np.max(peak)),q, float(max(abs(f.source_volume-s.chip_side**2*s.silicon_thickness)))])
        print(f"refinement {scale}: nodes={len(f.coords)}, tets={len(f.tets)}, mean={np.mean(m):.4f}C, max={max(peak):.4f}C, balance={(q-800):.3e}W",flush=True)
    a=np.asarray(records)
    save_csv(out/"mesh_study_aligned.csv",
             ["refinement","nodes","tetrahedra","mean_die_C","max_hotspot_C","Qout_W","max_die_volume_err_m3"],
             [a[:,i] for i in range(7)])
    quantities={"mean_die_C":richardson_gci(a[:,3],refinements),
                "max_hotspot_C":richardson_gci(a[:,4],refinements)}
    result={"method":"nested three-dimensional interface-aligned grids; 1.25 safety factor; not a validation", "refinements":refinements,
            "source":"Celik et al. 2008 doi:10.1115/1.2960953", "h_W_m2K":s.h_ref,
            "observables":quantities, "geometry_invariant":bool(np.max(a[:,6])<1e-17)}
    (out/"mesh_gci.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,4.5))
    ax.plot(1/a[:,0],a[:,3],"o-",label="Die-top area-weighted mean [C]")
    ax.plot(1/a[:,0],a[:,4],"s--",label="Mesh-node peak [C]")
    ax.set(xlabel="Normalized cell size h / h_coarsest",ylabel="Steady temperature [C]",
           title="Interface-aligned 3D refinement: numerical sensitivity")
    ax.invert_xaxis();ax.grid(alpha=.25);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/"mesh_study_aligned.png",dpi=170);plt.close(fig)
    print(json.dumps(result,indent=2))



def main(argv=None):
    p=argparse.ArgumentParser(description="OpenChipCool educational FEM / RC benchmark")
    commands=p.add_subparsers(dest="command",required=True)
    b=commands.add_parser("benchmark",help="Run FEM training, fit RC, validate on holdout, export Modelica")
    m=commands.add_parser("mesh-study",help="Check spatial resolution dependence of steady solutions")
    for sub in (b,m):
        sub.add_argument("--output",default="results")
        sub.add_argument("--copper-layers",type=int,default=3)
        sub.add_argument("--tim-layers",type=int,default=1)
        sub.add_argument("--die-layers",type=int,default=2)
    b.add_argument("--nx",type=int,default=None,help="Optional legacy unaligned nx; do not use for V&V")
    b.add_argument("--refinement",type=int,default=2,help="Aligned mesh scale, default 2")
    b.add_argument("--duration",type=float,default=6.)
    b.add_argument("--dt",type=float,default=.15)
    m.add_argument("--refinements",default="1,2,3",
                   help="Three geometry-aligned 3D refinements (default 1,2,3)")
    args=p.parse_args(argv)
    if args.command=="benchmark":benchmark(args)
    else:mesh_study(args)

if __name__=="__main__":main()
