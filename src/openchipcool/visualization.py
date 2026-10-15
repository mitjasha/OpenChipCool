"""ParaView-compatible legacy VTK export + transparent steady FEM field plots."""
from pathlib import Path
import numpy as np


def export_legacy_vtk(fem, temperature, path):
    """Pure Python ASCII VTK unstructured tetra mesh with nodal temperature."""
    path=Path(path)
    with path.open("w",encoding="utf-8") as f:
        f.write("# vtk DataFile Version 3.0\nOpenChipCool FEM temperature in Celsius\nASCII\nDATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(fem.coords)} float\n")
        for xyz in fem.coords:
            f.write(f"{xyz[0]:.9g} {xyz[1]:.9g} {xyz[2]:.9g}\n")
        f.write(f"CELLS {len(fem.tets)} {5*len(fem.tets)}\n")
        for tet in fem.tets:
            f.write("4 " + " ".join(str(x) for x in tet) + "\n")
        f.write(f"CELL_TYPES {len(fem.tets)}\n")
        f.write("10\n"*len(fem.tets))
        f.write(f"POINT_DATA {len(temperature)}\nSCALARS Temperature_C float 1\nLOOKUP_TABLE default\n")
        for v in temperature:
            f.write(f"{v:.8g}\n")


def plot_steady_field(fem, temperature, path, spec):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    x=np.unique(fem.coords[:,0]);y=np.unique(fem.coords[:,1]);z=np.unique(fem.coords[:,2])
    if len(x)*len(y)*len(z)!=len(temperature):
        raise ValueError("Expected Cartesian product vertex lattice")
    T=temperature.reshape((len(x),len(y),len(z)))
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout="constrained")
    top=axes[0].pcolormesh(1000*x,1000*y,T[:,:,-1].T,cmap="inferno",shading="auto")
    for xc in (-spec.chip_pitch/2,spec.chip_pitch/2):
        for yc in (-spec.chip_pitch/2,spec.chip_pitch/2):
            axes[0].add_patch(Rectangle((1000*(xc-spec.chip_side/2),1000*(yc-spec.chip_side/2)),
                    1000*spec.chip_side,1000*spec.chip_side,fill=False,edgecolor="cyan",linewidth=1.1,linestyle="--"))
    axes[0].set(title="3D FEM top-face temperature + die outlines",xlabel="X [mm]",ylabel="Y [mm]",aspect="equal")
    fig.colorbar(top,ax=axes[0],label="Temperature [°C]")
    mid=int(np.argmin(abs(y)))
    sl=axes[1].pcolormesh(1000*x,1000*z,T[:,mid,:].T,cmap="inferno",shading="auto")
    axes[1].axhline(spec.z_tim*1000,color="cyan",linestyle="--",linewidth=1,label="Copper / TIM")
    axes[1].axhline(spec.z_die*1000,color="white",linestyle=":",linewidth=1,label="TIM / die layer")
    axes[1].set(title=f"Cross-section y={y[mid]*1000:.1f} mm",xlabel="X [mm]",ylabel="Z [mm]")
    axes[1].legend(fontsize=7,loc="upper left")
    fig.colorbar(sl,ax=axes[1],label="Temperature [°C]")
    fig.savefig(path,dpi=170)
    plt.close(fig)
