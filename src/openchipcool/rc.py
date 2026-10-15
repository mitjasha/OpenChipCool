"""Eight-state positive-conductance physical RC network plus 1D coolant loop.

Four die junction nodes + four copper quadrant nodes. TIM is represented by
vertical resistances; lateral copper spreading by positive coupling resistors.
Values are first estimated from geometry and may be calibrated on *training*
FEM traces. Independently generated hold-out FEM traces quantify prediction error.
The low-order model approximates DIE AVERAGES, not spatial hotspot maxima.
"""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.optimize import least_squares
from .config import Spec, h_from_flow


@dataclass
class RCParameters:
    C_die: float
    C_plate: float
    G_vertical: float
    G_lateral: float
    G_conv_ref: float
    calibrated: bool = False
    training_rmse_K: float | None = None
    validation_rmse_K: float | None = None

    def as_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, values):
        return cls(**{k: v for k,v in values.items() if k in cls.__dataclass_fields__})


def from_geometry(s=Spec()):
    a = s.chip_side**2
    c_die = s.volumetric_heat_silicon*a*s.silicon_thickness
    c_plate = (s.volumetric_heat_copper*s.side**2*s.copper_thickness +
               s.volumetric_heat_tim*s.side**2*s.tim_thickness) / 4
    r_vertical = (s.silicon_thickness/(2*s.k_silicon*a) +
                  s.tim_thickness/(s.k_tim*a) +
                  s.copper_thickness/(2*s.k_copper*a))
    g_lateral = s.k_copper * s.copper_thickness  # approx square quadrant, A/L=thickness
    g_conv = s.h_ref*s.side**2/4
    return RCParameters(float(c_die),float(c_plate),float(1/r_vertical),float(g_lateral),float(g_conv))


def _laplacian(rc, conv_gain=1.0):
    K = np.zeros((8,8))
    for i in range(4):
        a,b = i,i+4
        g = rc.G_vertical
        K[a,a]+=g; K[b,b]+=g; K[a,b]-=g; K[b,a]-=g
    for i,j in [(0,1),(0,2),(1,3),(2,3)]:
        a,b = i+4,j+4
        g = rc.G_lateral
        K[a,a]+=g; K[b,b]+=g; K[a,b]-=g; K[b,a]-=g
    for i in range(4,8):
        K[i,i] += rc.G_conv_ref * conv_gain
    return K


def run_rc_fixed_inlet(times,powers,rc,inlet=30.0):
    """Same implicit Euler discretization as FEM (fixed h and inlet)."""
    times = np.asarray(times,dtype=float)
    powers = np.asarray(powers,dtype=float)
    dt = times[1]-times[0]
    if not np.allclose(np.diff(times),dt) or powers.shape!=(len(times),4):
        raise ValueError("Inconsistent time grid or power array")
    C = np.array([rc.C_die]*4+[rc.C_plate]*4)
    K = _laplacian(rc)
    LHS = np.diag(C/dt)+K
    x = np.full(8,inlet)
    states = np.empty((len(times),8)); states[0]=x
    b = np.array([0.]*4 + [rc.G_conv_ref*inlet]*4)
    for j in range(1,len(times)):
        x = np.linalg.solve(LHS,C/dt*x + np.r_[powers[j],np.zeros(4)] + b)
        states[j] = x
    return states


def fit_to_fem(times,powers,fem_die_means,s=Spec()):
    """Positive least squares calibration; no information from holdout used.

    Five parameters: vertical / lateral / convective conductances and die / base
    capacitances. NOT a validation; report held-out performance separately.
    """
    base=from_geometry(s)
    # Permit broad but physically positive adjustment of initial estimates.
    lo = np.array([.2,.05,.3,.2,.3])
    hi = np.array([5.,15.,4.,10.,5.])
    def convert(log_factors):
        a=np.exp(log_factors)
        return RCParameters(base.C_die*a[3], base.C_plate*a[4],
                            base.G_vertical*a[0],base.G_lateral*a[1],
                            base.G_conv_ref*a[2], calibrated=True)
    def residual(log_factors):
        r=convert(log_factors)
        x=run_rc_fixed_inlet(times,powers,r,s.coolant_initial)
        return (x[:,0:4]-fem_die_means).ravel()
    fit=least_squares(residual,np.zeros(5),bounds=(np.log(lo),np.log(hi)),max_nfev=80,
                      xtol=1e-7,ftol=1e-7,gtol=1e-7)
    result=convert(fit.x)
    result.training_rmse_K=float(np.sqrt(np.mean(fit.fun**2)))
    return result


def simulate_loop(times,powers,flow_lpm,rc,s=Spec(),initial=None):
    """Nine-state lumped 1D cooling-loop ODE via implicit Euler.

    The loop has one well-mixed coolant reservoir and a constant-UA radiator.
    This is NOT a distributed fluid network; pressure drop and piping are absent.
    h(flow) is anchored to digitized CFD Nusselt profiles from Qu/Mudawar
    (2002) and approximated outside their single Re case. Below Re=139 a
    deliberate, UNVALIDATED shutdown extrapolation is applied.
    """
    times=np.asarray(times,dtype=float)
    powers=np.asarray(powers,dtype=float)
    flows=np.broadcast_to(np.asarray(flow_lpm,dtype=float),times.shape)
    if powers.shape!=(len(times),4) or len(times)<2 or np.any(np.diff(times)<=0):
        raise ValueError("bad grid")
    C=np.array([rc.C_die]*4+[rc.C_plate]*4+[s.coolant_capacity])
    states=np.empty((len(times),9))
    states[0] = np.full(9,s.coolant_initial) if initial is None else np.asarray(initial,dtype=float)
    for j in range(1,len(times)):
        dt=times[j]-times[j-1]
        g=rc.G_conv_ref*h_from_flow(float(flows[j]),s)/s.h_ref
        K=np.zeros((9,9))
        K[:8,:8] = _laplacian(rc,conv_gain=0.)
        for i in range(4):
            p=i+4; q=8
            K[p,p]+=g; K[q,q]+=g; K[p,q]-=g; K[q,p]-=g
        K[8,8] += s.radiator_UA
        rhs=np.r_[powers[j],np.zeros(4),s.radiator_UA*s.ambient]
        states[j]=np.linalg.solve(K+np.diag(C/dt),rhs+(C/dt)*states[j-1])
    heat_to_coolant = np.array([
        4*rc.G_conv_ref*h_from_flow(f,s)/s.h_ref*(np.mean(state[4:8])-state[8])
        for f,state in zip(flows,states)])
    return {"times":times,"states":states,"die":states[:,:4],"plate":states[:,4:8],
            "coolant":states[:,8],"flows":flows,"heat_to_coolant":heat_to_coolant}
