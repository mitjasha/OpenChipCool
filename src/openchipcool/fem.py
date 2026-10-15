"""Actual 3D linear-tetrahedral Galerkin heat conduction FEM (SciPy sparse).

This is a deliberately small verification-oriented FEM implementation, not
FEniCS/Elmer, CFD or a substitute for a qualified industrial solver.
Boundary conditions: insulated sides/top; lumped Robin convection at z=0.
Mixed top silicon/epoxy and uniform TIM/copper base are simplifying assumptions.
"""
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import factorized
from .config import Spec, chip_centers


@dataclass
class FEM:
    coords: np.ndarray
    tets: np.ndarray
    K: sparse.csr_matrix
    C: np.ndarray  # lumped mass, J/K per vertex
    area: np.ndarray  # convective surface nodal areas in m^2
    Q: np.ndarray  # nodal load per 1 W per chiplet, shape (N,4)
    observation: np.ndarray  # surface-area-weighted node weights, shape (4,N)
    chip_top_nodes: list
    volume_by_material: dict
    source_volume: np.ndarray

    def matrix(self, h):
        return self.K + sparse.diags(h * self.area, format="csr")

    def steady(self, power=(200, 200, 200, 200), inlet=30., h=26000.):
        rhs = self.Q @ np.asarray(power, dtype=float) + h * self.area * inlet
        return sparse.linalg.spsolve(self.matrix(h), rhs)

    def average_die_temps(self, temperatures):
        return self.observation @ temperatures

    def peak_die_temps(self, temperatures):
        return np.array([np.max(temperatures[nodes]) for nodes in self.chip_top_nodes])

    def boundary_heat(self, temperature, inlet, h):
        return float(np.dot(h * self.area, temperature - inlet))


def build_fem(spec=Spec(), nx=14, ny=None, n_copper=3, n_tim=1, n_die=2, x_grid=None, y_grid=None):
    """Mesh contiguous cuboids into conforming 6-tet cells; integrate exactly.

    Geometric interfaces align to z-grid planes. Material coefficients are
    piecewise constant in each tetra, selected by its centroid.
    """
    if ny is None:
        ny = nx
    if min(nx, ny) < 6 or min(n_copper, n_tim, n_die) < 1:
        raise ValueError("Mesh too coarse")
    x = np.asarray(x_grid,dtype=float) if x_grid is not None else np.linspace(-spec.side / 2, spec.side / 2, nx + 1)
    y = np.asarray(y_grid,dtype=float) if y_grid is not None else np.linspace(-spec.side / 2, spec.side / 2, ny + 1)
    if len(x)!=nx+1 or len(y)!=ny+1 or np.any(np.diff(x)<=0) or np.any(np.diff(y)<=0):
        raise ValueError("Bad custom coordinate grids")
    z = np.concatenate([
        np.linspace(0, spec.z_tim, n_copper + 1),
        np.linspace(spec.z_tim, spec.z_die, n_tim + 1)[1:],
        np.linspace(spec.z_die, spec.height, n_die + 1)[1:],
    ])
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    coords = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    vertex = np.arange(coords.shape[0]).reshape((nx + 1, ny + 1, len(z)))
    offsets = np.array([[0,0,0], [1,0,0], [0,1,0], [1,1,0],
                        [0,0,1], [1,0,1], [0,1,1], [1,1,1]])
    cubes = np.indices((nx, ny, len(z)-1)).reshape(3,-1).T
    v = vertex[cubes[:, None, 0] + offsets[None, :, 0],
               cubes[:, None, 1] + offsets[None, :, 1],
               cubes[:, None, 2] + offsets[None, :, 2]]
    patterns = np.array([[0,1,3,7],[0,3,2,7],[0,2,6,7],
                         [0,6,4,7],[0,4,5,7],[0,5,1,7]])
    tets = v[:, patterns].reshape((-1, 4))
    vertices = coords[tets]
    J = np.transpose(vertices[:,1:,:] - vertices[:,:1,:], (0,2,1))
    det = np.linalg.det(J)
    volumes = np.abs(det) / 6
    if np.any(volumes <= 0):
        raise ValueError("Degenerate tetrahedron")
    invJ = np.linalg.inv(J)
    grad = np.empty((len(tets),4,3))
    grad[:,1:,:] = invJ
    grad[:,0,:] = -np.sum(invJ, axis=1)
    centroids = np.mean(vertices,axis=1)
    zz = centroids[:,2]
    cu = zz < spec.z_tim - 1e-12
    tim = (zz >= spec.z_tim - 1e-12) & (zz < spec.z_die - 1e-12)
    die_zone = zz >= spec.z_die - 1e-12
    k = np.where(cu, spec.k_copper, np.where(tim, spec.k_tim, spec.k_mold)).astype(float)
    cv = np.where(cu, spec.volumetric_heat_copper,
                  np.where(tim, spec.volumetric_heat_tim, spec.volumetric_heat_mold)).astype(float)
    chip = np.full(len(tets), -1, dtype=int)
    for i, (xc,yc) in enumerate(chip_centers(spec)):
        on_chip = die_zone & (np.abs(centroids[:,0]-xc) < spec.chip_side/2) & (np.abs(centroids[:,1]-yc) < spec.chip_side/2)
        chip[on_chip] = i
    valid = chip >= 0
    k[valid] = spec.k_silicon
    cv[valid] = spec.volumetric_heat_silicon
    # Element stiffness and lumped thermal capacity
    Ke = volumes[:,None,None] * k[:,None,None] * np.einsum("nia,nja->nij",grad,grad)
    rows = np.broadcast_to(tets[:,:,None],Ke.shape).reshape(-1)
    cols = np.broadcast_to(tets[:,None,:],Ke.shape).reshape(-1)
    K = sparse.coo_matrix((Ke.ravel(), (rows, cols)), shape=(len(coords),len(coords))).tocsr()
    C = np.bincount(tets.ravel(), weights=np.repeat(volumes*cv/4,4), minlength=len(coords))
    # Lumped Robin boundary surface term: integrate each square exactly for constants.
    area = np.zeros(len(coords))
    for i in range(nx):
        for j in range(ny):
            a = (x[i+1]-x[i])*(y[j+1]-y[j]) / 4
            area[vertex[i:i+2,j:j+2,0].ravel()] += a
    # Heat generation only within silicon. Normalize EACH spatial hotspot
    # discretization to exactly 1 W so all meshes conserve total input power.
    Q = np.zeros((len(coords),4))
    source_volume = np.zeros(4)
    for i, (xc,yc) in enumerate(chip_centers(spec)):
        mask = chip == i
        source_volume[i] = np.sum(volumes[mask])
        if source_volume[i] <= 0:
            raise ValueError(f"No die elements for chip {i}; refine mesh")
        # A displaced Gaussian hotspot plus background load.
        dx = centroids[mask,0] - (xc + (0.0018 if i%2==0 else -0.0018))
        dy = centroids[mask,1] - (yc + (0.0018 if i<2 else -0.0018))
        density = 0.35 + 0.65*np.exp(-(dx**2 + dy**2)/(2*(0.0021)**2))*5
        elem_weights = density*volumes[mask]
        elem_weights /= np.sum(elem_weights)
        Q[:,i] = np.bincount(tets[mask].ravel(),weights=np.repeat(elem_weights/4,4), minlength=len(coords))
    # Chip top-face averages via *area-weighted quadrature*, unlike the old
    # unweighted nodal mean. Each fully silicon top quad contributes A/4 to
    # its four corners (consistent integration for bilinear interpolation).
    top = np.isclose(coords[:,2],spec.height,atol=1e-10)
    observation = np.zeros((4,len(coords)))
    top_nodes=[]
    for i,(xc,yc) in enumerate(chip_centers(spec)):
        mask = top & (np.abs(coords[:,0]-xc) <= spec.chip_side/2+1e-12) & (np.abs(coords[:,1]-yc) <= spec.chip_side/2+1e-12)
        n = np.where(mask)[0]
        if not len(n):
            raise ValueError("Missing top observation nodes")
        top_nodes.append(n)
        w=np.zeros(len(coords))
        # Use actual mesh quadrilaterals; exclude cells crossing the chip edge
        # in non-aligned legacy meshes (mesh-study always uses aligned grids).
        for ix in range(nx):
            for iy in range(ny):
                xmid=.5*(x[ix]+x[ix+1]);ymid=.5*(y[iy]+y[iy+1])
                if abs(xmid-xc) < spec.chip_side/2 and abs(ymid-yc) < spec.chip_side/2:
                    a=(x[ix+1]-x[ix])*(y[iy+1]-y[iy])/4
                    w[vertex[ix:ix+2,iy:iy+2,-1].ravel()] += a
        if w.sum() <= 0:raise ValueError("Top-face area quadrature empty")
        observation[i]=w/w.sum()
    volumes_by_material = {
        "copper":float(np.sum(volumes[cu])),
        "TIM":float(np.sum(volumes[tim])),
        "silicon":float(np.sum(volumes[valid])),
        "mold":float(np.sum(volumes[die_zone & ~valid])),
    }
    return FEM(coords,tets,K,C,area,Q,observation,top_nodes,volumes_by_material,source_volume)


def run_fem(fem, times, powers, inlet=30., h=26000.):
    """Implicit Euler for fixed h and fixed inlet. `powers` has shape (steps,4)."""
    times = np.asarray(times, dtype=float)
    powers = np.asarray(powers, dtype=float)
    if len(times)<2 or powers.shape != (len(times),4):
        raise ValueError("times/powers dimensions")
    dt = np.diff(times)
    if not np.allclose(dt,dt[0]) or dt[0] <= 0:
        raise ValueError("Uniform, increasing time samples required")
    dt = dt[0]
    A = fem.matrix(h) + sparse.diags(fem.C/dt)
    solve = factorized(A.tocsc())
    u = np.full(len(fem.coords), inlet, dtype=float)
    averages = np.zeros((len(times),4))
    peaks = np.zeros_like(averages)
    averages[0] = fem.average_die_temps(u)
    peaks[0] = fem.peak_die_temps(u)
    for j in range(1,len(times)):
        rhs = fem.Q @ powers[j] + (fem.C/dt)*u + h*fem.area*inlet
        u = solve(rhs)
        averages[j] = fem.average_die_temps(u)
        peaks[j] = fem.peak_die_temps(u)
    return {"times":times,"means":averages,"peaks":peaks,"last_field":u}


def build_aligned_fem(spec=Spec(), refinement=1):
    """Nested 3D mesh, with exact chip edges and layer interfaces at all levels.

    Refinements 1, 2, 4 have constant spatial refinement ratio r=2 in x,y,z.
    Copper/TIM/die layers also refine; thermal and power-map geometry is
    identical at every mesh. Robin and heat-capacity lumping still apply.
    """
    if not isinstance(refinement,int) or refinement<1:
        raise ValueError("Positive integer refinement required")
    p=spec.chip_pitch/2; w=spec.chip_side/2
    edges=np.array([-spec.side/2,-p-w,-p+w,p-w,p+w,spec.side/2])
    # 2,4,2,4,2 segments at the coarsest level; interfaces remain exact.
    splits=np.array([2,4,2,4,2])*refinement
    grids=[]
    for a,b,n in zip(edges[:-1],edges[1:],splits):
        grids.extend(np.linspace(a,b,int(n)+1)[:-1])
    x=np.r_[grids,edges[-1]]
    return build_fem(spec,nx=len(x)-1,ny=len(x)-1,
                     n_copper=2*refinement,n_tim=refinement,n_die=2*refinement,
                     x_grid=x,y_grid=x)
