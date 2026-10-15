"""Structured conformal blockMesh generator for 21-channel copper/water CHT.

Generates *geometry and cellZones*, not a computed CFD solution.  The upper
polycarbonate cover is represented by an adiabatic external wall; experimental
inlet/outlet plenums, deep heater block and air slots are absent in this stage.

Source: Qu & Mudawar (2002), IJHMT 45, pp. 2549-2565, Table 1 & Fig. 2.
For the new 36-mm design all parameters (except channel cross section) are
explicit design hypotheses and NOT measurements from this paper.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import json
from pathlib import Path


@dataclass(frozen=True)
class Geometry:
    name: str = 'qu_mudawar_2002_core'
    n_channels: int = 21
    length_m: float = 0.044764  # Table 1, 44.764 mm in unit cell
    width_m: float = 0.010  # 1 cm overall width, Fig. 2
    channel_width_m: float = 231e-6
    channel_depth_m: float = 713e-6
    pitch_m: float = 467e-6  # 118 + 231 + 118 microns
    base_m: float = 5637e-6  # Hw2, Table 1
    tin_K: float = 288.15  # experiment 15 deg C
    nx: int = 20
    ny_channel: int = 3
    ny_wall: int = 2
    nz_base: int = 6
    nz_channel: int = 4

    def validate(self) -> None:
        if self.n_channels < 1 or min(self.length_m, self.width_m, self.channel_width_m,
                                       self.channel_depth_m, self.pitch_m, self.base_m) <= 0:
            raise ValueError('Geometry dimensions and number of channels must be positive')
        if self.channel_width_m >= self.pitch_m:
            raise ValueError('Channel width must be smaller than pitch')
        if self.n_channels*self.pitch_m >= self.width_m:
            raise ValueError('Channel array does not fit in specified plate width')
        if min(self.nx, self.ny_channel, self.ny_wall, self.nz_base, self.nz_channel) < 1:
            raise ValueError('Mesh subdivisions must be positive')


def geometry_profile(name='paper', level='coarse') -> Geometry:
    opts = {
        'coarse': dict(nx=20,ny_channel=3,ny_wall=2,nz_base=6,nz_channel=4),
        'medium': dict(nx=40,ny_channel=6,ny_wall=4,nz_base=12,nz_channel=8),
        'fine': dict(nx=80,ny_channel=12,ny_wall=8,nz_base=24,nz_channel=16),
    }
    if level not in opts:
        raise ValueError(f'Unknown resolution {level!r}: {list(opts)}')
    if name == 'paper':
        return Geometry(**opts[level])
    if name == 'design36':
        # The v0.2 108-channel plate adopted here.  The 36-mm design uses
        # 0.231 mm slots and a *new*, narrower ~0.10 mm web pitch.
        return Geometry(name='openchipcool_36mm_concept',n_channels=108,
                        length_m=.036,width_m=.036,
                        channel_width_m=231e-6,channel_depth_m=713e-6,
                        pitch_m=332e-6,base_m=.003,tin_K=303.15,
                        **opts[level])
    raise ValueError(f'Unknown geometry profile: {name}')


def _y_edges(g: Geometry):
    # center entire periodic microchannel array in the physical 1-cm width;
    # no guess about exact machined edge walls beyond published unit cell.
    margin = (g.width_m-g.n_channels*g.pitch_m)/2
    web_half = (g.pitch_m-g.channel_width_m)/2
    edges=[0.0, margin+web_half]
    kinds=['solid']
    for i in range(g.n_channels):
        edges.append(edges[-1]+g.channel_width_m)
        kinds.append('fluid')
        if i < g.n_channels-1:
            edges.append(edges[-1]+2*web_half)
        else:
            edges.append(g.width_m)
        kinds.append('solid')
    if abs(edges[-1]-g.width_m)>1e-11:
        raise AssertionError('Failed to close width')
    return edges, kinds


def geometry_metrics(g:Geometry)->dict:
    g.validate()
    edges,kinds=_y_edges(g)
    block_count=2*len(kinds)
    xy_cells=sum(g.ny_channel if k=='fluid' else g.ny_wall for k in kinds)
    fluid_cells=g.nx*g.ny_channel*g.nz_channel*g.n_channels
    copper_cells=g.nx*(xy_cells*g.nz_base+(len(kinds)-g.n_channels)*g.ny_wall*g.nz_channel)
    return dict(profile=g.name,channels=g.n_channels,length_m=g.length_m,width_m=g.width_m,
                pitch_m=g.pitch_m,depth_m=g.channel_depth_m,base_m=g.base_m,
                array_width_m=g.n_channels*g.pitch_m,edge_margin_m=(g.width_m-g.n_channels*g.pitch_m)/2,
                hydraulic_diameter_m=2*g.channel_width_m*g.channel_depth_m/(g.channel_width_m+g.channel_depth_m),
                area_channel_m2=g.channel_width_m*g.channel_depth_m,
                planform_area_m2=g.length_m*g.width_m,
                blocks=block_count,cells_fluid=fluid_cells,cells_copper=copper_cells,
                cells_total=fluid_cells+copper_cells,
                water_volume_m3=g.n_channels*g.length_m*g.channel_width_m*g.channel_depth_m,
                copper_volume_m3=g.width_m*g.length_m*(g.base_m+g.channel_depth_m)-
                       g.n_channels*g.length_m*g.channel_width_m*g.channel_depth_m,
                thermocouple_x_m=[.005,.016588,.028176,.039764] if g.name=='qu_mudawar_2002_core' else [],
                thermocouple_z_m=.003175 if g.name=='qu_mudawar_2002_core' else None,
                methodology='Full-length 21-channel CHT core, no plenums or top cover conduction, not full fixture')


def mesh_dictionary(g: Geometry)->str:
    """Produce complete blockMeshDict with 21 internal fluid ducts + connected Cu.

    Geometry shares exact face vertices across zones; every external surface is
    explicitly named, no opaque defaultFaces patches or missing blocks.
    """
    g.validate()
    ys,kinds=_y_edges(g)
    xs=[0.,g.length_m]
    zs=[0.,g.base_m,g.base_m+g.channel_depth_m]
    verts=[(x,y,z) for z in zs for y in ys for x in xs]
    ny=len(ys)
    def idx(x,y,z): return z*ny*2+y*2+x
    def blk(y,z):
        return (idx(0,y,z),idx(1,y,z),idx(1,y+1,z),idx(0,y+1,z),
                idx(0,y,z+1),idx(1,y,z+1),idx(1,y+1,z+1),idx(0,y+1,z+1))
    blocks=[]
    for z in range(2):
        for y,k in enumerate(kinds):
            zone='water' if z==1 and k=='fluid' else 'copper'
            blocks.append((blk(y,z),zone, g.nx,
                           g.ny_channel if k=='fluid' else g.ny_wall,
                           g.nz_base if z==0 else g.nz_channel))
    def f(y,z,face):
        b=blk(y,z)
        return [b[v] for v in {
            'x0':(0,3,7,4),'x1':(1,5,6,2),
            'y0':(0,4,5,1),'y1':(3,2,6,7),
            'z0':(0,1,2,3),'z1':(4,7,6,5),
        }[face]]
    surfaces={'inlet':[],'outlet':[],'copperEnds':[],
              'sideWalls':[],'heater':[],'cover':[]}
    for z in range(2):
        for y,k in enumerate(kinds):
            low=z==1 and k=='fluid'
            surfaces['inlet' if low else 'copperEnds'].append(f(y,z,'x0'))
            surfaces['outlet' if low else 'copperEnds'].append(f(y,z,'x1'))
            if y==0: surfaces['sideWalls'].append(f(y,z,'y0'))
            if y==len(kinds)-1: surfaces['sideWalls'].append(f(y,z,'y1'))
            if z==0: surfaces['heater'].append(f(y,z,'z0'))
            if z==1: surfaces['cover'].append(f(y,z,'z1'))
    out=['FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }',
         'scale 1;', 'vertices', '(']
    out += ['    (%.12g %.12g %.12g) // %d'%(x,y,z,i) for i,(x,y,z) in enumerate(verts)]
    out += [');','blocks','(']
    for b,zone,nx,nyc,nz in blocks:
        out.append('    hex (%s) %s (%d %d %d) simpleGrading (1 1 1)'%(
            ' '.join(map(str,b)), zone,nx,nyc,nz))
    out += [');','edges ();', 'boundary', '(']
    for nm,faces in surfaces.items():
        out += ['    %s { type %s; faces ('%(nm,'patch' if nm in ('inlet','outlet') else 'wall')]
        out += ['        (%s)'%' '.join(map(str,face)) for face in faces]
        out += ['    ); }']
    out += [');','mergePatchPairs ();']
    return '\n'.join(out)+'\n'


def write_geometry_case(g:Geometry, root:Path)->dict:
    root=Path(root)
    (root/'system').mkdir(parents=True,exist_ok=True)
    (root/'constant').mkdir(parents=True,exist_ok=True)
    (root/'system/blockMeshDict').write_text(mesh_dictionary(g))
    (root/'constant/regionProperties').write_text('''FoamFile { version 2.0; format ascii; class dictionary; object regionProperties; }
regions (fluid (water) solid (copper));
''')
    meta=geometry_metrics(g)
    (root/'geometry.json').write_text(json.dumps({**asdict(g), **meta},indent=2)+'\n')
    from .audit import audit_mesh_blocks
    (root/'mesh_audit.json').write_text(json.dumps(audit_mesh_blocks(g),indent=2)+'\n')
    return meta
