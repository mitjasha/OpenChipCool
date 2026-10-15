"""Independent mesh-block face bookkeeping (NOT an OpenFOAM checkMesh substitute)."""
from collections import Counter,defaultdict
import re
from .geometry import Geometry,mesh_dictionary

FACE_INDICES=((0,3,7,4),(1,5,6,2),(0,4,5,1),(3,2,6,7),
              (0,1,2,3),(4,7,6,5))


def audit_mesh_blocks(g:Geometry)->dict:
    s=mesh_dictionary(g)
    matches=re.findall(r'hex\s*\((\d+(?:\s+\d+){7})\)\s+(water|copper)\s+\((\d+)\s+(\d+)\s+(\d+)\)',s)
    if len(matches)!=2*(2*g.n_channels+1):
        raise RuntimeError('Missing blocks in generated mesh')
    allfaces=defaultdict(list)
    for i,(verts,zone,*_) in enumerate(matches):
        v=list(map(int,verts.split()))
        for ix in FACE_INDICES:
            f=tuple(sorted(v[j] for j in ix))
            allfaces[f].append(zone)
    if any(len(x)>2 for x in allfaces.values()):
        raise RuntimeError('Non-manifold three+ block face detected')
    face_patches={}
    for name in ('inlet','outlet','copperEnds','sideWalls','heater','cover'):
        pat=rf'\b{name}\s*\{{\s*type\s+\w+;\s*faces\s*\((.*?)\)\s*;\s*\}}'
        match=re.search(pat,s,re.S)
        if not match: raise RuntimeError(f'Missing external patch {name}')
        faces=re.findall(r'\((\d+(?:\s+\d+){3})\)',match.group(1))
        for x in faces:
            f=tuple(sorted(map(int,x.split())))
            if f in face_patches: raise RuntimeError(f'Double-labeled patch face {f}')
            face_patches[f]=name
    exterior={f for f,z in allfaces.items() if len(z)==1}
    if exterior!=set(face_patches):
        raise RuntimeError(f'Incorrect patch coverage: missing {len(exterior-set(face_patches))}, extra {len(set(face_patches)-exterior)}')
    interf=sum(1 for vals in allfaces.values() if len(vals)==2 and set(vals)=={'water','copper'})
    if interf!=3*g.n_channels:
        raise RuntimeError(f'Expected 3N block-face interfaces, got {interf}')
    return dict(status='block_topology_combinatorics_pass_NOT_OpenFOAM_checkMesh',
                blocks=len(matches),external_faces=len(exterior),internal_faces=len(allfaces)-len(exterior),
                copper_water_interface_block_faces=interf,
                patches={x:sum(v==x for v in face_patches.values()) for x in ('inlet','outlet','copperEnds','sideWalls','heater','cover')},
                warning='Does not check geometric normals, Jacobians, mesh non-orthogonality or solver coupling')
