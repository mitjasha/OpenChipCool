#!/usr/bin/env python3
"""Draw cross-sections of actual generated parametrized geometries, NOT CFD solutions."""
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from openchipcool.cfd.geometry import geometry_profile,_y_edges

def main():
 root=Path(__file__).resolve().parents[1]
 for nm,title in [('paper','Qu & Mudawar (2002): nominal 21 channels'),
                  ('design36','OpenChipCool proposal: 108 channels, 36 mm')]:
  g=geometry_profile(nm)
  edges,types=_y_edges(g)
  fig,ax=plt.subplots(figsize=(11,3.3),layout='constrained')
  widthmm=g.width_m*1e3; basemm=g.base_m*1e3; dmm=g.channel_depth_m*1e3
  ax.add_patch(Rectangle((0,0),widthmm,basemm,fc='#ad783d',ec='#6c4521',lw=.7))
  for i,kind in enumerate(types):
   xa=edges[i]*1e3;dx=(edges[i+1]-edges[i])*1e3
   ax.add_patch(Rectangle((xa,basemm),dx,dmm,fc='#5da9d1' if kind=='fluid' else '#ad783d',ec='none'))
  ax.set(xlim=(0,widthmm),ylim=(0,basemm+dmm+.4),xlabel='Plate width (mm)',
         ylabel='Height (mm)',title=f'{title} — copper / water cross-section, NOT computed CFD')
  ax.set_aspect('auto')
  fig.savefig(root/'results/cfd'/f'geometry_{nm}.png',dpi=180)
  plt.close(fig)
  print(root/'results/cfd'/f'geometry_{nm}.png')
if __name__=='__main__':main()
