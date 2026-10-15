#!/usr/bin/env python3
"""Plots manually estimated PAPER symbols and analytical energy — NEVER CFD."""
import csv
from pathlib import Path
import matplotlib.pyplot as plt
from openchipcool.cfd.reference import experimental_condition

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/cfd'; OUT.mkdir(exist_ok=True,parents=True)

def load(name):
    with (ROOT/'data/literature'/name).open(newline='') as f:
        return list(csv.DictReader(f))

fig5=load('qu_mudawar_2002_tout_fig5_visual_estimate.csv')
fig,ax=plt.subplots(figsize=(8.7,5.4))
for q in (100.,200.):
    sub=sorted((r for r in fig5 if float(r['qeff_W_cm2'])==q),key=lambda t:float(t['Re']))
    re=[float(r['Re']) for r in sub]
    vals=[float(r['delta_T_K_visual_approx']) for r in sub]
    uncertainty=[float(r['visual_reading_uncertainty_K']) for r in sub]
    ax.errorbar(re, vals,yerr=uncertainty,fmt='o',capsize=3,label=f'Fig. 5 symbols, estimated; {q:g} W/cm²')
    xx=[min(re)+i*(max(re)-min(re))/150 for i in range(151)]
    yy=[experimental_condition(x,q)['adiabatic_outlet_rise_K'] for x in xx]
    ax.plot(xx, yy, linestyle='--',label=f'Energy balance ONLY; {q:g} W/cm²')
ax.set_xlabel('Reynolds number, inlet')
ax.set_ylabel('Coolant outlet − inlet temperature (K)')
ax.set_title('Qu & Mudawar (2002), Fig. 5 — PAPER SYMBOLS + ENERGY BALANCE\nNOT CFD; graphical reading ±2–5 K (indicative)')
ax.legend(fontsize=8);ax.grid(alpha=.25);fig.tight_layout()
fig.savefig(OUT/'fig5_reference_energy_NOT_CFD.png',dpi=160)
plt.close(fig)

fig7=load('qu_mudawar_2002_tc_fig7_visual_estimate.csv')
for q in (100.,200.):
    fig,ax=plt.subplots(figsize=(8.7,5.2))
    modes=sorted({int(float(r['Re'])) for r in fig7 if float(r['qeff_W_cm2'])==q})
    for re in modes:
        a=sorted((r for r in fig7 if float(r['qeff_W_cm2'])==q and float(r['Re'])==re),
                 key=lambda r:int(r['thermocouple_number']))
        ax.errorbar([float(x['x_mm']) for x in a],
                    [float(x['Ttc_rise_K_visual_approx']) for x in a],
                    yerr=[float(x['visual_reading_uncertainty_K']) for x in a],
                    marker='o',capsize=2,label=f'Fig. 7 markers estimated Re={re}')
    ax.set_title(f'Qu & Mudawar (2002), Fig. 7 — {q:g} W/cm²\nPAPER SYMBOLS ONLY, visually estimated ±5–6 K; NOT CFD')
    ax.set_xlabel('Thermocouple x (mm)')
    ax.set_ylabel('Copper thermocouple temperature − inlet temperature (K)')
    ax.grid(alpha=.25);ax.legend(fontsize=9);fig.tight_layout()
    fig.savefig(OUT/f'fig7_reference_q{int(q)}_NOT_CFD.png',dpi=160)
    plt.close(fig)
print('Plotted PAPER graph estimates / energy balance only; NO CFD results.')
