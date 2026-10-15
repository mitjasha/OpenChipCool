"""Fail-closed experimental evidence checks for Qu & Mudawar 2002.

Nothing here runs CFD. The only accepted quantitative comparisons are from an
explicitly supplied *independent solver export* plus labelled approximate
published reference observations. Missing observables remain missing.
"""
from __future__ import annotations
from pathlib import Path
import csv, json, math
from .reference import experimental_condition

EXPERIMENT_DOI = '10.1016/S0017-9310(01)00337-4'
TC_X_MM = (5., 16.588, 28.176, 39.764)


def read_csv(path: Path) -> list[dict]:
    with Path(path).open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))


def pressure_budget(*, rho: float, u_channel: float, u_p2_in: float,
                    u_p2_out: float, u_p1_in: float, u_p1_out: float,
                    Kc1: float=1., Kc2: float=1., Ke2: float=1., Ke1: float=1.) -> dict:
    """Full Qu & Mudawar Eq. 18-22 pressure budget, pressure *loss* terms.

    The four plenum velocities ARE REQUIRED: the article does not provide a
    unique stand-alone geometric CAD model of both actual manifolds. Supply
    independent geometric estimates; do not choose areas to fit Fig. 4.
    Note kinetic energy differences can be negative but total losses should
    be audited. Channel loss is supplied independently by actual CFD.
    """
    nums = (rho,u_channel,u_p2_in,u_p2_out,u_p1_in,u_p1_out,Kc1,Kc2,Ke2,Ke1)
    if not all(math.isfinite(x) for x in nums) or rho <= 0 or min(nums[1:]) < 0:
        raise ValueError('Invalid physical parameters')
    return {
        'contraction_deep_to_shallow_Pa': .5*rho*(u_p2_in**2-u_p1_in**2+Kc1*u_p2_in**2),
        'contraction_shallow_to_channel_Pa': .5*rho*(u_channel**2-u_p2_in**2+Kc2*u_channel**2),
        'expansion_channel_to_shallow_Pa': .5*rho*(u_p2_out**2-u_channel**2+Ke2*u_channel**2),
        'expansion_shallow_to_deep_Pa': .5*rho*(u_p1_out**2-u_p2_out**2+Ke1*u_p2_out**2),
    }


def grid_convergence(coarse: float, medium: float, fine: float, r: float=2., fs: float=1.25) -> dict:
    """Three-grid Richardson/GCI estimate for uniform linear refinement.

    Valid only when solutions are non-oscillatory and in asymptotic regime.
    Does NOT prove asymptotic convergence; record actual cell sizes and repeat.
    """
    if r <= 1: raise ValueError('Need refinement ratio > 1')
    d32=coarse-medium
    d21=medium-fine
    if not math.isfinite(d21*d32) or d21*d32 <= 0:
        raise ValueError('No monotone convergence; do not report a GCI')
    if not abs(d21) or not abs(fine):
        raise ValueError('Zero difference or zero fine value')
    ratio=abs(d32/d21)
    if ratio <= 1:
        raise ValueError('No useful apparent convergence order (error not decreasing)')
    p=math.log(ratio)/math.log(r)
    extrapolated=fine+(fine-medium)/(r**p-1)
    gci_fine=fs*abs((fine-medium)/fine)/(r**p-1)
    return {'observed_order':p,'richardson_extrapolated':extrapolated,
            'GCI_fine_fraction':gci_fine,'GCI_fine_percent':100*gci_fine,
            'apparent_order_only':True,'refinement_ratio':r}


def validate_export(export_csv: Path, repo_root: Path) -> dict:
    """Review actual user-supplied CFD metrics against figure observations.

    The export must carry check flags and independent provenance to prevent
    mistaking a pressure-only analytical model for solved CHT.
    """
    rows=read_csv(export_csv)
    if not rows: raise ValueError('CFD export contains no rows')
    required={'case_id','mesh','Re','qeff_W_cm2','dp_core_Pa','dp_plenums_Pa',
              'Tin_C','Tout_massweighted_C','Ttc1_rise_K','Ttc2_rise_K',
              'Ttc3_rise_K','Ttc4_rise_K','mass_error_frac','energy_error_frac',
              'solver_log_path','mesh_water_check_path','mesh_copper_check_path'}
    if not required.issubset(rows[0]):
        raise ValueError('Missing required independently generated CFD columns: '+', '.join(sorted(required-set(rows[0]))))
    # A nonempty log filename is not proof of completion; check actual logs.
    out=[]
    ref_pressure=read_csv(Path(repo_root)/'data/literature/qu_mudawar_2002_pressure_approx.csv')
    ref_tc=read_csv(Path(repo_root)/'data/literature/qu_mudawar_2002_tc_fig7_visual_estimate.csv')
    ref_out=read_csv(Path(repo_root)/'data/literature/qu_mudawar_2002_tout_fig5_visual_estimate.csv')
    for row in rows:
        case=Path(row['case_id'])
        for fld in ('solver_log_path','mesh_water_check_path','mesh_copper_check_path'):
            log=Path(row[fld]); log=log if log.is_absolute() else Path(repo_root)/log
            if not log.is_file() or log.stat().st_size < 10:
                raise FileNotFoundError(f'Cannot establish real solver evidence: {fld} {log}')
            body=log.read_text(errors='replace')
            if fld=='solver_log_path':
                if '\nEnd\n' not in body and not body.rstrip().endswith('End'):
                    raise ValueError(f'OpenFOAM solver does not have completed End marker in {log}')
            elif 'Mesh OK' not in body:
                raise ValueError(f'Mesh check missing Mesh OK in {log}')
        q=float(row['qeff_W_cm2']); re=float(row['Re']); tin=float(row['Tin_C'])
        if not (100<=q<=200 and 130<=re<=1700):
            raise ValueError('Re/flux outside comparison domain')
        if float(row['mass_error_frac'])>=.001 or float(row['energy_error_frac'])>=.01:
            raise ValueError('Mass/energy closure fail (0.1%/1% proposed limits)')
        candidates=[x for x in ref_pressure if abs(float(x['qeff_W_cm2'])-q)<.1]
        if not candidates: raise ValueError('Missing Fig4 flux set')
        candidates.sort(key=lambda x:float(x['Re']))
        if re < float(candidates[0]['Re']) or re > float(candidates[-1]['Re']):
            raise ValueError('Pressure comparison outside digitized Re range')
        def lerp(a, b, key, re):
            ra=float(a['Re']); rb=float(b['Re'])
            if ra==rb: return float(a[key])
            return float(a[key])+(float(b[key])-float(a[key]))*(re-ra)/(rb-ra)
        for a,b in zip(candidates[:-1],candidates[1:]):
            if float(a['Re'])<=re<=float(b['Re']):
                dp_exp=lerp(a,b,'dp_bar_visual_approx',re)*1e5
                break
        dp_cfd=float(row['dp_core_Pa'])+float(row['dp_plenums_Pa'])
        if dp_cfd<=0 or float(row['dp_core_Pa'])<=0: raise ValueError('Invalid pressure budget')
        tc_candidates=[v for v in ref_tc if abs(float(v['qeff_W_cm2'])-q)<.1 and abs(float(v['Re'])-re)<=2]
        if len(tc_candidates)!=4:
            raise ValueError(f'Four independent Fig7 observations missing for q={q}, Re={re}')
        err_tc=[]
        for v in tc_candidates:
            k=int(v['thermocouple_number'])
            err_tc.append(float(row[f'Ttc{k}_rise_K'])-float(v['Ttc_rise_K_visual_approx']))
        tout_candidates=[v for v in ref_out if abs(float(v['qeff_W_cm2'])-q)<.1 and abs(float(v['Re'])-re)<=2]
        if len(tout_candidates)!=1:
            raise ValueError('Fig5 experimental observation unavailable for this case')
        delta_t=float(row['Tout_massweighted_C'])-tin
        exp_rise=float(tout_candidates[0]['delta_T_K_visual_approx'])
        out.append({'case_id':str(case),'Re':re,'flux_W_cm2':q,
                    'dp_core_Pa':float(row['dp_core_Pa']),
                    'dp_plenum_correction_Pa':float(row['dp_plenums_Pa']),
                    'dp_total_Pa':dp_cfd,'dp_paper_approx_Pa':dp_exp,
                    'dp_relative_error':(dp_cfd-dp_exp)/dp_exp,
                    'delta_T_fem_cfd_K':delta_t,'delta_T_fig5_approx_K':exp_rise,
                    'delta_T_difference_K':delta_t-exp_rise,
                    'thermocouple_error_K':err_tc,
                    'thermocouple_MAE_K':sum(map(abs,err_tc))/4,
                    'reference_status':'VISUAL_DIGITIZATION_UNCERTAIN',
                    'validation_status':'COMPARISON_ONLY_NOT_CERTIFIED'})
    return {'experiment':EXPERIMENT_DOI,'status':'COMPARISON_ONLY_NOT_CERTIFIED',
            'warnings':['Graphical digitizations, not raw experimental data.',
                        'Solver completion log does not establish iterative convergence.',
                        'Missing real fixture details or temperature-dependent viscosity invalidate full-experiment claim.'],
            'comparisons':out}
