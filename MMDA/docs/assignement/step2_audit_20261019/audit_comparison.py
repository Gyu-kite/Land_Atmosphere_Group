"""Read-only audit of one real monthly DA comparison; write local evidence only."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import nbformat

REPO = Path(__file__).resolve().parents[3]
ROOT = REPO / 'notebooks/MMDA_paper_rebuild_20261004/main/new_20261018/unified_scaled_daily'
CACHE = REPO / 'notebooks/MMDA_longterm_qcv4_figures/cache'
OUT = Path(__file__).parent
TARGET = ROOT / 'monthly_sm_climatology/smap_skill'
METHODS = ['OL', 'DA_single', 'DA_multi']
NAMES = ['ISMN'] + METHODS + ['SMAP']

def unique(pattern):
    files = list(CACHE.glob(pattern))
    if len(files) != 1:
        raise ValueError(f'Ambiguous or missing input: {pattern}: {files}')
    return files[0]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def anomaly(frame):
    """Independently reproduce the documented all-year circular 31-day mean."""
    ix = frame.index
    doy = ix.dayofyear.to_numpy() - 1
    doy -= (ix.is_leap_year & (ix.month > 2)).astype(int)
    leap = (ix.month == 2) & (ix.day == 29)
    clim = frame.loc[~leap].groupby(doy[~leap]).mean().reindex(range(365))
    smooth = np.full(clim.shape, np.nan)
    for d in range(365):
        window = clim.iloc[(np.arange(d-15, d+16) % 365)].to_numpy()
        counts = np.isfinite(window).sum(axis=0)
        smooth[d] = np.divide(np.nansum(window, axis=0), counts,
                              out=np.full(frame.shape[1], np.nan), where=counts >= 15)
    expected = smooth[doy].copy()
    expected[leap] = .5*(smooth[58]+smooth[59])
    return frame - expected

paths = {m: unique(pattern) for m, pattern in {
    'OL': 'OL_daily_common8_*.pkl', 'DA_single': 'EnKF_daily_common8_*.pkl',
    'DA_multi': 'MMDA_daily_common8_*.pkl', 'ISMN': 'ismn_qcv4_daily_good_*.pkl',
    'SMAP': 'smap_raw_scaled_daily_*.pkl'}.items()}
paths.update(stations=ROOT/'tables/station_metrics.csv',
             monthly=ROOT/'tables/station_year_month_metrics.csv',
             membership=ROOT/'monthly_sm_climatology/station_membership.csv',
             output=TARGET/'monthly_same_date_skill.csv',
             notebook=TARGET/'figures.ipynb', executed=TARGET/'figures.executed.ipynb',
             source=TARGET/'figures.py',
             cache_builder=REPO/'notebooks/mmda_longterm_qcv4_paper.py')
before = {k: sha(p) for k,p in paths.items()}
data = {m: pd.read_pickle(paths[m]) for m in METHODS+['ISMN']}
pair = pd.read_pickle(paths['SMAP'])
data['SMAP'] = pair['scaled']
stations = pd.read_csv(paths['stations'])
stored = pd.read_csv(paths['monthly'])
stored = stored.loc[stored.obs_only].set_index(['key','year','month'])
members = pd.read_csv(paths['membership'])
assert not stations.network.str.upper().eq('SNOTEL').any()
assert not stored.index.duplicated().any()
checks, date_rows = [], []
errors = {c: 0. for c in ['ACC', 'ubRMSE', 'TE_closure', 'screen_margin', 'group_output']}
recalculated = []
for number, s in enumerate(stations.itertuples()):
    raw = pd.DataFrame({m: data[m][s.key] for m in ['ISMN']+METHODS})
    assert raw.index.is_unique and raw.index.is_monotonic_increasing
    raw = raw.where(raw.notna().all(axis=1), np.nan)
    obs = data['SMAP'][s.key].reindex(raw.index).where(raw.ISMN.notna())
    screen = pd.DataFrame({'ISMN': raw.ISMN, 'OL': raw.OL, 'SMAP': obs}).dropna()
    assert len(screen) == s.n_screen
    margin = screen.ISMN.corr(screen.SMAP)-screen.ISMN.corr(screen.OL)
    assert bool(len(screen)>=60 and (screen.std(ddof=0)>0).all() and margin>-.2) == s.filtered
    if np.isfinite(s.margin_scaled):
        errors['screen_margin'] = max(errors['screen_margin'], abs(margin-s.margin_scaled))
    a = anomaly(raw)
    a['SMAP'] = anomaly(obs.to_frame('SMAP')).SMAP
    raw['SMAP'] = obs
    common = a[NAMES].dropna()
    date_rows.append(dict(key=s.key, model_common_days=int(raw[METHODS+['ISMN']].notna().all(axis=1).sum()),
        screen_days=len(screen), common_anomaly_obs_days=len(common),
        start=str(common.index.min().date()), end=str(common.index.max().date())))
    for period, g in common.groupby(common.index.to_period('M')):
        if len(g)<8:
            continue
        z=g-g.mean()
        sig=np.sqrt((z*z).mean())
        if not (sig>0).all():
            continue
        record=dict(key=s.key,year=period.year,month=period.month,N=len(g))
        old=stored.loc[(s.key,period.year,period.month)]
        assert int(old.N)==len(g)
        for m in METHODS+['SMAP']:
            acc=float(g[m].corr(g.ISMN))
            ub=float((raw.loc[g.index,m]-raw.loc[g.index,'ISMN']).std(ddof=0))
            ve=(sig[m]-sig.ISMN)**2
            ce=2*sig[m]*sig.ISMN*(1-acc)
            te=float(((z[m]-z.ISMN)**2).mean())
            errors['ACC']=max(errors['ACC'],abs(acc-old['ACC_'+m]))
            errors['ubRMSE']=max(errors['ubRMSE'],abs(ub-old['ubRMSE_'+m]))
            errors['TE_closure']=max(errors['TE_closure'],abs(ve+ce-te))
            record['ACC_'+m]=acc
            record['ubRMSE_'+m]=ub
        recalculated.append(record)
    if number%100==0:
        print('Daily reconstruction',number+1,'/',len(stations),flush=True)
rec=pd.DataFrame(recalculated)
assert len(rec)==len(stored)
rec.to_csv(OUT/'recalculated_station_year_month.csv',index=False)
pd.DataFrame(date_rows).to_csv(OUT/'date_support.csv',index=False)
cols=[f'{v}_{m}' for v in ['ACC','ubRMSE'] for m in METHODS+['SMAP']]
output=pd.read_csv(paths['output'])
group_counts=[]
for subset in ['all','filtered']:
    membership=members.loc[members.subset==subset]
    base=stations.loc[(stations.lc_code>0)&(stations.filtered if subset=='filtered' else True)]
    assert set(base.key)==set(membership.key)
    q1,q2=np.quantile(base.OL_mean_SM,[1/3,2/3])
    ref=base.assign(sm_group=np.where(base.OL_mean_SM<q1,0,np.where(base.OL_mean_SM<q2,1,2)))
    joined=ref[['key','sm_group']].merge(membership,on='key',suffixes=('_expected',''))
    assert (joined.sm_group_expected==joined.sm_group).all()
    q=rec.merge(membership[['key','sm_group']],on='key',validate='many_to_one')
    for group in [-1,0,1,2]:
        sample=q if group==-1 else q.loc[q.sm_group==group]
        sm=sample.groupby(['key','month'])[cols].mean().reset_index()
        avg=sm.groupby('month')[cols].mean()
        counts=sm.groupby('month').key.nunique()
        target=output.loc[(output.subset==subset)&(output.sm_group==group)].set_index('month')
        assert np.array_equal(target.N.to_numpy(),counts.reindex(target.index).to_numpy())
        delta=np.max(abs(avg[cols].to_numpy()-target.loc[avg.index,cols].to_numpy()))
        errors['group_output']=max(errors['group_output'],float(delta))
        group_counts.append(dict(subset=subset,sm_group=group,total_sites=sample.key.nunique(),
                                 minimum_month_N=int(counts.min()),maximum_month_N=int(counts.max())))
for value in errors.values():
    assert value<1e-10,errors
for key in ['notebook','executed']:
    nb=nbformat.read(paths[key],as_version=4)
    assert not [o for cell in nb.cells for o in cell.get('outputs',[]) if o.output_type=='error']
after={k:sha(p) for k,p in paths.items()}
assert before==after
pd.DataFrame(group_counts).to_csv(OUT/'sample_counts.csv',index=False)
audit=dict(input_files={k:dict(path=str(paths[k]),sha256=v) for k,v in before.items()},
    output_row_count=len(output),station_count=len(stations),filtered_count=int(stations.filtered.sum()),
    reconstructed_observation_months=len(rec),max_absolute_errors=errors,
    original_inputs_unchanged=True,
    limits=['No L2/DAOBS byte-level rebuild; existing daily caches are trusted inputs.',
            'Executed notebook was inspected, not rerun; its numerical output was independently recomputed.',
            'All-available-year climatology is descriptive evaluation, not out-of-sample forecast validation.',
            'Shared scaled SMAP is sourced from DA_single DAOBS, not separately from DA_multi.'])
(OUT/'audit.json').write_text(json.dumps(audit,indent=2))
print(json.dumps({k:v for k,v in audit.items() if k!='input_files'},indent=2))
