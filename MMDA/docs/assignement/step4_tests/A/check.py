from pathlib import Path
import pandas as pd,numpy as np,json
B=Path('/land1/user/gychoi/LIS/LISF_MMDA/notebooks');R=B/'MMDA_paper_rebuild_20261004/main/new_20261018/unified_scaled_daily';C=B/'MMDA_longterm_qcv4_figures/cache'
M=pd.read_csv(R/'tables/station_year_month_metrics.csv');S=pd.read_csv(R/'tables/station_metrics.csv');G=pd.read_csv(R/'monthly_sm_climatology/station_membership.csv');T=pd.read_csv(R/'monthly_sm_climatology/smap_skill/monthly_same_date_skill.csv')
D={n:pd.read_pickle(next(C.glob(p))) for n,p in {'OL':'OL_daily_common8*','DA_single':'EnKF_daily_common8*','DA_multi':'MMDA_daily_common8*','ISMN':'ismn_qcv4_daily_good*'}.items()};O=pd.read_pickle(next(C.glob('smap_raw_scaled_daily*')))['scaled']
def anomaly(f):
 ix=f.index;d=ix.dayofyear.to_numpy()-1-(ix.is_leap_year&(ix.month>2));keep=~((ix.month==2)&(ix.day==29)); b=f[keep].groupby(d[keep]).mean().reindex(range(365));v=b.to_numpy();sm=np.full(v.shape,np.nan)
 for j in range(365):
  w=v[[(j+h)%365 for h in range(-15,16)]];num=np.isfinite(w).sum(axis=0);sm[j]=np.divide(np.nansum(w,axis=0),num,out=np.full(v.shape[1],np.nan),where=num>=15)
 expected=sm[d].copy();leap=(ix.month==2)&(ix.day==29);expected[leap]=(sm[58]+sm[59])/2
 return f-expected
errs=[];nerr=[];recs=[];dates=[];screenerr=[];group_err=[]
for site in S.itertuples():
 raw=pd.DataFrame({k:v[site.key] for k,v in D.items()});raw.loc[raw.isna().any(axis=1),:]=np.nan
 obs=O[site.key].reindex(raw.index).where(raw.ISMN.notna());a=anomaly(raw);raw['SMAP']=obs;a['SMAP']=anomaly(obs.to_frame()).iloc[:,0]
 sc=raw[['ISMN','OL','SMAP']].dropna();valid=len(sc)>=60 and (sc.std(ddof=0)>0).all();f=bool(valid and sc.ISMN.corr(sc.SMAP)-sc.ISMN.corr(sc.OL)>-.2);screenerr.append(f!=site.filtered)
 for period,q in a.dropna().groupby(a.dropna().index.to_period('M')):
  if len(q)<8 or (q.std(ddof=0)<=0).any():continue
  row=dict(key=site.key,year=period.year,month=period.month,N=len(q));old=M[(M.key==site.key)&M.obs_only&(M.year==period.year)&(M.month==period.month)];assert len(old)==1
  nerr.append(len(q)-old.iloc[0].N);dates.extend(q.index)
  for method in ['OL','DA_single','DA_multi','SMAP']:
   row['ACC_'+method]=np.corrcoef(q.ISMN,q[method])[0,1];e=(raw.loc[q.index,method]-raw.loc[q.index,'ISMN']).to_numpy();row['ubRMSE_'+method]=np.sqrt(np.mean((e-e.mean())**2))
   for metric in ['ACC','ubRMSE']:errs.append(abs(row[metric+'_'+method]-old.iloc[0][metric+'_'+method]))
  recs.append(row)
Q=pd.DataFrame(recs);saved=M[M.obs_only];assert set(map(tuple,Q[['key','year','month']].values))==set(map(tuple,saved[['key','year','month']].values))
aggerr=[];stats={}
for subset in ['all','filtered']:
 g=G[G.subset==subset];st=S[(S.lc_code>0)&(S.filtered if subset=='filtered' else True)];qs=st.OL_mean_SM.quantile([1/3,2/3]).values
 assert set(g.key)==set(st.key)
 for r in g.itertuples():
  v=st.set_index('key').loc[r.key,'OL_mean_SM'];group_err.append(r.sm_group!=int(np.searchsorted(qs,v,side='right')))
 sample=Q[Q.key.isin(g.key)];sets=[set(sample[sample.month==m].key) for m in range(1,13)]
 stats[subset]={'membership':len(g),'monthly_N':[len(x) for x in sets],'fixed_all_months':len(set.intersection(*sets)),'adjacent_symmetric_difference':[len(sets[i]^sets[i+1]) for i in range(11)]}
 for group in [-1,0,1,2]:
  keys=g.key if group==-1 else g[g.sm_group==group].key;q=Q[Q.key.isin(keys)];cols=[c for c in Q if c.startswith(('ACC_','ubRMSE_'))];a=q.groupby(['key','month'])[cols].mean().groupby('month').mean()
  for method in ['SMAP','DA_single','DA_multi']:
   a['dACC_'+method]=a['ACC_'+method]-a.ACC_OL;a['dubRMSE_'+method]=a.ubRMSE_OL-a['ubRMSE_'+method]
  t=T[(T.subset==subset)&(T.sm_group==group)].set_index('month');aggerr.extend(np.abs(a-t[a.columns]).to_numpy().ravel())
result={'monthly_rows':len(Q),'date_min':str(min(dates)),'date_max':str(max(dates)),'N_range':[int(Q.N.min()),int(Q.N.max())],'N_max_error':int(max(map(abs,nerr))),'metric_max_error':float(max(errs)),'aggregate_max_error':float(np.nanmax(aggerr)),'screen_mismatches':sum(screenerr),'group_mismatches':sum(group_err),'composition':stats}
print(json.dumps(result,indent=2));Path('/tmp/mmda_fresh_audit/result.json').write_text(json.dumps(result,indent=2))
