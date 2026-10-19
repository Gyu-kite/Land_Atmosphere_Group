from pathlib import Path
import json
import pandas as pd
import numpy as np
R=Path(__file__).resolve().parents[3]
O=Path(__file__).parent
T=R/'notebooks/MMDA_paper_rebuild_20261004/main/new_20261018/unified_scaled_daily'
P=T/'monthly_sm_climatology/smap_skill'
nb=json.loads((P/'figures.ipynb').read_text())
code='\n\n'.join(''.join(c['source']).rstrip() for c in nb['cells'] if c['cell_type']=='code')
py=(P/'figures.py').read_text()
py='\n'.join(l for l in py.splitlines() if not l.startswith('# %%')).strip()
assert code.strip()==py
q=pd.read_csv(O/'recalculated_station_year_month.csv')
m=pd.read_csv(T/'monthly_sm_climatology/station_membership.csv')
f=pd.read_csv(P/'monthly_same_date_skill.csv')
rows=[]; sets=[]
for subset in ['all','filtered']:
    membership=m[m.subset==subset]
    z=q.merge(membership[['key','sm_group']],on='key')
    for g in [-1,0,1,2]:
        a=z if g==-1 else z[z.sm_group==g]
        keys={mo:set(a.loc[a.month==mo,'key']) for mo in range(1,13)}
        common=set.intersection(*keys.values())
        union=set.union(*keys.values())
        sets.append(dict(subset=subset,sm_group=g,union=len(union),all_12_months=len(common)))
        for mo in range(1,13):
            prev=keys[12 if mo==1 else mo-1]
            rows.append(dict(subset=subset,sm_group=g,month=mo,N=len(keys[mo]),added=len(keys[mo]-prev),removed=len(prev-keys[mo]),keys='|'.join(sorted(keys[mo]))))
errors={}
for method in ['SMAP','DA_single','DA_multi']:
    errors['dACC_'+method]=float(abs(f['dACC_'+method]-(f['ACC_'+method]-f.ACC_OL)).max())
    errors['dubRMSE_'+method]=float(abs(f['dubRMSE_'+method]-(f.ubRMSE_OL-f['ubRMSE_'+method])).max())
assert max(errors.values())<1e-12
pd.DataFrame(rows).to_csv(O/'monthly_site_composition.csv',index=False)
dates=pd.read_csv(O/'date_support.csv')
result=dict(notebook_matches_python=True,date_range=[dates.start.min(),dates.end.max()],years=[int(q.year.min()),int(q.year.max())],monthly_day_range=[int(q.N.min()),int(q.N.max())],site_composition=sets,max_delta_errors=errors,raw_ismn_acc_support='ISMN+models common daily mask; SMAP climatology uses its available observation days separately')
(O/'supplement.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
