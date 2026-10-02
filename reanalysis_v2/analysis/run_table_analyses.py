from pathlib import Path
import argparse,sys,json,itertools,hashlib
import numpy as np
import pandas as pd
from scipy import stats

p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--recount',type=Path,required=True);p.add_argument('--models',type=Path,required=True);p.add_argument('--qc',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
meta=pd.read_csv(a.data/'metadata/bkv_sample_metadata.csv').set_index('sample')
h=pd.read_csv(a.recount/'human_unique_counts.csv',index_col=0);v=pd.read_csv(a.recount/'viral_historical_counts.csv',index_col=0)
q=pd.read_csv(a.qc).set_index('sample')
sets={'all22':list(meta.index),'main_ge2000':list(h.columns[h.sum()>=2000]),'ge1000':list(h.columns[h.sum()>=1000]),'old_rule':list(h.columns[(h.sum()>=2000)&((h>=10).sum()>=30)])}
N=19999;seed=20260922

def bh(values):
 x=np.asarray(values,float);out=np.full(x.shape,np.nan);valid=np.isfinite(x);z=x[valid];ix=np.argsort(z);qv=np.minimum.accumulate((z[ix]*len(z)/np.arange(1,len(z)+1))[::-1])[::-1];back=np.empty(len(z));back[ix]=np.minimum(qv,1);out[valid]=back;return out

def perm(x,mask):
 # x samples x features; same group allocations reused across features.
 rng=np.random.default_rng(seed);n,k=len(mask),int(mask.sum());obs=x[mask].mean(0)-x[~mask].mean(0);hits=np.zeros(x.shape[1],int)
 for start in range(0,N,500):
  idx=np.argsort(rng.random((min(500,N-start),n)),axis=1)[:,:k]
  sums=x[idx].sum(1);d=sums/k-(x.sum(0)-sums)/(n-k)
  hits+=(np.abs(d)>=np.abs(obs)-1e-12).sum(0)
 return obs,(hits+1)/(N+1)

ab=pd.read_csv(a.data/'metaphlan_output/merged_abundance_table_species_phyloseq.csv',index_col=0)
tax=pd.read_csv(a.data/'metaphlan_output/merged_abundance_table_taxa_phyloseq.csv',index_col=0).loc[ab.index]
assert set(ab.columns)==set(meta.index) and tax.Kingdom.eq('k__Bacteria').all() and np.allclose(ab.sum(),100,atol=.001)
assert not ab.isna().any().any() and (ab>=0).all().all()
ab=ab.div(ab.sum())*100
ranks={};taxrows=[]
for rank in ['Class','Family','Genus','Species']:
 labels=tax[rank].fillna('Unclassified_'+rank).str.replace(r'^[a-z]__','',regex=True)
 ranks[rank]=ab.groupby(labels).sum()
 for name,ss in sets.items():
  x=ranks[rank][ss];x=x.loc[x.sum(axis=1)>0];mask=meta.loc[ss,'group'].eq('treatment').to_numpy()
  x.to_csv(a.out/f'{rank.lower()}_profiles_{name}.csv')
  desc=pd.DataFrame(index=x.index)
  for group,label in [(mask,'case'),(~mask,'control')]:
   z=x.iloc[:,group];desc[label+'_mean_pct']=z.mean(axis=1);desc[label+'_median_pct']=z.median(axis=1);desc[label+'_q1_pct']=z.quantile(.25,axis=1);desc[label+'_q3_pct']=z.quantile(.75,axis=1);desc[label+'_detected_n']=(z>0).sum(axis=1)
  desc['eligible_for_testing']=(x>0).sum(axis=1)>=3
  desc.to_csv(a.out/f'{rank.lower()}_descriptive_{name}.csv')
  for scale,eps in [('relative_pct',None),('CLR_0.001pct',.001),('CLR_0.0001pct',.0001),('CLR_0.01pct',.01)]:
   mat=x.to_numpy().T
   if eps is not None:
    mat=np.log(mat+eps);mat-=mat.mean(1,keepdims=True)
   # Presence in at least 3 samples limits unsupported sparse tests; full profiles exported.
   eligible=(x>0).sum(axis=1).to_numpy()>=3
   delta,pval=perm(mat[:,eligible],mask)
   for j,idx in enumerate(np.flatnonzero(eligible)):
    values=x.iloc[idx].to_numpy()
    taxrows.append(dict(strategy=name,rank=rank,scale=scale,taxon=x.index[idx],case_mean_pct=values[mask].mean(),control_mean_pct=values[~mask].mean(),case_median_pct=np.median(values[mask]),control_median_pct=np.median(values[~mask]),case_detected=int((values[mask]>0).sum()),control_detected=int((values[~mask]>0).sum()),effect_transformed=delta[j],pvalue=pval[j]))
t=pd.DataFrame(taxrows);t['FDR_all_ranks']=t.groupby(['strategy','scale']).pvalue.transform(bh);t['FDR_within_rank']=t.groupby(['strategy','scale','rank']).pvalue.transform(bh);t.to_csv(a.out/'bacterial_tests.csv',index=False)
viral=[]
for name,ss in sets.items():
 sf=pd.read_csv(a.models/name/'size_factors.csv').set_index('sample').size_factor.reindex(ss)
 mask=meta.loc[ss,'group'].eq('treatment').to_numpy()
 for scale,values in [('raw_log2',np.log2(v[ss]+1)),('human_sizefactor_log2',np.log2(v[ss].div(sf)+1)),('retained_read_CPM_log2',np.log2(v[ss].div(q.loc[ss,'retained_reads'])*1e6+1))]:
  delta,pval=perm(values.T.to_numpy(),mask)
  for i,mir in enumerate(v.index):viral.append(dict(strategy=name,scale=scale,miRNA=mir,delta=delta[i],pvalue=pval[i]))
vd=pd.DataFrame(viral);vd['FDR_two_viral']=vd.groupby(['strategy','scale']).pvalue.transform(bh);vd.to_csv(a.out/'viral_sensitivity.csv',index=False)
# Fixed historical candidates; each pair ranked and residualized on case/control.
candidates=['hsa-miR-345-5p','hsa-miR-361-5p','hsa-miR-193a-3p','hsa-miR-339-3p','hsa-miR-181c-3p','hsa-miR-96-5p','hsa-miR-151a-5p','hsa-miR-320b','hsa-miR-589-5p','hsa-let-7a-3p','bkv-miR-B1-3p','bkv-miR-B1-5p']
rows=[]
for name,ss in sets.items():
 sf=pd.read_csv(a.models/name/'size_factors.csv').set_index('sample').size_factor.reindex(ss)
 xx=pd.concat([h[ss],v[ss]]).loc[candidates].div(sf)
 yy=ranks['Class'][ss];yy=yy.loc[(yy>0).sum(axis=1)>=3]
 yy=np.log(yy+.001);yy=yy.sub(yy.mean(axis=0),axis=1)
 xr=stats.rankdata(xx.T.to_numpy(),axis=0);yr=stats.rankdata(yy.T.to_numpy(),axis=0)
 mask=meta.loc[ss,'group'].eq('treatment').to_numpy()
 for g in [mask,~mask]:xr[g]-=xr[g].mean(0);yr[g]-=yr[g].mean(0)
 xd=np.sqrt((xr*xr).sum(0));yd=np.sqrt((yr*yr).sum(0));den=xd[:,None]*yd[None,:]
 numer=np.einsum('ni,nj->ij',xr,yr)
 observed=np.divide(numer,den,out=np.full(den.shape,np.nan),where=den>0)
 hits=np.zeros(observed.shape,int);rng=np.random.default_rng(seed)
 for start in range(0,N,500):
  b=min(500,N-start);idx=np.tile(np.arange(len(ss)),(b,1))
  for g in [np.flatnonzero(mask),np.flatnonzero(~mask)]:idx[:,g]=g[np.argsort(rng.random((b,len(g))),axis=1)]
  numer=np.einsum('ni,bnj->bij',xr,yr[idx]);rp=np.divide(numer,den,out=np.zeros_like(numer),where=den>0)
  hits+=(np.abs(rp)>=np.abs(observed)-1e-12).sum(0)
 pv=(hits+1)/(N+1);pv[~np.isfinite(observed)]=np.nan
 for i,mir in enumerate(candidates):
  for j,cl in enumerate(yy.index):rows.append(dict(strategy=name,miRNA=mir,Class=cl,group_adjusted_rank_r=observed[i,j],pvalue=pv[i,j]))
c=pd.DataFrame(rows);c['FDR_all_pairs']=c.groupby('strategy').pvalue.transform(bh);c.to_csv(a.out/'group_conditioned_correlations.csv',index=False)
summary={'permutations':N,'seed':seed,'sets':sets,'bacterial_FDR05':t.groupby(['strategy','scale']).FDR_all_ranks.apply(lambda x:int((x<.05).sum())).to_dict(),'correlation_FDR05':c.groupby('strategy').FDR_all_pairs.apply(lambda x:int((x<.05).sum())).to_dict()}
summary['bacterial_FDR05']={str(k):val for k,val in summary['bacterial_FDR05'].items()}
(a.out/'table_analysis_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
