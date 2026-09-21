#!/usr/bin/env python3
"""Audit and extend the BK urine analysis from existing count/abundance tables.

Does not process raw reads, estimate absolute microbial loads, or alter inputs.
Run: python reanalyze.py --data-dir /path/to/BK_virus/data --out-dir ../
"""
import argparse
import hashlib
import itertools
import json
import platform
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import scipy
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def bh(p):
    p = np.asarray(p, float)
    out = np.full(p.shape, np.nan)
    keep = np.isfinite(p)
    pv = p[keep]
    if len(pv):
        order = np.argsort(pv, kind='stable')
        q = np.minimum.accumulate((pv[order]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1]
        z = np.empty_like(q); z[order] = np.minimum(q,1); out[keep] = z
    return out


def quantile_normalize(df, tie_method='legacy_min'):
    target = np.sort(df.to_numpy(), axis=0).mean(axis=1)
    lo = df.rank(method='min').astype(int)-1
    if tie_method == 'legacy_min':
        return pd.DataFrame(target[lo.to_numpy()], index=df.index, columns=df.columns)
    hi = df.rank(method='max').astype(int)-1
    cumulative = np.r_[0., np.cumsum(target)]
    vals = (cumulative[hi.to_numpy()+1]-cumulative[lo.to_numpy()])/(hi-lo+1).to_numpy()
    return pd.DataFrame(vals,index=df.index,columns=df.columns)


def compare(x, treatment, control):
    a=x[treatment].to_numpy(); b=x[control].to_numpy()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        p=stats.ttest_ind(a,b,axis=1,equal_var=False).pvalue
    delta=a.mean(axis=1)-b.mean(axis=1)
    return pd.DataFrame({'miRNA':x.index,'delta_log2_scale':delta,
                         'pvalue':np.nan_to_num(p,nan=1.)})


def correlations(x,y):
    rows=[]
    for mir in x.index:
        for tax in y.index:
            r,p=stats.spearmanr(x.loc[mir],y.loc[tax])
            rows.append((mir,tax,r,p))
    d=pd.DataFrame(rows,columns=['miRNA','Class','rho','pvalue'])
    d['fdr_168_or_actual_tests']=bh(d.pvalue)
    return d


def exact_mean_permutation(x, treatment_mask):
    """Two-sided absolute mean difference over every n-choose-k allocation.
    Observational groups require exchangeability; this is exploratory only.
    """
    a=x.to_numpy().T; n=len(a); k=int(treatment_mask.sum())
    obs=a[treatment_mask].mean(0)-a[~treatment_mask].mean(0)
    allocations=np.array(list(itertools.combinations(range(n), k)),dtype=int)
    extreme=np.zeros(a.shape[1],int)
    for start in range(0,len(allocations),512):
        sums=a[allocations[start:start+512]].sum(1)
        perm=sums/k-(a.sum(0)-sums)/(n-k)
        extreme+=(np.abs(perm)>=np.abs(obs)-1e-10).sum(0)
    return extreme/len(allocations)


def summarize(x,meta,do_test=False):
    c=meta.loc[meta.group=='control','sample'].tolist()
    t=meta.loc[meta.group=='treatment','sample'].tolist()
    out=pd.DataFrame(index=x.index)
    for name,ss in [('control',c),('case',t)]:
        out[name+'_n']=len(ss)
        out[name+'_mean_pct']=x[ss].mean(axis=1)
        out[name+'_median_pct']=x[ss].median(axis=1)
        out[name+'_q1_pct']=x[ss].quantile(.25,axis=1)
        out[name+'_q3_pct']=x[ss].quantile(.75,axis=1)
        out[name+'_detected_n']=(x[ss]>0).sum(axis=1)
    out['case_minus_control_pp']=out.case_mean_pct-out.control_mean_pct
    if do_test:
        out['permutation_p']=exact_mean_permutation(x[meta['sample']],meta.group.eq('treatment').to_numpy())
        out['fdr_within_rank']=bh(out.permutation_p)
    return out.sort_values('case_mean_pct',ascending=False)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--data-dir',type=Path,required=True)
    ap.add_argument('--out-dir',type=Path,required=True); args=ap.parse_args()
    data=args.data_dir.resolve(); root=args.out_dir.resolve()
    if root==data or data in root.parents:
        raise ValueError('Output must be outside the source data directory')
    out=root/'results'; figdir=root/'figures'; out.mkdir(parents=True,exist_ok=True);figdir.mkdir(exist_ok=True)
    manifest=[]
    def read(rel,index=False):
        p=data/rel; manifest.append({'relative_path':rel,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
        d=pd.read_csv(p,index_col=0 if index else None)
        if index and d.index.has_duplicates: raise ValueError(f'Duplicate feature IDs: {rel}')
        return d
    meta=read('metadata/bkv_sample_metadata.csv')
    if meta['sample'].duplicated().any() or set(meta.group)!={'treatment','control'}:
        raise ValueError('Invalid sample metadata')
    h=read('mirdeep/counts_matrix.csv',True)[meta['sample']]
    v=read('mirdeep/bkv_counts_matrix.csv',True)[meta['sample']]
    for matrix in [h,v]:
        if matrix.isna().any().any() or (matrix<0).any().any():raise ValueError('Missing/negative counts')
    qc=meta.set_index('sample').copy()
    qc['human_total']=h.sum();qc['human_features_ge10']=(h>=10).sum()
    qc['retained']=(qc.human_total>=2000)&(qc.human_features_ge10>=30)
    qc['viral_total_context_only']=v.sum()
    qc.to_csv(out/'sample_qc.csv')
    selected=meta.loc[meta['sample'].map(qc.retained)].copy()
    given=read('metadata/bkv_sample_metadata_16_qc.csv')
    assert set(selected['sample'])==set(given['sample']), 'QC selection differs from supplied subset'
    selected=given.copy(); samples=selected['sample'].tolist()
    t=selected.loc[selected.group=='treatment','sample'].tolist()
    c=selected.loc[selected.group=='control','sample'].tolist()
    hf=h[samples].loc[(h[samples]>=10).sum(axis=1)>=3]
    qn=quantile_normalize(np.log2(hf+1))
    vh=np.log2(v[samples]+1)
    hr=compare(qn,t,c); hr['feature_type']='human';hr['scale']='human_log2_quantile_legacy_min'
    vr=compare(vh,t,c); vr['feature_type']='viral';vr['scale']='viral_log2_raw_count_plus1'
    res=pd.concat([hr,vr],ignore_index=True)
    res['fdr_global_268']=bh(res.pvalue)
    res['nominal_candidate']=(res.pvalue<=.05)&(res.delta_log2_scale.abs()>=1)
    res.to_csv(out/'miRNA_all_tests_reproduced.csv',index=False)
    candidates=res[res.nominal_candidate].sort_values(['feature_type','pvalue']).copy()
    candidates.to_csv(out/'miRNA_candidates.csv',index=False)
    old=read('mirdeep/volcano_results_hybrid_human_qn_bkv_targeted_16_qc.csv').set_index('miRNA')
    check=res.set_index('miRNA').loc[old.index]
    assert np.allclose(check.delta_log2_scale,old.log2FC,atol=1e-9)
    assert np.allclose(check.pvalue,old.pvalue,atol=1e-9)
    # Sensitivities retain the same feature filter and patients.
    fixed=compare(quantile_normalize(np.log2(hf+1),'average_ties'),t,c)
    fixed['fdr_human']=bh(fixed.pvalue);fixed.to_csv(out/'human_average_ties_sensitivity.csv',index=False)
    rpm=np.log2(v[samples].div(h[samples].sum(),axis=1)*1e6+1)
    vs=compare(rpm,t,c);vs['fdr_two_viral_tests']=bh(vs.pvalue)
    vs.to_csv(out/'viral_per_million_human_miRNA_sensitivity.csv',index=False)
    # This denominator is human miRNA counts, NOT total sequencing reads or mL.
    depth=pd.DataFrame({'sample':samples,'group':selected.group.to_numpy(),'human_miRNA_reads':h[samples].sum().to_numpy()})
    depth.to_csv(out/'human_miRNA_depth.csv',index=False)
    abund=read('metaphlan_output/merged_abundance_table_species_phyloseq.csv',True)
    taxonomy=read('metaphlan_output/merged_abundance_table_taxa_phyloseq.csv',True)
    assert set(abund.index)==set(taxonomy.index),'Unmapped taxa'
    taxonomy=taxonomy.loc[abund.index]
    assert taxonomy.Kingdom.eq('k__Bacteria').all(),'Update denominator handling for nonbacterial taxa'
    assert np.allclose(abund.sum(),100,atol=.001)
    assert not abund.isna().any().any() and (abund>=0).all().all()
    abund=abund.div(abund.sum(),axis=1)*100
    rank_data={}; tested=[]
    for rank in ['Class','Family','Genus','Species']:
        labels=taxonomy[rank].fillna('Unclassified_'+rank).str.replace(r'^[a-z]__','',regex=True)
        x=abund.groupby(labels).sum(); xs=x[samples].loc[(x[samples]>0).any(axis=1)]
        rank_data[rank]=xs
        xs.to_csv(out/f'{rank.lower()}_per_sample_16.csv')
        table=summarize(xs,selected,True);table.index.name=rank
        table.to_csv(out/f'{rank.lower()}_summary_16.csv')
        z=table.reset_index().rename(columns={rank:'taxon'});z['rank']=rank;tested.append(z)
        alltable=summarize(x,meta,False);alltable.index.name=rank
        alltable.to_csv(out/f'{rank.lower()}_summary_all22.csv')
    multi=pd.concat(tested,ignore_index=True);multi['fdr_all_ranks']=bh(multi.permutation_p)
    multi.to_csv(out/'taxa_tests_all_ranks.csv',index=False)
    cls=rank_data['Class'].loc[rank_data['Class'].mean(axis=1).sort_values(ascending=False).index]
    names=candidates.miRNA.tolist()
    raw=pd.concat([np.log2(h[samples]+1),vh]).loc[names]
    corr=correlations(raw,cls)
    corr.to_csv(out/'correlations_original_scale.csv',index=False)
    normalized=pd.concat([qn,rpm]).loc[names]
    cs=correlations(normalized,cls);cs.to_csv(out/'correlations_normalized_sensitivity.csv',index=False)
    oldrho=read('mirdeep/correlation_miRNA_microbiome_all_classes_16_qc_9v7_rho.csv',True)
    rho=corr.pivot(index='miRNA',columns='Class',values='rho').loc[oldrho.index,oldrho.columns]
    assert np.allclose(rho,oldrho,atol=1e-9)
    # One group-aware descriptive sensitivity: correlations within each group.
    for name,ss in [('control',c),('case',t)]:
        rows=[]
        for mir in names:
            for tax in cls.index:
                a=normalized.loc[mir,ss];b=cls.loc[tax,ss]
                r=stats.spearmanr(a,b).statistic if a.nunique()>1 and b.nunique()>1 else np.nan
                rows.append((mir,tax,len(ss),r))
        pd.DataFrame(rows,columns=['miRNA','Class','n','rho_descriptive']).to_csv(out/f'correlations_within_{name}.csv',index=False)
    # Original paper candidate overlap, transcribed from its abstract, not raw original results.
    original_up=['hsa-miR-16-5p','hsa-miR-200c-3p','bkv-miR-B1-3p','hsa-let-7b-3p','hsa-miR-1269b','bkv-miR-B1-5p','hsa-miR-193a-3p','hsa-miR-944']
    original_down=['hsa-miR-134-5p','hsa-miR-4724-5p','hsa-miR-127-3p','hsa-miR-6500-3p','hsa-miR-507','hsa-miR-378b','hsa-miR-3911','hsa-miR-211-5p','hsa-miR-486-5p','hsa-miR-143-3p','hsa-miR-3195','hsa-miR-1307-5p','hsa-miR-29a-5p','hsa-miR-378f','hsa-miR-12136','hsa-miR-378g','hsa-miR-144-3p','hsa-miR-378a-3p','hsa-let-7i-5p','hsa-miR-204-5p','hsa-miR-146a-5p']
    orig=pd.DataFrame({'miRNA':original_up+original_down,'original_direction':['up']*len(original_up)+['down']*len(original_down)})
    orig.merge(res,on='miRNA',how='left').to_csv(out/'original_paper_candidate_comparison.csv',index=False)
    # Figures use explicit scales; the Sankey separates its two measurement domains.
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(12,5),gridspec_kw={'width_ratios':[1.3,1]})
    ss=candidates[candidates.feature_type=='human'].sort_values('delta_log2_scale')
    axes[0].barh(ss.miRNA,ss.delta_log2_scale,color=['#247b79' if x<0 else '#a54247' for x in ss.delta_log2_scale])
    axes[0].axvline(0,color='black',lw=.7);axes[0].set_xlabel('Difference in mean quantile-normalized log2(count + 1)')
    axes[0].set_title('A   Human miRNA candidates',loc='left',fontweight='bold')
    for i,(_,r) in enumerate(ss.iterrows()):axes[0].text(3.55,i,f'p={r.pvalue:.3f}',va='center',fontsize=9)
    axes[0].set_xlim(-3.1,4.9)
    for j,mir in enumerate(vh.index):
        for offset,ids,col in [(-.15,c,'#606060'),(.15,t,'#2965a0')]:
            axes[1].scatter(j+offset+np.linspace(-.05,.05,len(ids)),vh.loc[mir,ids],color=col,s=25,label=('Control' if offset<0 else 'BKPyV DNAemia') if j==0 else None)
    axes[1].set_xticks(range(len(vh)),[x.replace('bkv-miR-','') for x in vh.index]);axes[1].set_ylabel('log2(raw viral miRNA count + 1)')
    axes[1].set_title('B   Viral miRNA counts',loc='left',fontweight='bold');axes[1].legend(frameon=False)
    fig.tight_layout();fig.savefig(figdir/'Figure_1_miRNA.png',dpi=250);fig.savefig(figdir/'Figure_1_miRNA.svg');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(13,6),gridspec_kw={'width_ratios':[2.5,1]})
    palette=plt.get_cmap('tab20')(np.arange(len(cls)))
    order=c+t;bot=np.zeros(len(order))
    for i,cl in enumerate(cls.index):axes[0].bar(range(len(order)),cls.loc[cl,order],bottom=bot,color=palette[i],label=cl);bot+=cls.loc[cl,order].to_numpy()
    axes[0].set_xticks(range(len(order)),['C:'+s if s in c else 'T:'+s for s in order],rotation=60,ha='right');axes[0].set_ylabel('Bacterial relative abundance (%)');axes[0].set_ylim(0,100)
    axes[0].set_title('A   Individual samples',loc='left',fontweight='bold')
    means=pd.DataFrame({'Control':cls[c].mean(axis=1),'BKPyV DNAemia':cls[t].mean(axis=1)})
    bot=np.zeros(2)
    for i,cl in enumerate(cls.index):axes[1].bar([0,1],means.loc[cl],bottom=bot,color=palette[i]);bot+=means.loc[cl].to_numpy()
    axes[1].set_xticks([0,1],['Control\nn=7','DNAemia\nn=9']);axes[1].set_ylim(0,100);axes[1].set_title('B   Mean of sample percentages',loc='left',fontweight='bold')
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=4,frameon=False,fontsize=9)
    fig.tight_layout(rect=[0,.19,1,1]);fig.savefig(figdir/'Figure_2_bacteria.png',dpi=250);fig.savefig(figdir/'Figure_2_bacteria.svg');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(15,6),sharey=True)
    for ax,d,title in zip(axes,[corr,cs],['A   Original raw-count ranks','B   Normalization sensitivity']):
        matrix=d.pivot(index='miRNA',columns='Class',values='rho').loc[names,cls.index]
        im=ax.imshow(matrix,vmin=-1,vmax=1,cmap='RdBu_r',aspect='auto')
        ax.set_xticks(range(len(cls)),cls.index,rotation=60,ha='right',fontsize=8)
        ax.set_yticks(range(len(names)),[m.replace('bkv-miR-B1-3p','bkv-miR-B1-3p*') for m in names],fontsize=9)
        ax.set_title(title,loc='left',fontweight='bold')
    fig.tight_layout(rect=[0,0,.94,1]);cax=fig.add_axes([.95,.32,.015,.5]);fig.colorbar(im,cax=cax,label='Spearman rho')
    fig.savefig(figdir/'Supplement_Correlations.png',dpi=250);fig.savefig(figdir/'Supplement_Correlations.svg');plt.close(fig)
    from sankey_figure import draw_sankey
    draw_sankey(cls, raw, c, t, figdir, out)
    for rank in ['Family','Genus','Species']:
        x=rank_data[rank];top=x.mean(axis=1).nlargest(15).index
        fig,ax=plt.subplots(figsize=(9,6));y=np.arange(len(top))
        ax.barh(y-.18,x.loc[top,c].mean(axis=1),height=.36,label='Control',color='#777777')
        ax.barh(y+.18,x.loc[top,t].mean(axis=1),height=.36,label='BKPyV DNAemia',color='#2965a0')
        ax.set_yticks(y,[z.replace('_',' ') for z in top]);ax.invert_yaxis();ax.set_xlabel('Mean bacterial relative abundance (%)');ax.legend(frameon=False)
        ax.set_title(f'Top 15 {rank.lower()} taxa by overall mean abundance')
        fig.tight_layout();fig.savefig(figdir/f'Supplement_{rank}.png',dpi=250);fig.savefig(figdir/f'Supplement_{rank}.svg');plt.close(fig)
    summary={'n':len(samples),'case_n':len(t),'control_n':len(c),'human_tests':len(hr),'total_tests':len(res),
             'nominal_candidates':len(candidates),'global_fdr_significant':int((res.fdr_global_268<.05).sum()),
             'min_global_fdr':float(res.fdr_global_268.min()),'correlation_min_fdr':float(corr.fdr_168_or_actual_tests.min()),
             'normalized_correlation_min_fdr':float(cs.fdr_168_or_actual_tests.min()),
             'taxa_test_min_global_fdr':float(multi.fdr_all_ranks.min()),
             'viral_human_normalized':vs.to_dict('records'),
             'depth_group_median':depth.groupby('group').human_miRNA_reads.median().to_dict(),
             'reproduction_checks':'QC sample set, all 268 effect sizes/p values, and 168 correlation coefficients matched supplied exports',
             'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__},
             'inputs':manifest}
    (out/'run_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps({k:v for k,v in summary.items() if k!='inputs'},indent=2))


if __name__=='__main__':main()
