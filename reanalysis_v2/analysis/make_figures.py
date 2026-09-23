from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sankey_figure import draw_sankey

p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--recount',type=Path,required=True);p.add_argument('--metadata',type=Path,required=True);a=p.parse_args();root=a.root;fdir=root/'figures';fdir.mkdir(exist_ok=True);tables=root/'results/tables';models=root/'results/count_models'
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
meta=pd.read_csv(a.metadata).set_index('sample');sf=pd.read_csv(models/'main_ge2000/size_factors.csv').set_index('sample').size_factor;ss=sf.index.tolist();controls=[s for s in ss if meta.loc[s,'group']=='control'];cases=[s for s in ss if meta.loc[s,'group']=='treatment']
names=['hsa-miR-345-5p','hsa-miR-361-5p','hsa-miR-193a-3p','hsa-miR-339-3p','hsa-miR-181c-3p','hsa-miR-96-5p','hsa-miR-151a-5p','hsa-miR-320b','hsa-miR-589-5p','hsa-let-7a-3p']
r=pd.read_csv(models/'main_ge2000/deseq2_results.csv').set_index('miRNA').reindex(names)
fig,ax=plt.subplots(figsize=(9,5.8),layout='constrained');y=np.arange(len(names));ax.errorbar(r.log2FoldChange,y,xerr=1.96*r.lfcSE,fmt='o',color='#285b78',capsize=3);ax.axvline(0,color='gray',lw=1);ax.set(yticks=y,yticklabels=names,xlabel='Case versus control log2 fold change (95% Wald interval)',title='Historical candidates under revised counting and DESeq2');ax.invert_yaxis()
for i,mir in enumerate(names):
 if pd.isna(r.loc[mir,'log2FoldChange']):ax.text(.15,i,'Not tested: below feature-count filter',va='center',fontsize=9,color='#555555')
fig.savefig(fdir/'Figure_1_human_miRNA.png',dpi=200);fig.savefig(fdir/'Figure_1_human_miRNA.svg');plt.close(fig)
cls=pd.read_csv(tables/'class_profiles_all22.csv',index_col=0);order=meta.query("group=='control'").index.tolist()+meta.query("group=='treatment'").index.tolist();major=cls.mean(axis=1).nlargest(7).index;display=cls.loc[major].copy();display.loc['Other classes']=cls.drop(major).sum();colors=plt.get_cmap('tab10')(np.arange(len(display)))
fig,axs=plt.subplots(1,2,figsize=(13,5.5),gridspec_kw={'width_ratios':[3,1]},layout='constrained')
for ax,dat in [(axs[0],display[order]),(axs[1],pd.DataFrame({'Control':display[meta.query("group=='control'").index].mean(axis=1),'Case':display[meta.query("group=='treatment'").index].mean(axis=1)}))]:
 bottom=np.zeros(dat.shape[1])
 for (tax,row),color in zip(dat.iterrows(),colors):ax.bar(np.arange(dat.shape[1]),row,bottom=bottom,color=color,label=tax);bottom+=row.to_numpy()
 ax.set(ylim=(0,100),xticks=np.arange(dat.shape[1]),xticklabels=dat.columns,ylabel='Bacterial relative abundance (%)')
axs[0].tick_params(axis='x',rotation=90);axs[0].set_title('All 22 DNA profiles: controls then cases');axs[0].axvline(10.5,color='black',lw=.8);axs[1].set_title('Mean of individual profiles');fig.legend(handles=axs[0].get_legend_handles_labels()[0],labels=list(display.index),loc='outside lower center',ncol=4,frameon=False);fig.savefig(fdir/'Figure_2_bacteria.png',dpi=200);fig.savefig(fdir/'Figure_2_bacteria.svg');plt.close(fig)
h=pd.read_csv(a.recount/'human_unique_counts.csv',index_col=0);v=pd.read_csv(a.recount/'viral_historical_counts.csv',index_col=0)
normalized=np.log2(pd.concat([h,v])[ss].div(sf)+1).loc[names+['bkv-miR-B1-3p','bkv-miR-B1-5p']]
draw_sankey(cls[ss],normalized,controls,cases,fdir,tables)
fig,axs=plt.subplots(1,3,figsize=(13,4.6),layout='constrained')
for ax,(label,vals) in zip(axs,[('Raw counts',np.log2(v[ss]+1)),('Human size-factor scaled',np.log2(v[ss].div(sf)+1)),('DESeq2 size factors',None)]):
 if vals is None:
  for i,group in enumerate([controls,cases]):ax.scatter(np.full(len(group),i),sf[group],color=['#357599','#ad553d'][i])
  ax.set(xticks=[0,1],xticklabels=['Control','Case'],yscale='log',ylabel='Size factor')
 else:
  for j,mir in enumerate(v.index):
   for i,group in enumerate([controls,cases]):
    pos=j*3+i;ax.scatter(pos+np.linspace(-.12,.12,len(group)),vals.loc[mir,group],color=['#357599','#ad553d'][i],s=25)
  ax.set(xticks=[0,1,3,4],xticklabels=['3p Ctrl','3p Case','5p Ctrl','5p Case'],ylabel='log2(value + 1)')
 ax.set_title(label)
fig.savefig(fdir/'Supplement_viral_normalization.png',dpi=200);plt.close(fig)
# Directions across models and subsets: fixed historical human candidates.
strategies=['main_ge2000','all22','ge1000','old_rule'];mat=[]
for st in strategies:mat.append(pd.read_csv(models/st/'deseq2_results.csv').set_index('miRNA').reindex(names).log2FoldChange)
mat=pd.DataFrame(mat,index=strategies).T
fig,ax=plt.subplots(figsize=(8,5.5),layout='constrained');bound=np.nanmax(np.abs(mat.values));im=ax.imshow(mat,aspect='auto',cmap=plt.get_cmap('RdBu_r').with_extremes(bad='#dddddd'),vmin=-bound,vmax=bound);ax.set(yticks=np.arange(len(names)),yticklabels=names,xticks=np.arange(4),xticklabels=['≥2,000 (9/7)','All (11/11)','≥1,000 (9/8)','Old rule (8/7)'],title='DESeq2 effect estimates across inclusion strategies');fig.colorbar(im,ax=ax,label='Unshrunk log2 fold change');[ax.text(j,i,'Not tested',ha='center',va='center',fontsize=8,color='#555555') for i,j in zip(*np.where(mat.isna().to_numpy()))];fig.savefig(fdir/'Supplement_inclusion_sensitivity.png',dpi=200);plt.close(fig)
c=pd.read_csv(tables/'group_conditioned_correlations.csv').query("strategy=='main_ge2000'").pivot(index='miRNA',columns='Class',values='group_adjusted_rank_r').reindex(names+list(v.index));fig,ax=plt.subplots(figsize=(10,6),layout='constrained');im=ax.imshow(c,aspect='auto',vmin=-1,vmax=1,cmap='RdBu_r');ax.set(yticks=np.arange(len(c)),yticklabels=c.index,xticks=np.arange(c.shape[1]),xticklabels=c.columns,title='Group-conditioned rank associations; none survives FDR');ax.tick_params(axis='x',rotation=45);fig.colorbar(im,ax=ax,label='Residual-rank correlation');fig.savefig(fdir/'Supplement_correlations.png',dpi=200);plt.close(fig)
# Supplemental full-rank summaries, all22.
for rank in ['family','genus','species']:
 x=pd.read_csv(tables/f'{rank}_profiles_all22.csv',index_col=0);z=x.loc[x.mean(axis=1).nlargest(15).index];y=np.arange(len(z));fig,ax=plt.subplots(figsize=(9,6),layout='constrained');ax.barh(y-.2,z[meta.query("group=='control'").index].mean(axis=1),height=.4,label='Control');ax.barh(y+.2,z[meta.query("group=='treatment'").index].mean(axis=1),height=.4,label='Case');ax.set(yticks=y,yticklabels=z.index,xlabel='Mean bacterial relative abundance (%)',title=f'All 22 profiles: {rank}');ax.invert_yaxis();ax.legend();fig.savefig(fdir/f'Supplement_{rank}.png',dpi=200);plt.close(fig)
print('Figures complete')
