from pathlib import Path
import argparse,json,shutil
import pandas as pd
p=argparse.ArgumentParser();p.add_argument('--profiles',type=Path,required=True);p.add_argument('--mapping',type=Path,required=True);p.add_argument('--metadata',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
m=pd.read_csv(a.mapping,sep='\t',dtype=str);meta=pd.read_csv(a.metadata);m['sample']='S'+m.sample_title
assert m['sample'].is_unique and m.run_accession.is_unique and set(m['sample'])==set(meta['sample'])
assert (m.sample_alias=='R'+m.sample_title).all()
tax={};data={}
for row in m.itertuples():
 f=a.profiles/(row.run_accession+'_metaphlan4_profile_output.txt');values={}
 for line in f.read_text().splitlines():
  if line.startswith('#'):continue
  fields=line.split('\t');ranks=fields[0].split('|')
  if not ranks[-1].startswith('s__'):continue
  assert len(ranks)==7 and ranks[0]=='k__Bacteria'
  name=ranks[-1];assert name not in values
  assert name not in tax or tax[name]==ranks
  values[name]=float(fields[2]);tax[name]=ranks
 data[row.sample]=pd.Series(values)
x=pd.DataFrame(data).fillna(0).reindex(columns=meta['sample']);assert ((x.sum()-100).abs()<.001).all()
labels=sorted(tax);ids={label:'T'+str(i+1) for i,label in enumerate(labels)}
x=x.reindex(labels);x.index=[ids[k] for k in labels];x.index.name='TAXA'
t=pd.DataFrame([tax[k] for k in labels],index=x.index,columns=['Kingdom','Phylum','Class','Order','Family','Genus','Species'])
(a.out/'metadata').mkdir(parents=True,exist_ok=True);(a.out/'metaphlan_output').mkdir(exist_ok=True)
x.to_csv(a.out/'metaphlan_output/merged_abundance_table_species_phyloseq.csv');t.to_csv(a.out/'metaphlan_output/merged_abundance_table_taxa_phyloseq.csv')
meta.to_csv(a.out/'metadata/bkv_sample_metadata.csv',index=False);m.merge(meta,on='sample').to_csv(a.out/'metadata/verified_accession_sample_mapping.csv',index=False)
print('Rebuilt',len(x),'species across',x.shape[1],'samples')
