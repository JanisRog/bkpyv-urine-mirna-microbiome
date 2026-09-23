from pathlib import Path
import argparse,re,json,hashlib
import pandas as pd
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--profiles',type=Path,required=True);p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
old=pd.read_csv(a.data/'metaphlan_output/merged_abundance_table_species_phyloseq.csv',index_col=0)
tax=pd.read_csv(a.data/'metaphlan_output/merged_abundance_table_taxa_phyloseq.csv',index_col=0).loc[old.index]
old.index=tax.Species.str.replace('s__','',regex=False);old=old.groupby(level=0).sum()
new={};qc=[];hashes={}
for f in sorted(a.profiles.glob('*_metaphlan4_profile_output.txt')):
 acc=f.name.split('_')[0];text=f.read_text();lines=text.splitlines();species={}
 for line in lines:
  if line.startswith('#'):continue
  fields=line.split('\t');label=fields[0].split('|')[-1]
  if label.startswith('s__'):species[label[3:]]=float(fields[2])
 new[acc]=pd.Series(species,dtype=float)
 aln=a.profiles/(acc+'_con_alnstats.txt');at=aln.read_text()
 def count(label):
  match=re.search(r'^(\d+) \+ (\d+) '+re.escape(label)+r'(?:\s|$)',at,re.M)
  if not match:raise ValueError((acc,label))
  return int(match[1])+int(match[2])
 primary=count('primary');mapped=count('primary mapped');both=count('with itself and mate mapped');singles=count('singletons')
 processed=int(re.search(r'#(\d+) reads processed',text)[1])
 # For paired primary alignments, one mapped singleton implies one unmapped mate.
 expected_both_unmapped=primary-both-2*singles
 qc.append(dict(accession=acc,primary_read_records=primary,primary_host_mapped=mapped,host_mapped_pct=100*mapped/primary,unmapped_primary=primary-mapped,mapped_singletons=singles,expected_both_unmapped_read_records=expected_both_unmapped,metaphlan_reported_reads=processed,processed_matches_both_unmapped=processed==expected_both_unmapped,profile_species=len(species),species_sum=sum(species.values()),database=lines[0][1:]))
 for file in [f,aln]:hashes[file.name]=hashlib.sha256(file.read_bytes()).hexdigest()
new=pd.DataFrame(new).fillna(0);new.to_csv(a.out/'downloaded_species_by_accession.csv');pd.DataFrame(qc).to_csv(a.out/'dna_read_qc_by_accession.csv',index=False)
# Similarity is a diagnostic, not an accession-to-patient crosswalk.
features=old.index.union(new.index);o=old.reindex(features).fillna(0);n=new.reindex(features).fillna(0)
d=np.abs(n.values.T[:,None,:]-o.values.T[None,:,:]).sum(axis=2)
pd.DataFrame(d,index=new.columns,columns=old.columns).to_csv(a.out/'profile_L1_distance_to_historical.csv')
summary={'profile_count':len(new.columns),'historical_species':len(old),'downloaded_species':len(new),'overlap_species_names':len(old.index.intersection(new.index)),'closest_L1_range':[float(d.min(axis=1).min()),float(d.min(axis=1).max())],'exact_profile_matches_at_0_001_L1':int((d.min(axis=1)<.001).sum()),'note':'No clinical sample identities inferred from similarity.'}
(a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(a.out/'input_checksums.json').write_text(json.dumps(hashes,indent=2)+'\n');print(json.dumps(summary,indent=2));print(pd.DataFrame(qc).to_string(index=False))
