#!/usr/bin/env python3
"""Exploratory group-conditioned human-miRNA × bacterial-genus heatmap."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
from scipy import stats


N_PERMUTATIONS = 19_999
SEED = 20260922
PSEUDOCOUNT_PCT = 0.001


def bh(pvalues):
    values = np.asarray(pvalues,dtype=float)
    out = np.full(values.shape,np.nan)
    good = np.isfinite(values)
    x = values[good]
    order = np.argsort(x)
    adjusted = np.minimum.accumulate(
        (x[order]*len(x)/np.arange(1,len(x)+1))[::-1]
    )[::-1]
    back = np.empty(len(x))
    back[order] = np.minimum(adjusted,1)
    out[good] = back
    return out


def association_matrix(human, genera, samples, group, sf, mirnas, taxon_names,
                       pseudocount=PSEUDOCOUNT_PCT, normalized=True):
    # Spearman-style ranks on host-normalized miRNA counts and genus CLR values.
    # The CLR reference uses all genera, not just the displayed subset.
    mirna_values = human.loc[mirnas,samples]
    if normalized:
        mirna_values = mirna_values.div(sf.loc[samples],axis=1)
    genus_values = genera[samples].astype(float)
    log_genera = np.log(genus_values+pseudocount)
    clr = log_genera.sub(log_genera.mean(axis=0),axis=1).loc[taxon_names]
    xr = stats.rankdata(mirna_values.T.to_numpy(),axis=0).astype(float)
    yr = stats.rankdata(clr.T.to_numpy(),axis=0).astype(float)
    mask = group.loc[samples].eq("treatment").to_numpy()
    for selected in (mask,~mask):
        xr[selected] -= xr[selected].mean(axis=0)
        yr[selected] -= yr[selected].mean(axis=0)
    den = np.sqrt((xr*xr).sum(axis=0))[:,None]*np.sqrt((yr*yr).sum(axis=0))[None,:]
    observed = np.divide(np.einsum("ni,nj->ij",xr,yr),den,
                         out=np.full(den.shape,np.nan),where=den>0)

    rng = np.random.default_rng(SEED)
    hits = np.zeros(observed.shape,dtype=int)
    for start in range(0,N_PERMUTATIONS,500):
        batch = min(500,N_PERMUTATIONS-start)
        idx = np.tile(np.arange(len(samples)),(batch,1))
        for selected in (np.flatnonzero(mask),np.flatnonzero(~mask)):
            idx[:,selected] = selected[np.argsort(rng.random((batch,len(selected))),axis=1)]
        numer = np.einsum("ni,bnj->bij",xr,yr[idx])
        permuted = np.divide(numer,den,out=np.zeros_like(numer),where=den>0)
        hits += (np.abs(permuted)>=np.abs(observed)-1e-12).sum(axis=0)
    pvalue = (hits+1)/(N_PERMUTATIONS+1)
    pvalue[~np.isfinite(observed)] = np.nan
    qvalue = bh(pvalue.ravel()).reshape(pvalue.shape)
    return observed,pvalue,qvalue


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,required=True)
    parser.add_argument("--models",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error(f"Output exists: {args.out}")
    args.out.mkdir(parents=True)
    meta = pd.read_csv(args.root/"verified_data/metadata/bkv_sample_metadata.csv").set_index("sample")
    human = pd.read_csv(args.root/"inputs/recount/human_unique_counts.csv",index_col=0)
    viral = pd.read_csv(args.root/"inputs/recount/viral_historical_counts.csv",index_col=0)
    counts = pd.concat([human,viral])
    genera = pd.read_csv(args.root/"results/tables/genus_profiles_all22.csv",index_col=0)
    assert len(meta)==22 and set(meta.index)==set(human.columns)==set(genera.columns)
    de = pd.read_csv(args.models/"main_ge2000/deseq2_results.csv")
    selected = de[(de.feature_type.eq("human")) & (de.pvalue<.05) &
                  (de.log2FoldChange.abs()>=1)].sort_values("pvalue")
    human_mirnas = selected.miRNA.tolist()
    mirnas = human_mirnas + ["bkv-miR-B1-3p","bkv-miR-B1-5p"]
    assert len(human_mirnas)==15 and len(mirnas)==17 and len(set(mirnas))==17

    main_sf = pd.read_csv(args.root/"results/count_models/main_ge2000/size_factors.csv").set_index("sample")
    samples = main_sf.index.tolist()
    assert len(samples)==16 and set(meta.loc[samples,"group"])=={"control","treatment"}
    taxon_names = sorted(genera.index[(genera[samples]>0).sum(axis=1)>=3],
                         key=lambda name:(-(genera.loc[name,samples]>0).sum(),name))
    assert len(taxon_names)==13
    detection = (genera.loc[taxon_names,samples]>0).sum(axis=1)
    obs,pval,qval = association_matrix(counts,genera,samples,meta.group,
                                       main_sf.size_factor,mirnas,taxon_names)
    raw_obs,raw_p,raw_q = association_matrix(counts,genera,samples,meta.group,
                                             main_sf.size_factor,mirnas,taxon_names,
                                             normalized=False)
    rows = []
    for i,mirna in enumerate(mirnas):
        for j,genus in enumerate(taxon_names):
            rows.append({"subset":"main_9_cases_7_controls","miRNA":mirna,
                         "genus":genus,"genus_detected_n":int(detection[genus]),
                         "group_adjusted_rank_r":obs[i,j],"permutation_p":pval[i,j],
                         "BH_q_across_all_displayed_pairs":qval[i,j]})
    result = pd.DataFrame(rows)
    result.to_csv(args.out/"Figure_2_miRNA_genus_association_tests.csv",index=False)
    assert len(result)==221
    raw_rows = []
    for i,mirna in enumerate(mirnas):
        for j,genus in enumerate(taxon_names):
            raw_rows.append({"subset":"main_9_cases_7_controls","miRNA":mirna,
                             "genus":genus,"genus_detected_n":int(detection[genus]),
                             "group_adjusted_rank_r":raw_obs[i,j],
                             "permutation_p":raw_p[i,j],
                             "BH_q_across_all_displayed_pairs":raw_q[i,j]})
    pd.DataFrame(raw_rows).to_csv(args.out/"Figure_2_raw_count_sensitivity.csv",index=False)

    # Keep the 15 main-subset candidates fixed in the all-22 sensitivity.
    all_sf = pd.read_csv(args.root/"results/count_models/all22/size_factors.csv").set_index("sample")
    all_samples = all_sf.index.tolist()
    all_obs,all_p,all_q = association_matrix(counts,genera,all_samples,meta.group,
                                              all_sf.size_factor,mirnas,taxon_names)
    all_rows = []
    for i,mirna in enumerate(mirnas):
        for j,genus in enumerate(taxon_names):
            all_rows.append({"subset":"all_11_cases_11_controls","miRNA":mirna,
                             "genus":genus,"group_adjusted_rank_r":all_obs[i,j],
                             "permutation_p":all_p[i,j],"BH_q_across_all_displayed_pairs":all_q[i,j]})
    pd.DataFrame(all_rows).to_csv(args.out/"Supplement_all22_miRNA_genus_associations.csv",index=False)

    sensitivity = []
    for pseudo in (.0001,.001,.01):
        rr,pp,qq = association_matrix(counts,genera,samples,meta.group,
                                      main_sf.size_factor,mirnas,taxon_names,pseudo)
        sensitivity.append({"pseudocount_percentage_points":pseudo,
                            "finite_pairs":int(np.isfinite(pp).sum()),
                            "FDR_below_0_05":int(np.nansum(qq<.05)),
                            "minimum_q":float(np.nanmin(qq))})
    pd.DataFrame(sensitivity).to_csv(args.out/"Supplement_pseudocount_sensitivity.csv",index=False)

    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#e0e0e0")
    fig, ax = plt.subplots(figsize=(12.5,9.1))
    fig.subplots_adjust(left=.19,right=.88,top=.91,bottom=.27)
    norm = TwoSlopeNorm(vmin=-1,vcenter=0,vmax=1)
    image = ax.imshow(obs,cmap=cmap,norm=norm,interpolation="nearest",aspect="auto")
    ax.set_xticks(np.arange(len(taxon_names)),taxon_names,rotation=58,
                  ha="right",fontsize=10)
    display_mirnas = human_mirnas + ["bkv-miR-B1-3p*","bkv-miR-B1-5p"]
    ax.set_yticks(np.arange(len(mirnas)),display_mirnas,fontsize=10)
    ax.set_xticks(np.arange(-.5,len(taxon_names),1),minor=True)
    ax.set_yticks(np.arange(-.5,len(mirnas),1),minor=True)
    ax.grid(which="minor",color="white",linewidth=1.0)
    ax.axhline(14.5,color="#334155",linewidth=1.5)
    ax.tick_params(which="minor",bottom=False,left=False)
    ax.tick_params(axis="both",length=0,pad=5)
    ax.set_title("Figure 2  Exploratory urinary miRNA–genus associations",
                 fontsize=15,fontweight="bold",pad=16)
    for spine in ax.spines.values():
        spine.set_visible(False)
    bar = fig.add_axes([.9,.33,.018,.49])
    cb = fig.colorbar(image,cax=bar,ticks=[-1,-.5,0,.5,1])
    cb.set_label("Group-adjusted rank correlation",labelpad=10)
    n_supported = int(np.nansum(qval<.05))
    raw_supported = int(np.nansum(raw_q<.05))
    fig.text(.53,.035,f"9 cases + 7 controls  |  {N_PERMUTATIONS:,} within-group permutations  |  "
             f"FDR q < 0.05: {n_supported}/{len(result)}",
             ha="center",fontsize=10,color="#475569")
    fig.savefig(args.out/"Figure_2_miRNA_genus_heatmap.png",dpi=240)
    fig.savefig(args.out/"Figure_2_miRNA_genus_heatmap.svg")
    plt.close(fig)
    (args.out/"Figure_2_caption.txt").write_text(
        "Figure 2. Exploratory associations between 15 human miRNAs meeting the "
        "Figure 1 descriptive rule (unadjusted P<0.05 and |log2 fold change|≥1), "
        "two prespecified viral miRNAs, and 13 bacterial genera detected in at least "
        "three of the 16 paired samples (nine blood-BKPyV-DNAemia cases, seven controls). "
        "The horizontal rule separates human from viral features. No human miRNA met "
        "FDR<0.05 for differential abundance; the viral 3p sequence had q=0.057 and "
        "the viral 5p primary P and q were filtered by Cook's distance. "
        "All miRNA counts were divided by host-derived "
        "DESeq2 size factors before ranking. Each cell shows the correlation after "
        "removing case/control group means from both "
        "ranked variables. Genus percentages were centered-log-ratio transformed "
        "across all observed genera with a 0.001 "
        "percentage-point pseudocount. Two-sided P values use 19,999 label "
        "permutations restricted within groups; Benjamini–Hochberg q values cover "
        "all 221 displayed pairs. No pair had q<0.05. "
        "Raw-count, all-22-recipient and pseudocount checks are supplied separately. "
        "Selection of miRNAs from the same cohort makes these associations exploratory. "
        "The asterisk marks the BKPyV/JCPyV-shared 3p sequence; viral 3p counts "
        "retain the BK-only assignment and cannot establish species origin. "
        "These relative, low-biomass taxonomic "
        "profiles cannot establish bacterial viability or direct miRNA interactions.\n"
    )
    summary = {"main_samples":len(samples),"cases":int(meta.loc[samples,"group"].eq("treatment").sum()),
               "controls":int(meta.loc[samples,"group"].eq("control").sum()),
               "human_miRNAs":len(human_mirnas),"viral_miRNAs":2,
               "miRNAs":len(mirnas),"genera":len(taxon_names),"pairs":len(result),
               "finite_pairs":int(result.permutation_p.notna().sum()),
               "raw_FDR_below_0_05":raw_supported,
               "raw_min_q":float(np.nanmin(raw_q)),
               "main_FDR_below_0_05":n_supported,
               "main_min_q":float(result.BH_q_across_all_displayed_pairs.min()),
               "all22_FDR_below_0_05":int(np.nansum(all_q<.05)),
               "pseudocount_sensitivity":sensitivity}
    (args.out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    main()
