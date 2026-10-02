"""Independent consistency checks on saved numerical outputs; no expected discoveries."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control


def check_bh(frame, p_column, q_column):
    mask = frame[p_column].notna()
    assert frame.loc[~mask, q_column].isna().all()
    p = frame.loc[mask, p_column].to_numpy()
    assert ((p >= 0) & (p <= 1)).all()
    np.testing.assert_allclose(frame.loc[mask, q_column], false_discovery_control(p), atol=1e-12)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--recount', type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    counts = pd.read_csv(args.recount / 'human_unique_counts.csv', index_col=0)
    meta = pd.read_csv(root / 'verified_data/metadata/bkv_sample_metadata.csv').set_index('sample')
    assert meta.index.is_unique and counts.index.is_unique and counts.columns.is_unique
    assert set(meta.index) == set(counts.columns)
    assert np.isfinite(counts.to_numpy()).all() and (counts >= 0).all().all()
    np.testing.assert_array_equal(counts, np.floor(counts))
    summary = pd.read_csv(root / 'results/count_models/model_summary.csv')
    for row in summary.itertuples():
        folder = root / 'results/count_models' / row.strategy
        sf = pd.read_csv(folder / 'size_factors.csv').set_index('sample')
        assert (sf.size_factor > 0).all() and np.isfinite(sf.size_factor).all()
        np.testing.assert_array_equal(sf.raw_total, counts[sf.index].sum())
        selected = counts[sf.index]
        keep = selected.index[(selected >= 10).sum(axis=1) >= 3]
        deseq = pd.read_csv(folder / 'deseq2_results.csv').set_index('miRNA')
        edger = pd.read_csv(folder / 'edger_TMMwsp_results.csv').set_index('miRNA')
        assert set(deseq.index) == set(edger.index) == set(keep)
        check_bh(deseq, 'pvalue', 'padj')
        check_bh(edger, 'PValue', 'FDR')
        normalized = pd.read_csv(folder / 'deseq2_normalized_counts.csv').set_index('miRNA')
        np.testing.assert_allclose(normalized, counts.loc[normalized.index, normalized.columns].div(sf.size_factor), rtol=1e-12)
        np.testing.assert_allclose(deseq.wald_lower95, deseq.log2FoldChange-1.96*deseq.lfcSE)
        assert row.n == len(sf) and row.cases == meta.loc[sf.index].group.eq('treatment').sum()
        assert row.deseq2_FDR05 == (deseq.padj < .05).sum()
        assert row.edger_FDR05 == (edger.FDR < .05).sum()
    tabs = root / 'results/tables'
    bacterial = pd.read_csv(tabs / 'bacterial_tests.csv')
    for (_, _), group in bacterial.groupby(['strategy', 'scale']):
        check_bh(group, 'pvalue', 'FDR_all_ranks')
        for _, rank in group.groupby('rank'):
            check_bh(rank, 'pvalue', 'FDR_within_rank')
    for filename, keys, q in [
        ('group_conditioned_correlations.csv', ['strategy'], 'FDR_all_pairs'),
        ('viral_sensitivity.csv', ['strategy', 'scale'], 'FDR_two_viral')]:
        frame = pd.read_csv(tabs / filename)
        if filename == 'group_conditioned_correlations.csv':
            finite_r = frame.group_adjusted_rank_r.dropna()
            assert np.isfinite(finite_r).all() and (finite_r.abs() <= 1 + 1e-12).all()
        for _, group in frame.groupby(keys):
            check_bh(group, 'pvalue', q)
        finite = frame.pvalue.dropna()
        np.testing.assert_allclose(finite * 20000, np.round(finite * 20000), atol=1e-9)
    for path in tabs.glob('*_profiles_*.csv'):
        frame = pd.read_csv(path, index_col=0)
        assert (frame >= 0).all().all() and np.isfinite(frame.to_numpy()).all()
        np.testing.assert_allclose(frame.sum(), 100, atol=1e-8)
        stem = path.name.replace('_profiles_', '_descriptive_')
        desc = pd.read_csv(tabs / stem, index_col=0)
        for group, label in [('treatment', 'case'), ('control', 'control')]:
            samples = [s for s in frame.columns if meta.loc[s, 'group'] == group]
            np.testing.assert_allclose(desc[label+'_mean_pct'], frame[samples].mean(axis=1), atol=1e-12)
    print('PASS: counts, model normalization, feature filters, effect intervals, BH families, permutation grid, bacterial closure and group summaries')


if __name__ == '__main__':
    main()
