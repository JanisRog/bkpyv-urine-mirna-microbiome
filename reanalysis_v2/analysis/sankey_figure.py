"""Two-scale descriptive Sankey; link widths retain original numerical definitions."""
import xml.etree.ElementTree as ET
ET.register_namespace('', 'http://www.w3.org/2000/svg')
ET.register_namespace('xlink', 'http://www.w3.org/1999/xlink')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle


def draw_sankey(cls, raw, controls, cases, figdir, out):
    order = sorted(controls, key=lambda s: int(s[1:])) + sorted(cases, key=lambda s: int(s[1:]))
    # Display-only grouping: preserve every sample's total and full analysis tables.
    original_cls = cls[order].copy()
    low = original_cls.max(axis=1) < 1.0
    mapping = pd.DataFrame({
        'Class': original_cls.index,
        'max_sample_pct': original_cls.max(axis=1).to_numpy(),
        'display_class': ['Other bacterial classes' if flag else name
                          for name, flag in zip(original_cls.index, low)],
    })
    mapping.to_csv(out / 'sankey_class_display_mapping.csv', index=False)
    cls = original_cls.loc[~low].copy()
    if low.any():
        cls.loc['Other bacterial classes'] = original_cls.loc[low].sum(axis=0)
    assert np.allclose(cls.sum(axis=0), original_cls.sum(axis=0))
    # DESeq2 human-size-factor normalized log2-count deviations; descriptive only.
    delta = raw.sub(raw[controls].mean(axis=1), axis=0)[order]
    weight = delta.abs().where(delta.abs() >= .5, 0.)
    delta.to_csv(out / 'sankey_signed_log2_deviations.csv')
    fig, ax = plt.subplots(figsize=(18, 12))
    fig.subplots_adjust(left=.015, right=.985, top=.875, bottom=.15)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
    original_palette = dict(zip(original_cls.index, plt.get_cmap('tab20')(np.arange(len(original_cls)))))
    palette = [original_palette[name] if name in original_palette else '#8a8a8a' for name in cls.index]
    sample_y = np.linspace(.95, .05, len(order))
    # A single width conversion within each domain; never rescale each sample.
    left_scale = .038 / 100
    right_scale = min(.038 / max(weight.sum(axis=0).max(), 1),
                      .052 / max(weight.sum(axis=1).max(), 1))
    class_total = cls[order].sum(axis=1)
    heights = class_total * left_scale
    gap = (0.94 - heights.sum()) / max(len(cls)-1, 1)
    class_y = {}; top = .97
    for name, height in heights.items():
        class_y[name] = top-height/2
        top -= height+gap
    mir_y = dict(zip(raw.index, np.linspace(.95, .05, len(raw))))
    records = []; tips = {}
    def ribbon(x0, x1, y0, y1, height, color, label):
        k = len(tips); gid = f'sankey_link_{k}'; tips[gid] = label
        dx=(x1-x0)*.48
        verts=[(x0,y0),(x0+dx,y0),(x1-dx,y1),(x1,y1),
               (x1,y1+height),(x1-dx,y1+height),(x0+dx,y0+height),(x0,y0+height),(x0,y0)]
        codes=[Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.LINETO,
               Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.CLOSEPOLY]
        patch=PathPatch(Path(verts,codes),facecolor=color,alpha=.40,edgecolor='none',gid=gid)
        ax.add_patch(patch)
    def bar(x, y, height, color):
        ax.add_patch(Rectangle((x-.003,y-height/2),.006,height,facecolor=color,edgecolor='none',zorder=4))
    sample_left = {s:y-.019 for s,y in zip(order,sample_y)}
    for i, name in enumerate(cls.index):
        yy=class_y[name]; cursor=yy-heights[name]/2
        for s in order:
            value=float(cls.loc[name,s]); height=value*left_scale
            if value > 0:
                ribbon(.205,.448,cursor,sample_left[s],height,palette[i],f'{name} → {s}: {value:.4f}% bacterial relative abundance')
                if name == 'Other bacterial classes':
                    # Qualitative tracing aid over the quantitatively scaled ribbon.
                    # Its fixed dashed stroke is not an abundance-width encoding.
                    x0, x1 = .205, .448
                    dx = (x1-x0)*.48
                    y0, y1 = cursor+height/2, sample_left[s]+height/2
                    guide = PathPatch(Path([(x0,y0),(x0+dx,y0),(x1-dx,y1),(x1,y1)],
                                           [Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4]),
                                      facecolor='none',edgecolor='#454545',lw=1.0,
                                      linestyle=(0,(4,3)),zorder=6,gid=f'other_guide_{s}')
                    ax.add_patch(guide)
                    tips[f'other_guide_{s}'] = f'Other bacterial classes → {s}: {value:.4f}%; dashed guide has no quantitative width meaning'
                records.append(dict(domain='bacterial',source=name,target=s,value=value,unit='bacterial_percent',signed_delta=np.nan))
                cursor+=height; sample_left[s]+=height
        bar(.202,yy,heights[name],palette[i])
        label = name
        if name == 'Other bacterial classes':
            positive = cls.loc[name][cls.loc[name] > 0]
            detail = (f'{positive.iloc[0]:.4f}% in {positive.index[0]}'
                      if len(positive) == 1 else f'{len(positive)} samples')
            label += '\n' + detail + ' (dashed guide)'
        ax.text(.194,yy,label,ha='right',va='center',fontsize=10)
    mir_cursor={m:mir_y[m]-weight.loc[m].sum()*right_scale/2 for m in raw.index}
    for s, yy in zip(order,sample_y):
        total=weight[s].sum()*right_scale; cursor=yy-total/2
        for m in raw.index:
            value=float(weight.loc[m,s]); signed=float(delta.loc[m,s])
            if value > 0:
                col=('#2563eb' if signed>=0 else '#93c5fd') if m.startswith('bkv-') else ('#c43d3d' if signed>=0 else '#238b45')
                height=value*right_scale
                ribbon(.552,.803,cursor,mir_cursor[m],height,col,f'{s} → {m}: signed Δ={signed:.4f}, width |Δ|={value:.4f} log2-count units')
                records.append(dict(domain='miRNA',source=s,target=m,value=value,unit='absolute_log2_normalized_count_deviation',signed_delta=signed))
                cursor+=height; mir_cursor[m]+=height
        col='#777777' if s in cases else '#c8c8c8'
        bar(.451,yy,.038,col);bar(.549,yy,total,col)
        ax.plot([.46,.54],[yy,yy],color='#b5b5b5',lw=.7,ls=':',zorder=1)
        ax.text(.5,yy,('Case ' if s in cases else 'Ctrl ')+s,ha='center',va='center',fontsize=9,
                bbox=dict(facecolor='white',edgecolor='none',pad=1.5),zorder=5)
    for m in raw.index:
        bar(.806,mir_y[m],weight.loc[m].sum()*right_scale,'#495057')
        ax.text(.814,mir_y[m],m+('*' if m=='bkv-miR-B1-3p' else ''),va='center',fontsize=10)
    ax.text(.202,1.015,'Bacterial classes',ha='center',fontweight='bold',fontsize=12)
    ax.text(.5,1.015,'Same samples',ha='center',fontweight='bold',fontsize=12)
    ax.text(.806,1.015,'Historical miRNA candidates',ha='center',fontweight='bold',fontsize=12)
    fig.suptitle('Figure 3   Descriptive Sankey of urinary profiles',x=.5,y=.98,fontsize=17,fontweight='bold')
    fig.text(.5,.945,'Two independently scaled link sets joined by sample identity',ha='center',fontsize=12)
    fig.text(.5,.918,'Other: each class stays below 1% in every sample; dashed guide shows connection only, not abundance',ha='center',fontsize=10)
    # The scale examples use precisely the same axes-coordinate width conversions.
    for x0,x1,scale,label in [(.23,.31,left_scale,'10 percentage points'),(.64,.72,right_scale,'10 normalized log2 units')]:
        val=10
        ax.add_patch(Rectangle((x0,-.065),x1-x0,val*scale,facecolor='#606060',clip_on=False))
        ax.text((x0+x1)/2,-.083,label,ha='center',va='top',fontsize=9,clip_on=False)
    fig.text(.21,.062,'Left: width = bacterial relative abundance (%)',fontsize=10,ha='center')
    fig.text(.72,.062,'Right: width = |sample log2(normalized count + 1) − control mean|; |Δ| ≥ 0.5',fontsize=10,ha='center')
    fig.text(.5,.036,'Right-link colors: human above / below control mean = red / green; viral above / below = dark / light blue.',ha='center',fontsize=10)
    fig.text(.5,.015,'Widths are comparable only within each side. Links show shared sample membership, not interactions. *Shared BKPyV/JCPyV 3p sequence.',ha='center',fontsize=10)
    stem=figdir/'Figure_3_sankey'
    fig.savefig(stem.with_suffix('.png'),dpi=220)
    fig.savefig(stem.with_suffix('.svg'))
    plt.close(fig)
    # Self-contained browser version adds exact-value tooltips to the SVG ribbons.
    tree=ET.parse(stem.with_suffix('.svg')); svg=tree.getroot()
    for element in svg.iter():
        gid=element.attrib.get('id')
        if gid in tips:
            title=ET.Element('{http://www.w3.org/2000/svg}title');title.text=tips[gid];element.insert(0,title)
    tree.write(stem.with_suffix('.svg'),encoding='utf-8',xml_declaration=True)
    svgtext=ET.tostring(svg,encoding='unicode')
    stem.with_suffix('.html').write_text('<!doctype html><html><head><meta charset="utf-8"><title>BKPyV descriptive Sankey</title><style>body{margin:16px;font-family:Arial}svg{width:100%;height:auto}svg g[id^="sankey_link"]:hover{filter:brightness(.65)}</style></head><body><p>Hover over a ribbon for its exact value. The two sides have independent units and scales.</p>'+svgtext+'</body></html>')
    pd.DataFrame(records).to_csv(out/'sankey_links.csv',index=False)
    assert np.allclose(pd.DataFrame(records).query("domain == 'bacterial'").groupby('target').value.sum(),100)
    assert np.isclose(sum(r['value'] for r in records if r['domain']=='miRNA'),weight.to_numpy().sum())
