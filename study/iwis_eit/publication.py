"""Saved-result summaries and figures only: no FEM, reconstruction, or RNG."""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from .monte_carlo import summary
from .figure_config import configure_figure, figure_size_inches, boxplot_y_limits, circle_label


def table_rows(data,metadata):
    s=summary(data['single_correlations'],data['unmixed_correlations'])
    rows=[]
    for i,name in enumerate(data['inclusion_names']):
        row=dict(position=metadata['position'],noise_level=metadata['noise_level'],lambda_value=metadata['lambda_value'],
                 mode=metadata['mode'],trial_count=metadata['trial_count'],inclusion=circle_label(name))
        for method in ['single','unmixed']:
            for statistic in ['median','q1','q3']:
                row[f'{method}_{statistic}']=s[method][statistic][i]
        rows.append(row)
    return rows


def markdown_table(rows):
    lines=['| d (m) | Noise η | Inclusion | Single median [Q1, Q3] | Unmixed median [Q1, Q3] | n |',
           '| --- | --- | --- | --- | --- | --- |']
    for row in rows:
        cells=[f"{row['position']:g}",f"{row['noise_level']:g}",row['inclusion']]
        for method in ['single','unmixed']:
            cells.append(f"{row[method+'_median']:.3f} [{row[method+'_q1']:.3f}, {row[method+'_q3']:.3f}]")
        cells += [str(row['trial_count'])]
        lines.append('| '+' | '.join(cells)+' |')
    return '\n'.join(lines)


def export_tables(rows,directory):
    directory=Path(directory)
    directory.mkdir(parents=True,exist_ok=True)
    with (directory/'summary.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (directory/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')
    header=f"# Computed results — {rows[0]['mode']}, λ={rows[0]['lambda_value']:g}\n\n"
    note='\n\nBrackets are Q1/Q3 (linear interpolation); undefined correlations are omitted from quantiles.\n'
    (directory/'summary.md').write_text(header+markdown_table(rows)+note)


def condition_title(m):
    return f"{m['mode'].upper()} | d={m['position']:g} m | η={m['noise_level']:g} | λ={m['lambda_value']:g} | n={m['trial_count']}"


def boxplots(data,m,*,config):
    fig,ax=plt.subplots(figsize=figure_size_inches(config['figure_sizes']['boxplot']),layout='constrained')
    colors=plt.rcParams['axes.prop_cycle'].by_key()['color'][:2]
    for offset,key,color in [(-.18,'single_correlations',colors[0]),(.18,'unmixed_correlations',colors[1])]:
        values=[v[np.isfinite(v)] for v in data[key].T]
        bp=ax.boxplot(values,positions=np.arange(1,4)+offset,widths=.28,patch_artist=True,
                      medianprops=dict(color='black'),flierprops=dict(marker='.',markersize=3))
        for patch in bp['boxes']:patch.set(facecolor=color)
    ax.set(xticks=[1,2,3],xticklabels=[circle_label(name) for name in data['inclusion_names']],
           ylabel='Absolute Pearson correlation |r|',title='Single versus Unmixed')
    y_limits=boxplot_y_limits(config,m['position'])
    if y_limits is not None:ax.set_ylim(y_limits)
    ax.legend(handles=[Patch(facecolor=c,label=n) for c,n in zip(colors,['Single','Unmixed'])],loc='lower left')
    ax.grid(axis='y',alpha=.2)
    fig.suptitle(condition_title(m))
    fig.supxlabel('Boxes: Q1/Q3; line: median; whiskers: 1.5 IQR; points: outliers.',fontsize=8)
    return fig


def combined_boxplot(conditions,*,config):
    """Group by circle, then noise; each noise has adjacent Single/Unmixed boxes."""
    if not conditions:
        raise ValueError('Select at least one noise condition')
    first_data,first_metadata=conditions[0]
    names=first_data['inclusion_names']
    noises=[metadata['noise_level'] for _,metadata in conditions]
    if len(set(noises))!=len(noises):
        raise ValueError('Select each noise level once')
    for data,metadata in conditions:
        if (metadata['position']!=first_metadata['position'] or
                metadata['lambda_value']!=first_metadata['lambda_value'] or
                metadata['trial_count']!=first_metadata['trial_count'] or
                not np.array_equal(data['inclusion_names'],names)):
            raise ValueError('Combined plot requires matching distance, lambda, trials and inclusion order')
    fig,ax=plt.subplots(figsize=figure_size_inches(config['figure_sizes']['combined_boxplot']),layout='constrained')
    colors=plt.rcParams['axes.prop_cycle'].by_key()['color'][:2]
    group_spacing=.8*len(noises)+1.0
    centers=np.arange(len(names))*group_spacing
    noise_offsets=(np.arange(len(noises))-(len(noises)-1)/2)*.8
    for circle,center in enumerate(centers):
        for noise_offset,(data,metadata) in zip(noise_offsets,conditions):
            pair_center=center+noise_offset
            for method_offset,key,color in [(-.15,'single_correlations',colors[0]),
                                             (.15,'unmixed_correlations',colors[1])]:
                values=data[key][:,circle]
                values=values[np.isfinite(values)]
                bp=ax.boxplot([values],positions=[pair_center+method_offset],widths=.24,
                              patch_artist=True,manage_ticks=False,
                              medianprops=dict(color='black'),
                              flierprops=dict(marker='.',markersize=3))
                bp['boxes'][0].set(facecolor=color)
            ax.text(pair_center,-.025,f"η={metadata['noise_level']:g}",
                    transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=9)
    ax.set(xticks=centers,xticklabels=[circle_label(name) for name in names],
           xlim=(centers[0]-.4*len(noises)-.25,centers[-1]+.4*len(noises)+.25),
           ylabel='Absolute Pearson correlation |r|',title='Single versus Unmixed across noise levels')
    y_limits=boxplot_y_limits(config,first_metadata['position'])
    if y_limits is not None:ax.set_ylim(y_limits)
    ax.tick_params(axis='x',length=0,pad=29,labelsize=11)
    for boundary in (centers[:-1]+centers[1:])/2:
        ax.axvline(boundary,color='0.85',linewidth=.8)
    ax.grid(axis='y',alpha=.2)
    ax.legend(handles=[Patch(facecolor=c,label=n) for c,n in zip(colors,['Single','Unmixed'])],loc='lower left')
    fig.suptitle(f"{first_metadata['mode'].upper()} | d={first_metadata['position']:g} m | "
                 f"λ={first_metadata['lambda_value']:g} | n={first_metadata['trial_count']} per noise level")
    fig.supxlabel('Within each circle: noise-level pairs, with Single on the left and Unmixed on the right.',fontsize=9)
    return fig


def export_figure(fig,stem,config):
    stem=Path(stem)
    configure_figure(fig,config)
    with plt.rc_context({'svg.fonttype':'none','svg.hashsalt':'IWIS2026'}):
        fig.savefig(Path(str(stem)+'.svg'),dpi=config['png_dpi'],metadata={'Date':None})
    fig.savefig(Path(str(stem)+'.png'),dpi=config['png_dpi'],metadata={'Software':'IWIS2026'})
    return [Path(str(stem)+'.svg'),Path(str(stem)+'.png')]
