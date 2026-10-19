"""Publication-sized stacked panels from verified saved scientific metrics."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'data/results_v2/frontiers_strengthening_2026_10_16/doublets'


def draw():
    decision=pd.read_csv(OUT/'doublet_decision_summary.csv',float_precision='round_trip')
    conditional=pd.read_csv(OUT/'triplet_conditional_summary.csv',float_precision='round_trip')
    plt.rcParams.update({'font.size':11,'pdf.fonttype':42,'ps.fonttype':42,
        'lines.linewidth':2.5,'axes.linewidth':2.5,'xtick.major.width':2.5,'ytick.major.width':2.5})
    # Two-column production width: 180 mm; embedded labels remain above 8 pt.
    width=180/25.4
    fig,axes=plt.subplots(2,1,figsize=(width,width*1.15))
    fig.subplots_adjust(left=.29,right=.97,bottom=.13,top=.95,hspace=.35)
    groups=[(['repaired_product','historical_product','repaired_max','repaired_mean','raw_logistic_product','raw_knn_product','training_prevalence_mean'],
             ['Repaired product','Historical product','Repaired maximum','Repaired mean','Raw logistic product','Raw neighbor product','Training susceptibility']),
            (['historical_product','repaired_product','historical_max','historical_mean','raw_logistic_product','raw_knn_product','training_prevalence_mean'],
             ['Historical product','Repaired product','Historical maximum','Historical mean','Raw logistic product','Raw neighbor product','Training susceptibility'])]
    for index,(methods,labels) in enumerate(groups):
        source,metric=(decision,'gain_vs_uniform') if index==0 else (conditional,'pair_weighted_conditional_auroc')
        scale=100 if index==0 else 1
        for offset,resampling,color,marker in [(-.12,'host_cluster','#21618C','o'),(.12,'crossed','#666666','s')]:
            part=source.loc[source.resampling.eq(resampling)&source.metric.eq(metric)].set_index('method').loc[methods]
            axes[index].errorbar(part.estimate*scale,np.arange(len(methods))+offset,
                xerr=np.vstack([part.estimate-part.ci_low,part.ci_high-part.estimate])*scale,
                fmt=marker,color=color,capsize=3,markersize=4,label='Host cluster' if resampling=='host_cluster' else 'Crossed host/composition')
        axes[index].set_yticks(np.arange(len(methods)),labels)
        axes[index].invert_yaxis()
        axes[index].axvline(0 if index==0 else .5,color='#888888',linestyle='--',linewidth=2.5)
        axes[index].spines[['top','right']].set_visible(False)
        axes[index].tick_params(labelsize=10)
        axes[index].set_title('(A) Choice among measured doublets' if index==0 else '(B) Host discrimination within the same triplet',fontsize=12,pad=15)
        axes[index].set_xlabel('Activity gain over uniform choice\n(percentage points)' if index==0 else 'Within-composition concordance',fontsize=11)
    axes[1].set_xlim(-.03,1.03)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.62,.004),ncol=1,frameon=False,fontsize=10)
    for extension,kwargs in [('png',{'dpi':300}),('pdf',{}),('tiff',{'dpi':300,'pil_kwargs':{'compression':'tiff_lzw'}})]:
        fig.savefig(OUT/f'measured_menu_and_conditional_host_portrait.{extension}',**kwargs)
    path=OUT/'measured_menu_and_conditional_host_portrait.tiff'
    with Image.open(path) as source:
        rgb=source.convert('RGB')
    rgb.save(path,compression='tiff_lzw',dpi=(300,300))
    plt.close(fig)


if __name__=='__main__':draw()
