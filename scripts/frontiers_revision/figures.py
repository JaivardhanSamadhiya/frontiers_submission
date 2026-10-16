"""Publication figures from verified saved results, without model fitting."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, precision_recall_curve
from PIL import Image
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/frontiers_revision_2026_10_14/figures'
AUDIT=ROOT/'data/results_v2/frontiers_revision_2026_10_14/cocktail_reporting'


def save(fig,name):
    for line in fig.findobj(Line2D):line.set_linewidth(max(2.2,line.get_linewidth()))
    for ax in fig.axes:
        title=ax.get_title()
        if len(title)>2 and title[0] in 'ABCD' and title[1]==' ':
            ax.set_title(f'({title[0]}) {title[1:].strip()}')
        for spine in ax.spines.values():spine.set_linewidth(2.2)
        for collection in ax.collections:
            if hasattr(collection,'set_linewidths'):
                collection.set_linewidths(np.maximum(2.2,collection.get_linewidths()))
    fig.savefig(OUT/(name+'.png'),dpi=400,bbox_inches='tight')
    fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight')
    with Image.open(OUT/(name+'.png')) as source:
        width=round(180/25.4*400)
        rgb=source.convert('RGB');height=round(rgb.height*width/rgb.width)
        rgb=rgb.resize((width,height),Image.Resampling.LANCZOS)
        rgb.save(OUT/(name+'.tiff'),compression='tiff_lzw',dpi=(400,400))
    plt.close(fig)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,
                         'axes.labelsize':10,'xtick.labelsize':9,'ytick.labelsize':9,
                         'legend.fontsize':9,'axes.spines.top':False,'axes.spines.right':False})
    plt.rcParams['pdf.fonttype']=42
    data=json.loads((ROOT/'data/results_v2/significance_results.json').read_text())['skill_vs_chance']
    regimes=['loso_species','host_cluster','phage_cluster','combined_unseen']
    labels=['Eligible host species','Host sequence cluster','Phage sequence cluster','Both entity clusters']
    fig,axes=plt.subplots(2,1,figsize=(6.7,5.8),layout='constrained')
    for model,shift,color in [('GBM',-.1,'#2468a0'),('GNN',.1,'#ae5630')]:
        scores=[next(r['auc'] for r in data if r['regime']==s and r['model']==model) for s in regimes]
        axes[0].plot(scores,np.arange(4)+shift,'o-' if model=='GBM' else 's--',label=model,color=color)
    axes[0].set(yticks=range(4),yticklabels=labels,xlim=(.5,1),ylim=(-.4,3.6),xlabel='Pooled AUROC',title='A  Historical test cohorts')
    axes[0].invert_yaxis();axes[0].legend(loc='lower left');axes[0].grid(axis='x',alpha=.2)
    sensitivity=pd.read_csv(ROOT/'data/results_v2/canonical_sensitivity/summary.csv')
    order=['shared_entities_allowed','host_cluster_excluded','phage_cluster_excluded','both_clusters_excluded']
    for i,name in enumerate(order):
        v=sensitivity[sensitivity.regime==name].auroc
        axes[1].plot(v,np.full(len(v),i),'o',color='#2468a0',alpha=.7)
        axes[1].plot([v.min(),v.max()],[i,i],color='#2468a0')
    axes[1].set(yticks=range(4),yticklabels=['Shared entities','Host excluded','Phage excluded','Both excluded'],xlim=(.5,1),ylim=(-.4,3.6),xlabel='Pooled AUROC',title='B  Canonical cutoff 0.05, fixed cohort')
    axes[1].invert_yaxis();axes[1].grid(axis='x',alpha=.2)
    save(fig,'Figure_1_validation')
    frame=pd.read_csv(AUDIT/'host_scores.csv')
    y=frame.productive.to_numpy(int)
    methods=[('historical_product','Historical product','#2468a0','-'),
             ('maximum_constituent','Maximum constituent','#bd6728','--'),
             ('logistic_component','Logistic component','#28846b','-.'),
             ('host_phylogroup','Host phylogroup','#666666',':')]
    fig,axes=plt.subplots(1,2,figsize=(6.7,4.1))
    for method,label,color,style in methods:
        fpr,tpr,_=roc_curve(y,frame[method])
        precision,recall,_=precision_recall_curve(y,frame[method])
        axes[0].plot(fpr,tpr,label=label,color=color,ls=style,lw=1.7)
        axes[1].step(recall,precision,where='post',color=color,ls=style,lw=1.7)
    axes[0].plot([0,1],[0,1],ls=':',color='.7',lw=1)
    axes[1].axhline(.82,ls='--',color='.25',lw=1)
    axes[1].text(.08,.79,'Prevalence = 0.820',va='top',fontsize=9)
    axes[0].set(xlabel='False positive rate',ylabel='True positive rate',title='A  ROC',xlim=(0,1),ylim=(0,1.02))
    axes[1].set(xlabel='Recall',ylabel='Precision',title='B  Precision-recall',xlim=(0,1),ylim=(0,1.02))
    for ax in axes:ax.grid(alpha=.15)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False)
    fig.subplots_adjust(left=.09,right=.985,bottom=.25,top=.9,wspace=.34)
    save(fig,'Figure_2_cocktails')
    fi=pd.read_csv(ROOT/'data/results_v2/feature_importance.csv').iloc[:6].sort_values('importance')
    fi['feature']=fi.feature.replace({'Phage genome composition':'Phage genome descriptors',
        'Host genome composition':'Host genome descriptors',
        'Nucleotide homology':'Nucleotide match descriptors',
        'Composition distance':'Composition comparisons',
        'Protein homology':'Translated-word proxies','CRISPR matching':'CRISPR-like proxies'})
    fig,ax=plt.subplots(figsize=(6.7,3.6),layout='constrained')
    ax.barh(fi.feature,fi.importance,color=['#6e91a8']*len(fi))
    ax.set(xlabel='Fraction of normalized XGBoost gain',xlim=(0,.72))
    for y,v in enumerate(fi.importance):ax.text(v+.008,y,f'{100*v:.1f}%',va='center',fontsize=10)
    ax.grid(axis='x',alpha=.2)
    save(fig,'Figure_3_importance')
    edge=pd.read_csv(ROOT/'data/results_v2/feature_importance.csv').iloc[6:].sort_values('importance')
    fig,ax=plt.subplots(figsize=(6.7,6.7),layout='constrained')
    ax.barh(edge.feature,edge.importance,color='#6e91a8')
    ax.set(xlabel='Fraction of normalized total XGBoost gain')
    ax.grid(axis='x',alpha=.2)
    save(fig,'Supplementary_Figure_S1_pair_importance')
    triage=pd.read_csv(AUDIT/'triage_ranked_records.csv')
    selected=triage[(triage.cohort=='historical_100') & (triage.method=='historical_product')]
    fig,axes=plt.subplots(1,2,figsize=(6.7,3.8),layout='constrained')
    axes[0].step(selected['rank'],selected['precision'],where='post',color='#2468a0',label='Historical score')
    axes[0].axhline(.82,color='.5',ls='--',label='Random expectation')
    axes[0].set(xlabel='Records screened',ylabel='Precision',xlim=(0,100),ylim=(0,1.04),title='A  Precision at rank')
    axes[1].step(selected['rank'],selected['recall'],where='post',color='#2468a0')
    axes[1].plot([0,100],[0,1],color='.5',ls='--')
    axes[1].set(xlabel='Records screened',ylabel='Fraction of positives found',xlim=(0,100),ylim=(0,1.04),title='B  Observed-positive recovery')
    axes[0].legend(loc='lower left',frameon=False)
    for ax in axes:ax.grid(alpha=.15)
    save(fig,'Supplementary_Figure_S2_triage')
    calibration=pd.read_csv(AUDIT.parent/'remaining_checks/recalibration_summary.csv')
    budgets=pd.read_csv(AUDIT.parent/'remaining_checks/recalibration_label_budgets.csv')
    fig,axes=plt.subplots(1,2,figsize=(6.7,4.0))
    for method,label,color,style in [('platt','Platt','#2468a0','-'),('isotonic','Isotonic','#bd6728','--'),('target_prior','Target prior','#666666',':')]:
        fractions=sorted(calibration.fraction.unique())
        sizes=[budgets[budgets.fraction==fraction].train_rows.mean() for fraction in fractions]
        for ax,key in zip(axes,['ece','brier']):
            means=[];low=[];high=[]
            for fraction in fractions:
                values=calibration[(calibration.fraction==fraction)&(calibration.method==method)][key]
                means.append(values.mean());low.append(max(0.,values.mean()-values.min()));high.append(max(0.,values.max()-values.mean()))
            ax.errorbar(sizes,means,yerr=[low,high],color=color,ls=style,marker='o',markersize=3,capsize=2,lw=1.4,label=label)
    for ax,key,title in zip(axes,['ece','brier'],['A  Calibration error','B  Brier score']):
        raw=calibration[calibration.method=='uncalibrated'][key].mean()
        ax.axhline(raw,color='#9d3543',ls='-.',label='Raw frozen score')
        ax.set(xlabel='Calibration records per fold (mean)',ylabel='10-bin ECE' if key=='ece' else 'Brier score',title=title)
        ax.set_xticks([100,400,800]);ax.grid(alpha=.15)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False)
    fig.subplots_adjust(left=.09,right=.985,bottom=.28,top=.9,wspace=.32)
    save(fig,'Supplementary_Figure_S3_recalibration')
    print('Rendered six final figure sets from preserved results.')


if __name__=='__main__':main()
