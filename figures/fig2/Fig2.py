# -*- coding: utf-8 -*-
"""Genus-level schematic of long-read adoption and assembly contiguity.
Genus-level analogue of build_figure.py (which is order-level).
Input : Supplementary File 1 v2.xlsx  (sheet '1. NCBI and JGI data')
Output: Schematic_genus_summary.pdf / .png
Run   : & "C:\\Users\\Krolle\\AppData\\Local\\anaconda3\\python.exe" build_genus_figure.py
"""
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

plt.rcParams.update({
    'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':0.6,
    'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'
})

SRC = r'V:\Fungal_long_read_review\Nature_Genetics_Submission\Supplementary File 1 v2.xlsx'
OUTDIR = r'V:\Fungal_long_read_review'

df = pd.read_excel(SRC, sheet_name='1. NCBI and JGI data')
df['N50']=df['Assembly_Stats_Contig_N50']/1e6
df['is_LR']=df['Seq_type']=='long_read'

def agg(g):
    n=len(g); nlr=int(g.is_LR.sum())
    lr=g.loc[g.is_LR,'N50']; sr=g.loc[~g.is_LR,'N50']
    return pd.Series({'phylum':g['Current_phylum'].mode().iloc[0],'n':n,'nlr':nlr,
        'pctLR':100*nlr/n,'lrN50':lr.median() if nlr else np.nan,
        'srN50':sr.median() if (n-nlr) else np.nan})

o=df.groupby('Current_genus',group_keys=True).apply(agg,include_groups=False)
o['fold']=o['lrN50']/o['srN50']
o=o.dropna(subset=['fold'])                      # need both categories present

# Selection: 18 most-sequenced genera (with both categories) ...
base=o[o.n>=20].sort_values('n',ascending=False).head(18)
# ... plus striking high-gain / early-diverging outliers to show the decoupling
force=['Coemansia','Mortierella','Rhizopus']
extra=o.loc[[g for g in force if g in o.index and g not in base.index]]
sel=pd.concat([base,extra])
sel=sel.sort_values('fold',ascending=True)       # largest fold gain at top
names=list(sel.index)
y=np.arange(len(names))

PH={'Ascomycota':'#0072B2','Basidiomycota':'#E69F00'}
def pcol(p): return PH.get(p,'#009E73')
PURPLE='#6a3d9a'; PINK='#e7298a'

fig=plt.figure(figsize=(9.2,7.3))
gs=fig.add_gridspec(1,2,width_ratios=[2.15,1.0],wspace=0.05)
axA=fig.add_subplot(gs[0]); axB=fig.add_subplot(gs[1],sharey=axA)

# ---- Panel a: dumbbell SR -> LR contig N50 ----
axA.set_xscale('log')
for i,name in enumerate(names):
    r=sel.loc[name]
    axA.annotate('',xy=(r.lrN50,i),xytext=(r.srN50,i),
        arrowprops=dict(arrowstyle='-|>',color='#9e9e9e',lw=1.1,
                        shrinkA=3.5,shrinkB=3.5,mutation_scale=9))
axA.scatter(sel.srN50,y,s=34,color=PINK,zorder=3,edgecolor='white',linewidth=0.6)
axA.scatter(sel.lrN50,y,s=34,color=PURPLE,zorder=3,edgecolor='white',linewidth=0.6)
for i,name in enumerate(names):
    r=sel.loc[name]
    axA.text(r.lrN50*1.35,i,f'{r.fold:.0f}×',va='center',ha='left',fontsize=6.7,color='#333333')
axA.axvline(1.0,color='#bdbdbd',ls=':',lw=0.8,zorder=0)
axA.text(1.02,len(names)-0.3,'1 Mbp',fontsize=6,color='#999999',rotation=90,va='top')
axA.set_xlim(0.002,60); axA.set_ylim(-0.8,len(names)-0.2)
axA.set_yticks(y); axA.set_yticklabels(names,fontsize=7.5,fontstyle='italic')
for tick,name in zip(axA.get_yticklabels(),names):
    tick.set_color(pcol(sel.loc[name,'phylum']))
axA.set_xlabel('Median contig N50 (Mbp, log scale)',fontsize=8)
axA.set_title('a  Contiguity jump from short- to long-read',loc='left',fontsize=9.5,fontweight='bold',pad=6)
for s in ['top','right']: axA.spines[s].set_visible(False)
axA.tick_params(axis='x',labelsize=7)
axA.grid(axis='x',color='#eeeeee',lw=0.5,zorder=0)

# ---- Panel b: long-read adoption % ----
axB.barh(y,sel.pctLR,height=0.62,color=[pcol(p) for p in sel.phylum],zorder=3)
for i,name in enumerate(names):
    axB.text(sel.loc[name,'pctLR']+1.5,i,f'{sel.loc[name,"pctLR"]:.0f}%',va='center',ha='left',fontsize=6.7,color='#333333')
axB.axvline(20.0,color='#757575',ls='--',lw=0.8,zorder=2)
axB.text(20,len(names)-0.2,'20% overall',fontsize=6,color='#555555',ha='center',va='bottom')
axB.set_xlim(0,100)
axB.set_xlabel('Long-read adoption (%)',fontsize=8)
axB.set_title('b  Adoption',loc='left',fontsize=9.5,fontweight='bold',pad=6)
plt.setp(axB.get_yticklabels(),visible=False)
axB.tick_params(axis='y',length=0)
for s in ['top','right','left']: axB.spines[s].set_visible(False)
axB.tick_params(axis='x',labelsize=7)
axB.grid(axis='x',color='#eeeeee',lw=0.5,zorder=0)

# ---- legends ----
legA=[Line2D([0],[0],marker='o',color='w',markerfacecolor=PINK,markersize=7,label='Short-read median N50'),
      Line2D([0],[0],marker='o',color='w',markerfacecolor=PURPLE,markersize=7,label='Long-read median N50'),
      Line2D([0],[0],color='#9e9e9e',lw=1.1,label='n× = fold gain')]
axA.legend(handles=legA,loc='lower right',fontsize=6.8,frameon=False,handletextpad=0.4)
legB=[Patch(facecolor=PH['Ascomycota'],label='Ascomycota'),
      Patch(facecolor=PH['Basidiomycota'],label='Basidiomycota'),
      Patch(facecolor='#009E73',label='Early-diverging phyla')]
axB.legend(handles=legB,loc='lower right',fontsize=6.6,frameon=False,handletextpad=0.5,title='Phylum',title_fontsize=6.8)

# ---- header ----
fig.suptitle('Long-read contiguity gains and adoption across the most-sequenced fungal genera',
             fontsize=11.0,fontweight='bold',x=0.5,y=0.985)
fig.text(0.5,0.948,
    '18 most-sequenced genera plus three high-gain early-diverging outliers (Supplementary File 1 v2).',
    ha='center',fontsize=7.6,color='#444444')
fig.text(0.5,0.928,
    'As at order level, the genera that gain most in contiguity are often the least sequenced with long reads.',
    ha='center',fontsize=7.6,color='#444444',fontstyle='italic')

# callout: decoupling
axA.annotate('Largest gains yet low adoption:\nCoemansia 137× at 1% LR;\nPyricularia (rice blast) 108×',
    xy=(1.15,len(names)-1.05),xytext=(3.0,len(names)-1.9),
    fontsize=6.7,color='#7a0177',ha='left',va='center',
    bbox=dict(boxstyle='round,pad=0.3',fc='#fdf2fb',ec='#c994c7',lw=0.6),
    arrowprops=dict(arrowstyle='-',color='#c994c7',lw=0.7))

fig.subplots_adjust(left=0.145,right=0.965,top=0.855,bottom=0.115)
for ext in ['pdf','png']:
    fig.savefig(rf'{OUTDIR}\Schematic_genus_summary.{ext}',dpi=300 if ext=='png' else None)

print('rows:',len(names))
print(sel[['phylum','n','nlr','pctLR','srN50','lrN50','fold']].round(2).to_string())
print('saved Schematic_genus_summary.pdf/.png to',OUTDIR)
