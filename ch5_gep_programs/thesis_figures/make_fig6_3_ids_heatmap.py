import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family']='DejaVu Sans'
rows=['Lung','Kidney','Breast','Pancreas','Sarcoma','Ovarian','Brain','Liver/Biliary','Skin','Prostate','Colorectal']
v=np.array([[.12,.21,.38],[.07,.20,.65],[.07,.17,.52],[.12,.16,.32],[.28,.15,.21],[.11,.14,.36],[.06,.12,.27],[.08,.12,.25],[.10,.11,.31],[.12,.10,.18],[.03,.09,.31]])
cols=['cNMF–\ncNMF','cNMF–\nGeneNMF','GeneNMF–\nGeneNMF']
fig,ax=plt.subplots(figsize=(4.3,4.5))
im=ax.imshow(v,cmap='Blues',vmin=0,vmax=0.7,aspect='auto')
for i in range(v.shape[0]):
    for j in range(3):
        ax.text(j,i,f'{v[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if v[i,j]>0.4 else '#1a1a1a')
ax.set_xticks(range(3)); ax.set_xticklabels(cols,fontsize=8)
ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows,fontsize=8.5)
ax.set_xticks(np.arange(-.5,3,1),minor=True); ax.set_yticks(np.arange(-.5,len(rows),1),minor=True)
ax.grid(which='minor',color='white',lw=1.5); ax.tick_params(which='minor',length=0); ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
ax.xaxis.tick_top()
cb=fig.colorbar(im,ax=ax,fraction=0.05,pad=0.03); cb.set_label('Median IDS',fontsize=8.5); cb.ax.tick_params(labelsize=8)
fig.tight_layout(); fig.savefig('c13/ch5_fig3.png',dpi=400)
