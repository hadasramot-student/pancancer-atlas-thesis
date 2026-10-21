from PIL import Image
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams['font.family']='DejaVu Sans'
def crop(path,box):
    im=Image.open(path).convert('RGB'); c=im.crop(box); a=np.asarray(c.convert('L'))<245
    ys,xs=np.where(a); return np.asarray(c.crop((xs.min()-4,ys.min()-4,xs.max()+5,ys.max()+5)))
f1='v123/content/images/ch3_fig1.png'; f2='v123/content/images/ch3_fig2.png'
p1=[crop(f1,(200,175,600,710)),crop(f1,(370,780,700,1316))]
p2=[crop(f2,(150,180,500,715)),crop(f2,(340,785,690,1323))]
W=5.39
def panels(fig,axs,ims):
    for ax,im,t in zip(axs,ims,['Uncorrected','scANVI']):
        ax.imshow(im,interpolation='lanczos'); ax.set_axis_off(); ax.set_title(t,fontsize=9.5,fontweight='bold',pad=3)
# Fig 4.1
fig,axs=plt.subplots(1,2,figsize=(W*0.67,2.85)); panels(fig,axs,p1)
plt.subplots_adjust(left=0.01,right=0.99,top=0.92,bottom=0.01,wspace=0.08)
fig.savefig('c5/ch3_fig1_main.png',dpi=400); plt.close()
# Fig 4.2
labs=['Activated CD4+ T cell','Activated CD8+ T cell','CD4+ T cell','CD8+ cytotoxic T cell','CD8+ memory T cell','CD8+ T cell','Gamma-delta T cell','Naive CD4+ T cell','Regulatory T cell','T cell (unspecified)','Other']
cols=['#e41a1c','#377eb8','#4daf4a','#ff7f00','#984ea3','#a65628','#f781bf','#66c2a5','#fc8d62','#8da0cb','#bbbbbb']
fig=plt.figure(figsize=(W,2.85))
a1=fig.add_axes([0.0,0.01,0.33,0.91]); a2=fig.add_axes([0.34,0.01,0.33,0.91]); panels(fig,[a1,a2],p2)
h=[Line2D([],[],marker='o',ls='',ms=5.5,color=c) for c in cols]
fig.legend(h,labs,loc='center left',ncol=1,fontsize=8,frameon=False,title='Cell type',title_fontsize=8.5,
           handletextpad=0.2,labelspacing=0.45,bbox_to_anchor=(0.68,0.5),alignment='left')
fig.savefig('c5/ch3_fig2_main.png',dpi=400); plt.close()
