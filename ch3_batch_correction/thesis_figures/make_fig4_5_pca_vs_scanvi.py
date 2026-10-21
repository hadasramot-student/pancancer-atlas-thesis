from PIL import Image
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams['font.family']='DejaVu Sans'
src=Image.open('v117docx/word/media/image5.png').convert('RGB')
def crop(box):
    c=src.crop(box); a=np.asarray(c.convert('L'))<240
    ys,xs=np.where(a); return np.asarray(c.crop((xs.min()-8,ys.min()-8,xs.max()+9,ys.max()+9)))
pa=crop((120,300,1800,1780)); pb=crop((1900,300,3620,1780))
print(pa.shape,pb.shape)
W=5.39
fig=plt.figure(figsize=(W,3.05))
a1=fig.add_axes([0.015,0.13,0.475,0.80]); a2=fig.add_axes([0.515,0.13,0.475,0.80])
for ax,im,t in [(a1,pa,'a   Uncorrected (PCA)'),(a2,pb,'b   scANVI')]:
    ax.imshow(im,interpolation='lanczos'); ax.set_axis_off(); ax.set_title(t,fontsize=9.5,fontweight='bold',loc='left',pad=3)
cols=['#e4572e','#3b7dd8','#5bb867']
# sample real colors
arr=np.asarray(src)
labs=['3CA (87,773)','ScBaseCount (98,946)','Breastscape (13,281)']
h=[Line2D([],[],marker='o',ls='',ms=6,color=c) for c in cols]
fig.legend(h,labs,loc='lower center',ncol=3,fontsize=8.5,frameon=False,handletextpad=0.2,columnspacing=1.4,title='Source repository (cells)',title_fontsize=8.5,alignment='center',bbox_to_anchor=(0.5,0.0))
fig.savefig('c8/ch2_fig3.png',dpi=400)
