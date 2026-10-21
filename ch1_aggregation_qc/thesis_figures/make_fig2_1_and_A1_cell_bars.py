import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
C={'Breast':'#fdbf6f','Lung':'#cab2d6','Ovarian':'#b09c85','Colorectal':'#a6cee3','Other models':'#91d1c2','Unknown':'#80c680',
   'Pancreas':'#ffdc91','Skin':'#dc0000','Prostate':'#8491b4','Kidney':'#7e6148','Brain':'#3c5488','Liver/Biliary':'#f39b7f',
   'Head and Neck':'#e64b35','Neuroendocrine':'#00a087','Sarcoma':'#4dbbd5'}
pre=[('Breast',5.30,35.2,98),('Lung',2.09,13.9,45),('Ovarian',2.00,13.3,26),('Colorectal',1.33,8.8,35),('Other models',0.94,6.2,33),
     ('Unknown',0.57,3.8,29),('Pancreas',0.51,3.4,19),('Skin',0.46,3.0,17),('Prostate',0.45,3.0,20),('Kidney',0.41,2.7,9),
     ('Brain',0.28,1.9,14),('Liver/Biliary',0.27,1.8,10),('Head and Neck',0.24,1.6,5),('Neuroendocrine',0.14,0.9,4),('Sarcoma',0.08,0.5,3)]
post=[('Breast',2.78,34.2,59),('Lung',1.08,13.3,24),('Colorectal',0.98,12.1,20),('Unknown',0.57,7.1,29),('Pancreas',0.39,4.8,12),
      ('Ovarian',0.39,4.7,12),('Kidney',0.37,4.6,7),('Skin',0.37,4.6,10),('Other models',0.31,3.8,4),('Prostate',0.23,2.8,10),
      ('Liver/Biliary',0.19,2.3,7),('Brain',0.18,2.2,3),('Neuroendocrine',0.13,1.6,4),('Sarcoma',0.08,1.0,3),('Head and Neck',0.07,0.9,1)]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.linewidth':0.8,'xtick.labelsize':8.5,'ytick.labelsize':9.5})
def draw(data,out,xmax,ticks):
    fig,ax=plt.subplots(figsize=(5.39,3.75))
    n=len(data); ys=list(range(n))[::-1]
    for y,(name,m,p,proj) in zip(ys,data):
        ax.barh(y,m,height=0.72,color=C[name],edgecolor='none')
        ax.text(m+xmax*0.012,y,f"{m:.2f}M ({p:.1f}%), {proj} project{'s' if proj!=1 else ''}",va='center',ha='left',fontsize=8.5,color='#222222')
    ax.set_yticks(ys); ax.set_yticklabels([d[0] for d in data])
    ax.set_xlim(0,xmax); ax.set_xticks(ticks); ax.set_ylim(-0.6,n-0.4)
    ax.set_xlabel('Number of cells (millions)',fontsize=9.5)
    for s in ('top','right'): ax.spines[s].set_visible(False)
    ax.tick_params(axis='y',length=0,pad=4); ax.tick_params(axis='x',length=3)
    fig.tight_layout(pad=0.3)
    fig.savefig(out,dpi=400); plt.close(fig)
draw(pre,'figs/ch1_fig1.png',9.2,[0,1,2,3,4,5,6,7,8,9])
draw(post,'figs/ch1_fig2.png',4.8,[0,0.5,1,1.5,2,2.5,3,3.5,4,4.5])
print('ok')
