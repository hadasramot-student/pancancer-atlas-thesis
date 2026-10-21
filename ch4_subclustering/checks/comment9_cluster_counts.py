import csv,glob,os,collections
D='/groups/ofircohen-group/users/hadasra/temp_objects/clustered_lineages'
OUT='/groups/ofircohen-group/users/hadasra/thesis_checks/comment9_2026-10-21/cluster_counts.csv'
RES=['0.15','0.20','0.25','0.30','0.50','0.75','1.00','2.00']
rows=[]
for f in sorted(glob.glob(D+'/*_cluster_summary.tsv')):
    obj=os.path.basename(f).replace('_cluster_summary.tsv','')
    s=list(csv.DictReader(open(f),delimiter='\t'))[0]
    sc=collections.defaultdict(list)
    for r in csv.DictReader(open(f'{D}/{obj}_scimilarity_summary.tsv'),delimiter='\t'):
        sc[r['resolution'].replace('louvain_res','')].append(r)
    o={'object':obj,'lineage':obj.split('_')[0],'n_cells':s['n_cells'],'default_res_stored':s.get('default_res','')}
    for r in RES:
        o['clusters_res'+r]=s.get('res_'+r, str(len(sc[r])) if sc.get(r) else 'NA')
    c=sc.get('0.50',[])
    ids=sorted(int(x['cluster']) for x in c); ns=[int(x['n_cells']) for x in c]; cf=[float(x['confidence_pct']) for x in c]
    o['res0.50_sci_rows']=len(c)
    o['res0.50_max_id']=max(ids)
    o['res0.50_id_gaps']=len(set(range(max(ids)+1))-set(ids))
    o['res0.50_zero_cell']=sum(n==0 for n in ns)
    o['res0.50_clusters_ge10cells']=sum(n>=10 for n in ns)
    o['res0.50_clusters_lt10cells']=sum(n<10 for n in ns)
    o['res0.50_mean_conf']=round(sum(cf)/len(cf),1)
    o['res0.50_n_conf_gt80']=sum(x>80 for x in cf)
    o['in_thesis_75']=obj!='Epithelial_Brain'
    for rr in ['0.30']:
        c3=sc.get(rr,[])
        o['res0.30_clusters_ge10cells']=sum(int(x['n_cells'])>=10 for x in c3) if c3 else 'NA'
    o['flag_lt1000cells']=int(s['n_cells'])<1000
    o['flag_le3clusters_res0.50']=int(s['res_0.50'])<=3
    rows.append(o)
w=csv.DictWriter(open(OUT,'w'),fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(len(rows))
