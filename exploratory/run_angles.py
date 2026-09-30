import hub_reduction as h,itertools,json,os,gc,time,traceback
from pathlib import Path
out=str(Path(__file__).resolve().parent/'symmetry_angles.json')
data,M,T=h.setup();results=[]
for tags in itertools.combinations_with_replacement(['pp','pm','mp','mm','x'],3):
 # Minimise the compressed dimension by optimally choosing the final branch.
 tags=min(sorted(set(itertools.permutations(tags))),key=lambda t:data[t[0]]['e']*data[t[1]]['e']*data[t[2]]['d'])
 try:
  result,V=h.triple_angle(data,M,T,tags);del V;gc.collect();results.append(result)
  with open(out,'w') as f:json.dump({'status':'FLOATING_EXPLORATION_NOT_CERTIFIED','results':results},f,indent=2)
 except Exception as e:
  traceback.print_exc();raise
print('DONE',max(results,key=lambda v:v['theta_estimate']),flush=True)
