import certify_hub as h,itertools,json,time,os,gc
from fractions import Fraction as F
from validated import require
from paths import EVIDENCE
out=EVIDENCE/'validated_all.json'
data,M,em=h.load();results=[];start=time.time()
for tags in itertools.combinations_with_replacement(['pp','pm','mp','mm','x'],2):
 tags=min(set(itertools.permutations(tags)),key=lambda t:(data[t[0]]['e']*data[t[1]]['d'],t))
 results.append(h.pair(data,M,em,tags));gc.collect()
 with open(out,'w') as f:json.dump({'status':'IN_PROGRESS','results':results},f,indent=2)
 require(results[-1]['status']=='PASS','pair certificate rejected')
for tags in itertools.combinations_with_replacement(['pp','pm','mp','mm','x'],3):
 tags=min(set(itertools.permutations(tags)),key=lambda t:(data[t[0]]['e']*data[t[1]]['e']*data[t[2]]['d'],t))
 results.append(h.triple(data,M,em,tags));gc.collect()
 with open(out,'w') as f:json.dump({'status':'IN_PROGRESS','results':results},f,indent=2)
 require(results[-1]['status']=='PASS','triple certificate rejected')
require(len(results)==50,'coverage count')
with open(out,'w') as f:json.dump({'status':'PASS','pair_sectors':15,'triple_sectors':35,'results':results,'elapsed_seconds':time.time()-start,'scope':'sector-dependent two-plus-one gluing; no global common grouping claimed'},f,indent=2)
print('ALL 50 VALIDATED',time.time()-start,flush=True)
