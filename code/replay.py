#!/usr/bin/env python3
"""Rebuild and verify the v0.3 proof. This never accepts stored scalar receipts.

Default: full deterministic replay. --quick omits the 50 main-sector tests, so
its successful completion is explicitly NOT a full kagome proof replay.
"""
import argparse,os,sys,subprocess,json,time,platform,hashlib,shutil
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--quick',action='store_true');a=p.parse_args()
root=Path(__file__).resolve().parent.parent;code=root/'code';work=root/'.work';logs=root/'logs'/'replay';evidence=root/'evidence'/'new'
for d in [work,logs,evidence]:d.mkdir(parents=True,exist_ok=True)
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',AKLT_WORK=str(work),AKLT_EVIDENCE=str(evidence),PYTHONHASHSEED='0')
# Remove authoritative receipts so a stale result cannot stand in for a stage.
expected=['block_gap_validated','star_inertia','bowtie_certified','structural_checks','validator_tests']
if not a.quick:expected+=['validated_all']
for name in expected:
 f=evidence/(name+'.json')
 if f.exists():f.unlink()
t=time.time();steps=[]
def call(name,cmd):
 print('REPLAY',name,flush=True);t0=time.time()
 with open(logs/(name+'.log'),'w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env,cwd=root)
 steps.append({'step':name,'returncode':r.returncode,'elapsed_seconds':time.time()-t0})
 (evidence/'replay_receipt.json').write_text(json.dumps({'status':'IN_PROGRESS' if r.returncode==0 else 'REJECT','full':not a.quick,'steps':steps},indent=2))
 if r.returncode:raise SystemExit('REJECT: '+name+'; see '+str(logs/(name+'.log')))
def py(name):return [sys.executable,'-u',str(code/'certified'/name)]
call('compile_exact_inertia',['g++','-O3','-std=c++17',str(code/'certified'/'inertia.cpp'),'-lgmpxx','-lgmp','-o',str(work/'inertia')])
# All exact boundary and projector inputs are reconstructed, never copied from receipts.
call('rebuild_exact_tensors',py('rational_hub.py'))
if not a.quick:call('all_50_angles',py('run_validated_all.py'))
call('block_gap',py('block_gap.py'))
call('star_inertia',py('star_inertia.py'))
call('bowtie_exact',py('bowtie_certified.py'))
call('structure_and_geometry',py('structural_checks.py'))
call('validator_adversarial_tests',py('validator_tests.py'))
call('validator_adversarial_tests_optimized',[sys.executable,'-O','-u',str(code/'certified'/'validator_tests.py')])
for name in ['exact_core','exact_star','exact_triangle','exact_cross','recoupling']:
 call(name,[sys.executable,'-u',str(code/'exact_legacy'/(name+'.py'))])
 call(name+'_optimized',[sys.executable,'-O','-u',str(code/'exact_legacy'/(name+'.py'))])
call('environment',[sys.executable,'-c','import sys,platform,numpy,scipy,sympy; print(sys.version); print(platform.platform()); print("NumPy",numpy.__version__,"SciPy",scipy.__version__,"SymPy",sympy.__version__); numpy.show_config()'])
for name in expected:
 obj=json.loads((evidence/(name+'.json')).read_text())
 if obj.get('status')!='PASS':raise SystemExit('REJECT: invalid authoritative receipt '+name)
 if name=='block_gap_validated' and (len(obj['results'])!=4 or any(r['status']!='PASS' for r in obj['results'])):raise SystemExit('REJECT: block coverage')
 if name=='validated_all' and (len(obj['results'])!=50 or any(r['status']!='PASS' for r in obj['results'])):raise SystemExit('REJECT: angle coverage')
receipt={'status':'PASS_FULL' if not a.quick else 'PASS_QUICK_NOT_FULL_PROOF','full':not a.quick,'steps':steps,'elapsed_seconds':time.time()-t,'python':sys.version,'platform':platform.platform(),'thread_environment':{k:env[k] for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','PYTHONHASHSEED']},'source_sha256':{str(x.relative_to(root)):hashlib.sha256(x.read_bytes()).hexdigest() for x in sorted(code.rglob('*')) if x.suffix in ['.py','.cpp']}}
(evidence/'replay_receipt.json').write_text(json.dumps(receipt,indent=2));print(receipt['status'],time.time()-t,flush=True)
