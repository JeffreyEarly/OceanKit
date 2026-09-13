"""Paired fresh-process screens of the actual native retained advection operator."""
from pathlib import Path
import argparse,hashlib,json,subprocess,statistics,time,math
import numpy as np
p=argparse.ArgumentParser();p.add_argument('archive',type=Path);p.add_argument('--blocks',type=int,default=4);args=p.parse_args();a=args.archive.resolve();out=a/'stage-screen';out.mkdir(exist_ok=False)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
roles=['baseline','memory-only','simd'];binaries={r:a/r/'NativeStage' for r in roles};frozen={str(v):digest(v) for v in binaries.values()}
for attempt in range(6):
 time.sleep(10);ps=subprocess.check_output(['ps','-axo','pid=,pcpu=,comm='],text=True);(out/f'preflight-{attempt}.txt').write_text(ps)
 if not any(float(s.split(None,2)[1])>10 for s in ps.splitlines()):break
else:raise RuntimeError('Host busy; no process altered')
save(out/'protocol.json',{'kind':'exploratory native-stage screen, excludes MM/assembly/integration','blocks':args.blocks,'perProcessWarmups':2,'perProcessSamples':9,'binarySHA256':frozen,'modeSHA256':{str(a/f'modes-{nz}.txt'):digest(a/f'modes-{nz}.txt') for nz in [28,129]},'order':'rotate three roles per block, reverse every other block'})
rows=[]
for nz in [28,129]:
 for targets in [3,4]:
  for block in range(args.blocks):
   order=roles[block%3:]+roles[:block%3]
   if block%2:order=list(reversed(order))
   reports={};paths={}
   for role in order:
    prefix=out/f'nz{nz}-t{targets}-block{block}-{role}';paths[role]=prefix.with_suffix('.bin')
    result=subprocess.run([str(binaries[role]),str(a/f'modes-{nz}.txt'),str(targets),str(prefix)],check=True,capture_output=True,text=True)
    (prefix.with_suffix('.stderr')).write_text(result.stderr);r=json.loads(result.stdout);reports[role]=r;save(prefix.with_suffix('.json'),r)
   baseline=np.fromfile(paths['baseline'],dtype=np.float64)
   comparisons={}
   for role in roles[1:]:
    candidate=np.fromfile(paths[role],dtype=np.float64);assert baseline.shape==candidate.shape
    # Scale each physical field and each split target component separately.
    m=int((a/f'modes-{nz}.txt').read_text().split()[3]);sizes=[256*256*nz]*4+[nz*m]*targets*2;start=0;errors=[]
    for size in sizes:
     b=baseline[start:start+size];c=candidate[start:start+size];error=float(np.max(np.abs(b-c))/max(1.,float(np.max(np.abs(b)))));assert error<=2e-12,(nz,targets,role,error);errors.append(error);start+=size
    assert start==len(baseline)
    for key in ['columnInverses','rowInverses','reusedColumns']:assert reports[role][key]==reports['baseline'][key]
    comparisons[role]={'bitwiseEqual':bool(np.array_equal(baseline.view(np.uint64),candidate.view(np.uint64))),'maximumScaledError':max(errors),'workspaceBytes':reports[role]['workspaceBytes'],'ratioToBaseline':statistics.median(reports[role]['seconds'])/statistics.median(reports['baseline']['seconds'])}
   row={'nz':nz,'targets':targets,'block':block,'medianSeconds':{r:statistics.median(v['seconds']) for r,v in reports.items()},'comparisons':comparisons,'outputSHA256':{r:digest(v) for r,v in paths.items()}};rows.append(row);save(out/'pairs.json',rows);print(row,flush=True)
   for path in paths.values():path.unlink()
summary=[]
for nz in [28,129]:
 for targets in [3,4]:
  group=[r for r in rows if r['nz']==nz and r['targets']==targets]
  summary.append({'nz':nz,'targets':targets,'memoryToBaseline':math.exp(statistics.mean(math.log(r['medianSeconds']['memory-only']/r['medianSeconds']['baseline']) for r in group)),'simdToMemory':math.exp(statistics.mean(math.log(r['medianSeconds']['simd']/r['medianSeconds']['memory-only']) for r in group)),'simdToBaseline':math.exp(statistics.mean(math.log(r['medianSeconds']['simd']/r['medianSeconds']['baseline']) for r in group))})
assert all(digest(Path(p))==h for p,h in frozen.items());save(out/'summary.json',summary);print(summary)
