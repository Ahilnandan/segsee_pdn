import json
from multiprocessing import Pool
from sweep_dummy import job
if __name__=='__main__':
    jobs=[(wl,pol,TH,Qb,Ir,'bridge') for wl in ('mm','qsort','dhrystone') for TH in (60,100) for Qb,Ir in ((0,0),(250,25)) for pol in ('A','TD')]
    with Pool(12) as p, open('sweep_bridge.jsonl','w') as f:
        for r in p.imap_unordered(job,jobs):
            f.write(json.dumps(r)+'\n'); f.flush(); print(r['wl'],r['pol'],r['TH'],r['Qb'],'done',flush=True)
