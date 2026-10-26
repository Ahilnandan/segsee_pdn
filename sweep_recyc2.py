import json
from multiprocessing import Pool
from sweep_recyc import job
if __name__=='__main__':
    jobs=[(wl,'TD',Qb,Ir,150,0.8,'harvest') for wl in ('mm','qsort','dhrystone') for Qb in (500,1000) for Ir in (0,10,25)]
    with Pool(9) as p, open('sweep_recyc2.jsonl','w') as f:
        for r in p.imap_unordered(job,jobs): f.write(json.dumps(r)+'\n'); f.flush()
