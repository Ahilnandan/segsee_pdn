import json
from multiprocessing import Pool
def job(a):
    wl,cb,h,dm=a
    import pulsim_chargepump as pc, cosim
    pc.CB_F=cb*1e-9; cosim.CB_F=cb*1e-9          # build_circuit reads pc.CB_F; token sensing reads cosim.CB_F
    from sweep import setup
    inst,prods,cons,age,pin=setup(wl)
    r=cosim.cosim(inst,prods,cons,age,pin,ir_mA=25,harvest=h,dummy=dm,THF=(200 if dm else 150))
    r.update(wl=wl,cb=cb,harv=h,dummy=dm or 'none'); return r
if __name__=='__main__':
    jobs=[(wl,cb,h,dm) for wl in ('qsort','dhrystone','mm') for (cb,h,dm) in ((1.0,False,None),(1.0,True,None),(0.5,False,'fill'),(0.5,True,'fill'))]
    with Pool(8) as p, open('cosim2.jsonl','w') as f:
        for r in p.imap_unordered(job,jobs): f.write(json.dumps(r)+'\n'); f.flush(); print(r,flush=True)
