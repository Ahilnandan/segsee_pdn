import numpy as np, pdn3
n=3000; I=np.tile(np.array([60.,25.,30.,5.]),(n,1))
for s in range(500,n,600): I[s:s+60]*=2.7
P=dict(A_x=[0.1*x for x in pdn3.SHARE])
def show(tag,spec=None,pred=None,p=P):
    m=pdn3.simulate(I,p,spec,pred=pred); a=pdn3.area(dict(pdn3.DEFAULT,**p),spec)
    print('%-58s max %5.1f p99 %5.1f over %5.1f | buck loss %6.1f%% fsw/reg %5.0f MHz | area %.3f'%(tag,m['max'],m['p99'],m['over_max'],100*m['loss_buck'],m['f_sw_MHz']/4,a['total']))
show('passive 0.1 mm2')
for L in (5e-9,10e-9,20e-9):
  for Vin in (1.3,1.5):
    for dw in (2,4):
      show('buck L=%2.0fnH Vin=%.1f dwell=%d Ron=0.5 Ipk=200mA*share th=15'%(L*1e9,Vin,dw),[dict(buck=dict(L=L,Vin=Vin,Ron=0.5,Ipk=0.2*s,th=0.015,dwell=dw)) for s in pdn3.SHARE])
pred=np.zeros((n,4),np.uint8)
for s in range(500,n,600): pred[s-3:s+20]=1
show('buck 10nH 1.3V + perfect predictor lead 3 (th 40mV)',[dict(buck=dict(L=10e-9,Vin=1.3,Ron=0.5,Ipk=0.2*s,th=0.04,pred=1,dwell=2)) for s in pdn3.SHARE],pred)
