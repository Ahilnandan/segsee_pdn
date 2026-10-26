import numpy as np, pdn3
n=3000; I=np.tile(np.array([60.,25.,30.,5.]),(n,1))
for s in range(500,n,600): I[s:s+60]*=2.7
P=dict(A_x=[0.1*x for x in pdn3.SHARE])
def show(tag,spec=None,pred=None,p=P):
    m=pdn3.simulate(I,p,spec,pred=pred); a=pdn3.area(dict(pdn3.DEFAULT,**p),spec)
    print('%-52s max %5.1f p99 %5.1f over %5.1f mV | buck loss %5.1f%% fsw %5.0f MHz | ad %5.1f%% | area %.3f mm2'%(tag,m['max'],m['p99'],m['over_max'],100*m['loss_buck'],m['f_sw_MHz'],100*m['loss_ad'],a['total']))
show('passive 0.1 mm2 decap')
show('passive 0.4 mm2 decap',p=dict(A_x=[0.4*x for x in pdn3.SHARE]))
for L,Ron,Ipk,th in ((2e-9,1.0,0.08,0.015),(2e-9,0.5,0.15,0.015),(2e-9,0.5,0.15,0.025),(5e-9,0.5,0.15,0.015)):
    show('buck L=%.0fnH Ron=%.1f Ipk=%.0fmA th=%.0fmV'%(L*1e9,Ron,Ipk*1e3,th*1e3),[dict(buck=dict(L=L,Ron=Ron,Ipk=Ipk,th=th))]*4)
for Ca,G in ((1e-9,4),(2e-9,4),(2e-9,8)):
    show('active decap Ca=%.0fnF G=%d'%(Ca*1e9,G),[dict(ad=dict(Ca=Ca,G=G))]*4)
pred=np.zeros((n,4),np.uint8)
for s in range(500,n,600): pred[s-3:s+20]=1
show('buck 2nH 0.5ohm 150mA + perfect predictor lead 3',[dict(buck=dict(L=2e-9,Ron=0.5,Ipk=0.15,th=0.015,pred=1))]*4,pred)
