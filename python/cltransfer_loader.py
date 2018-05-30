import numpy as np
import pylab as pl
from scipy.interpolate import interp1d
from scipy.integrate import quad

camb_g_l5= np.loadtxt("/Users/darshkodwani/Desktop/CAMB-0.1.6.1/pycamb/grow_trnsf_l5.txt",unpack=True)
camb_d_l5= np.loadtxt("/Users/darshkodwani/Desktop/CAMB-0.1.6.1/pycamb/decay_trnsf_l5.txt",unpack=True)
camb_ks = np.loadtxt("/Users/darshkodwani/Desktop/CAMB-0.1.6.1/pycamb/list_ks.txt",unpack=True)

def load_cltf(fname):

    F=open(fname, 'rb')

    tt_size, l_size, q_size = np.fromfile(F, dtype=np.intc, count=3)

    l = np.fromfile(F, dtype=np.intc, count=l_size)
    q = np.fromfile(F, dtype=np.float64, count=q_size)

    data=np.fromfile(F, dtype=np.float64).reshape(tt_size, l_size, q_size)

    return l, q, data




if __name__=='__main__':
    #
    root='/Users/darshkodwani/Documents/Darsh/Research/Sound_modes/class_dcmode/'

    # scalar #
    folder='output/'
    gname=['decay_scalarcltransfer_ad.dat', 'decay_scalarcltransfer_ad.dat']
	
    folder='output/'
    dname=['decay_scalarcltransfer_addcs.dat', 'decay_scalarcltransfer_addcs.dat']

    # tesnor #
    #folder='output/'
    #fname=['decay_tensor/cltransfer_ad.dat', 'decay_tensor/cltransfer_ten.dat', \
    #       'decay_tensor/cltransfer_addct.dat']



    # for scalar, the first few common outputs are [t2, e, t0, t1, b ...]
    # for tensor, the common outputs are [t2, e, b ]


    gl, gq, gt=load_cltf(root+folder+gname[0])
    dl, dq, dt=load_cltf(root+folder+dname[0])
    print gt.shape
    print 'gl=', dl
    print 'gq=', gq
    

    #tt_idx = -1
    llist=[0, 1, 2, 3, 4, 5, 10, 20, 50]    
    
    for li in llist:
        pl.semilogx(dq, gt[0,li,:]+gt[1,li,:]+gt[2,li,:]+gt[3,li,:]+gt[4,li,:],label= li)   # t0, t1, t2
        #pl.semilogx(llist, dt[0,li,:])
    
    pl.show()
    
    def power(k,kp,ns,As):
        return As*(k/kp)**(ns-1)
    
    def cls(l,k,kp,As,ns):
        T_total= (gt[0,l,:]+gt[2,l,:]+gt[3,l,:])**2
        tr_inter = interp1d(gq,T_total,kind='cubic') 
        cl_int = quad(lambda x: tr_inter(x)*power(x,kp,ns,As)/x, min(gq),max(gq))[0]
        return cl_int
    
    #cls_list=[]
    #for ell in llist:
        #cls_list.append(cls(ell,gq,0.05,10**(-9),0.96))
    
    #cls_interp=interp1d(llist,cls_list, kind='cubic')   
     
    #pl.loglog(llist,cls_list)
    #ls = np.linspace(0,50,50)
    #pl.plot(ls,cls_interp(ls))
    #pl.show()
    
    #pl.semilogx(gq,gt[0,15,:]+gt[2,15,:]+gt[3,15,:],label= 'class grow l=5')
    #pl.semilogx(camb_ks,camb_g_l5,label ='camb grow l=5')
    #pl.semilogx(camb_ks,camb_d_l5,label ='camb decay l=5')
    #pl.semilogx(dq,dt[0,15,:]+dt[2,15,:]+dt[3,15,:],label= 'class decay l=5')
    #pl.xlim([5e-6, 0.5])
    #pl.legend()
    #pl.show()
