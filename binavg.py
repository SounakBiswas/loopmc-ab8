import numpy as np
V=-2.0

c=0
fname="./outfiles/binfile_mon_v%.2f.dat"%(V)
data=np.loadtxt(fname)
data=data[:,2:]
lines=np.shape(data)[0]
avg=np.sum(data,axis=0)/(1.0*lines)
err=np.sum(data**2,axis=0)/(1.0*lines)
err=np.sqrt((err-avg**2)/(1.0*lines))
np.savetxt("c%dtmdensities_v%.2f.dat"%(c,V),avg)
#np.savetxt("c%dtmderrs.dat"%(c),err)


fname="./outfiles/binfile_dim_v%.2f.dat"%(V)
data=np.loadtxt(fname)
data=data[:,2:]
lines=np.shape(data)[0]
avg=np.sum(data,axis=0)/(1.0*lines)
err=np.sum(data**2,axis=0)/(1.0*lines)
err=np.sqrt((err-avg**2)/(1.0*lines))
np.savetxt("c%dtddensities_v%.2f.dat"%(c,V),avg)
#np.savetxt("c%dtdderrs.dat"%(c),err)
