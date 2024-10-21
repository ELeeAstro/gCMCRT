import numpy as np

#wls = 1.12
#wle = 1.17

#wls = 1.21
#wle = 1.32

#wls = 1.35
#wle = 1.43

#wls = 1.54
#wle = 1.60

wls = 1.62
wle = 1.69

nwl = 100

wl = np.linspace(wls, wle, nwl)

T = np.ones(nwl)

#print(wl)
#print(T)

for i in range(nwl):
    print(wl[i]*1e4,T[i])
