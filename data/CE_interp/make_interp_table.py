import numpy as np


fname = 'FastChem_1x_cond.txt'
#fout = 'interp_FastChem_10x_cond.txt'
fout = 'gCMCRT_RCE_1x_cond.txt'

#fname = 'FastChem_1x.txt'
#fout = 'mini_chem_IC_FastChem_1x.txt'

# Read first line to get species names
with open(fname) as f:
    keyword = np.array(f.readline().split())

# Now read data using read
data = np.loadtxt(fname)

Pg = data[:,0]
Tg = data[:,1]
mu = data[:,4]

nl = len(Pg)

print(nl)

# Get a list of names of molecules as in the file
#f_list = ['H2','H','H1-','e-','He','H2O1','C1O1','C1O2','C1H4','O1Ti1','O1V1','Fe','Fe1+','Fe1H1','O1Si1','Na','K','H3N1','H2S1','C1H1N1_1','Cl1H1','H3P1','H1O1','C2H2','F1H1','H1S1','V','V1+','Ti','Ti1+']
#f_list = ['H1O1','H2','H2O1','H','C1O1','C1O2','O','C1H4','C2H2','H3N1','N2','C1H1N1_1','He']
f_list = ['H2','He','H2O1','C1H4','C1O1','C1O2','H3N1','N2','Na','K']


# The normal text for each species for output
#f_names = ['H2','H','H-','e-','He','H2O','CO','CO2','CH4','TiO','VO','Fe','Fe+','FeH','SiO','Na','K','NH3','H2S','HCN','HCl','PH3','OH','C2H2','HF','SH','V','V+','Ti','Ti+']
#f_names = ['OH','H2','H2O','H','CO','CO2','O','CH4','C2H2','NH3','N2','HCN','He']
f_names = ['H2','He','H2O','CH4','CO','CO2','NH3','N2','Na','K']

nm = len(f_names)

VMR = np.zeros((nl,nm))
for m in range(nm):
    ind = np.where(keyword == f_list[m])[0][0]
    print(ind, keyword[ind], f_list[m], f_names[m])
    VMR[:,m] = data[:,ind]

# Prepare interpolation table
P = np.unique(Pg)
T = np.unique(Tg)

nP = len(P)
nT = len(T)

output = open(fout,'w')
output.write(str(nT) + ' ' + str(nP) + ' ' + str(nl) + ' ' + str(nm) + '\n')
output.write(' '.join(map(str, f_names)) + '\n')
output.write(' '.join(map(str, T)) + '\n')
output.write(' '.join(map(str, P)) + '\n')

for i in range(nl):
    output.write(str(mu[i]) + ' ' +  ' '.join(map(str, VMR[i,:])) + '\n')
output.close()
