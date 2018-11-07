from MPRAGE_FLASH_lookup import generate_lookup
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

#Sequence parameters for a GE BRAVO sequence (basically MPRAGE) and GE SPGR sequence (basically FLASH):
#T1W-MPRAGE: TI=1, ES=0.007964, TD=1.1, flip angle = 12, pe1 = 100, centric encoding
#PDW-FLASH: ES=0.007916, flip angle = 4, pe1=180

T1_lookup = generate_lookup(1,0.007964,1.1,12,100,'centric',0.007908,4,180)

#To pandas df and CSV so we don't lose it.
lookup_df = pd.DataFrame(T1_lookup, columns=["T1","B1err","Ratio"])
lookup_df.to_csv("T1_Lookup_Table.csv")

# Create graphs:
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(T1_lookup[:, 1], T1_lookup[:, 2], T1_lookup[:, 0], cmap=plt.cm.jet, linewidth=0.2, antialiased=True)
ax.set_xlabel('B1')
ax.set_ylabel('Ratio')
ax.set_zlabel('T1')
ax.view_init(azim=210)
plt.show()

#How to read from CSV:
T1_Lookup = pd.read_csv("T1_Lookup_Table.csv", usecols=["T1","B1err","Ratio"])

#Load in Ratios with nibabel and for loops:

import nibabel as nib
import os

path=r'C:\Users\laagi_000\Documents\Laagi\New folder\Sample'
input_dirs=next(os.walk(path))[1]
lh_ratios = np.zeros(shape=(len(input_dirs),32492))
rh_ratios = np.zeros(shape=(len(input_dirs),32492))
for x in range(0,len(input_dirs)):
    subj = input_dirs[x]
    lh_path = '{0}\\{1}\\32k\\lh.MyelinMap.32k_fs_LR.func.gii'.format(path,subj)
    rh_path = '{0}\\{1}\\32k\\rh.MyelinMap.32k_fs_LR.func.gii'.format(path, subj)
    lh_ratio = nib.load(lh_path)
    rh_ratio = nib.load(rh_path)
    lh_ratios[x,:] = lh_ratio.darrays[0].data
    rh_ratios[x,:] = rh_ratio.darrays[0].data

#Change a random imported Gifti image's data with the mean data
lh_ratio.darrays[0].data = np.mean(lh_ratios,0,dtype='float32')
rh_ratio.darrays[0].data = np.mean(rh_ratios,0,dtype='float32')

#Save the gifti image
nib.save(lh_ratio,'{}\\lh_ratio_mean.func.gii'.format(path))
nib.save(rh_ratio,'{}\\rh_ratio_mean.func.gii'.format(path))


#Interpolate with griddata
lh_r1 = interpolate.griddata(np.array([T1_Lookup["B1err"],T1_Lookup["Ratio"]]).T,
                          T1_Lookup["T1"],
                          (lh_b1s.flatten()/160,lh_ratios.flatten()))

rh_r1 = interpolate.griddata(np.array([T1_Lookup["B1err"],T1_Lookup["Ratio"]]).T,
                          T1_Lookup["T1"],
                          (rh_b1s.flatten()/160,rh_ratios.flatten()))

#Reshape back into original:

lh_r1s = lh_r1.reshape(10,32492)
rh_r1s = rh_r1.reshape(10,32492)

lh_ratio.darrays[0].data = 1/np.mean(lh_r1s,0,dtype='float32')
rh_ratio.darrays[0].data = 1/np.mean(rh_r1s,0,dtype='float32')

nib.save(lh_ratio,'{}\\lh_R1_mean.func.gii'.format(path))
nib.save(rh_ratio,'{}\\rh_R1_mean.func.gii'.format(path))