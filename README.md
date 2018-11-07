## Using the code

Using a T1W-MPRAGE, PDW-FLASH and B1+ map.

Generating lookup table:
```
# Make sure you have numpy
import numpy as np
# For plotting:
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# pandas to read and write to csv
import pandas as pd
from MPRAGE_FLASH_lookup import generate_lookup

# Sequence parameters: TI, ES, TD, alpha, pe1, enc_scheme,ES_PD,alpha_PD,pe1_PD
# This takes about 5 minutes, can be made faster or slower depending on how
# dense you want to sample the function
T1_lookup = generate_lookup(1,0.007916,1.1,12,100,'centric',0.007916,4,180)

# Plotting lookup table:
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(T1_lookup[:, 1], T1_lookup[:, 2], T1_lookup[:, 0], cmap=plt.cm.jet, linewidth=0.2, antialiased=True)
ax.set_xlabel('B1')
ax.set_ylabel('Ratio')
ax.set_zlabel('T1')
ax.view_init(azim=210)
plt.show()

# Also can save to CSV if you don't want to lose it
lookup_df = pd.DataFrame(T1_lookup, columns=["T1","B1err","Ratio"])
lookup_df.to_csv("T1_Lookup_Table.csv")

```
![](imgs/Lookup_Table.png)

To make the sampling finer or coarser, you can go into MPRAGE_FLASH_lookup.py and change the step-size or range of T1 and B1 values.
Now we can look up ratio image intensities and apply B1+ correction them:
If you have freesurfer run, and map your data onto the surface folliwing HCP pipelines you can load in your gifti files, run it through the look up Table, and export it as R1 maps:

![](imgs/Ratio_Image.JPG) ![](imgs/B1_map.JPG)


We can use this surface data as input into the lookup table:

```
# Import gifti data using nibabel
import nibabel as nib
import os

# Loop through and get all ratio surface data
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

#Loop again for B1 maps.
path=r'C:\Users\laagi_000\Documents\Laagi\New folder\Sample'
input_dirs=next(os.walk(path))[1]
lh_b1s = np.zeros(shape=(len(input_dirs),32492))
rh_b1s = np.zeros(shape=(len(input_dirs),32492))
for x in range(0,len(input_dirs)):
    subj = input_dirs[x]
    lh_path = '{0}\\{1}\\32k\\lh.B1.32k_fs_LR.func.gii'.format(path,subj)
    rh_path = '{0}\\{1}\\32k\\rh.B1.32k_fs_LR.func.gii'.format(path, subj)
    lh_b1 = nib.load(lh_path)
    rh_b1 = nib.load(rh_path)
    lh_b1s[x,:] = lh_b1.darrays[0].data/150
    rh_b1s[x,:] = lh_b1.darrays[0].data/150

```

Now we can use our lookup table and find values of T1 for every R1 and B1+.

```
from scipy import interpolate
T1_Lookup = pd.read_csv("T1_Lookup_Table.csv", usecols=["B1err","ratio","T1"])
#Interpolate with griddata
lh_t1 = interpolate.griddata(np.array([T1_Lookup["B1err"],T1_Lookup["Ratio"]]).T,
                          T1_Lookup["T1"],
                          (lh_b1s.flatten()/160,lh_ratios.flatten()))

rh_t1 = interpolate.griddata(np.array([T1_Lookup["B1err"],T1_Lookup["Ratio"]]).T,
                          T1_Lookup["T1"],
                          (rh_b1s.flatten()/160,rh_ratios.flatten()))

#Reshape back into original:

lh_t1s = lh_r1.reshape(10,32492)
rh_t1s = rh_r1.reshape(10,32492)

lh_r1s = 1/lh_t1s
rh_t1s = 1/lh_t1s

#Mean and export to gifti, using template of imported gifti: Basically Changing a random imported Gifti image's data with the mean R1 data

lh_ratio.darrays[0].data = 1/np.mean(lh_r1s,0,dtype='float32')
rh_ratio.darrays[0].data = 1/np.mean(rh_r1s,0,dtype='float32')

nib.save(lh_ratio,'{}\\lh_R1_mean.func.gii'.format(path))
nib.save(rh_ratio,'{}\\rh_R1_mean.func.gii'.format(path))
```
![](imgs/T1_Map.JPG) ![](imgs/R1_Map.JPG) ![](imgs/Ratio_Image.JPG)

All brain images were made by exporting to gifti and viewing in connectome-workebnch.

You can also add a parameter for inversion efficiency (the theta in wang_mprage).

###Deciphering sequence parameters (Siemens & GE):
For example GE's opti setting for their BRAVO sequence combines TI and TD. They have a pos_start_ir which is TI, and opti minus that gives TD.
TR in GE protocol and DICOM refers to echo spacing, whereas on Siemens, TR is the repetition time.

##Quantitative MRI
Map values of T1, a property of tissue. Several methods exist to do T1 mapping across the whole brain.
Popular ones include. Paper. We can also use variable flip angle method. Generally, two SPGR/FLASH images are collected.
Keep everything exactly the same, but change flip angle. One T1W image and PDW image. We can do a similar method, using
MPRAGE and FLASH images, based on this paper: Bock et al.

The principle is, the signal intensity can be worked out mathematically as a function of imaging parameters. For example the code uses MPRAGE equations published from Wang et al. and Protti et al. The Bock et al. equations are also in the code and they give the same result. FLASH also known as SPGR equations have also been worked out.

# B1+ & Inversion efficiency

##Imaging parameters:
TI - time between inversion pulse and beginning of acquisition block
ES - time between successive excitation pulses
TD - time delay from end of acquisition to next inversion
flip angle -
Phase encoding steps - Number of phase encodes in the inner loop
Phase encoding scheme - We assume two possibilities, centric (k-space centre acquired first) or linear (k-space centre acquired after N/2 steps)

If you have an adibiatic pulse, your sequence should be insensitive to B1+ inhomogenity.
These signal equations are idealized, many other factors play a role such as receiver gain, bandwidth,
We hope the division takes care of many things, but there things it doesn't take care of like imperfect spoiling. We can also incorporate an inversion efficiency term, however we don't know how to map that across the brain.


Using a lookup table, we assume many things. We consider k-space centre to be intensity values. K-space centre does determine most of the contrast, but not for high spatial frequencies.
