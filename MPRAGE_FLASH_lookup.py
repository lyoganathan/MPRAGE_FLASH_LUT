import numpy as np
from equations import bock_mprage, wang_mprage, flash_steady_state
# Uses MPRAGE and FLASH signal equations from equations.py at different values of T1 and B1+ error to calculate what the
# ratio intensity of MPRAGE/FLASH (T1W/PDW) should be

########################################################################################################################
# #This block tests that the three mprage equations give the same answer, and also visualize the effect of B1
# #Also visualize HC and LC curves
# #For values of T1, how does the ratio + B1err curve look like:
# T1=np.arange(0.8,1.5,0.1)
# B1err=np.arange(0.5,1.5,0.01)
# ratio_wang = np.zeros([len(T1),len(B1err)])
# ratio_wang_HC = np.zeros([len(T1),len(B1err)])
# ratio_bock = np.zeros([len(T1),len(B1err)])
# ratio_kim = np.zeros([len(T1),len(B1err)])
# for i in range(0,len(T1)):
#     for j in range(0,len(B1err)):
#         #T1,TI,ES,TD,flip angle,phase encoding steps,M0,phase encoding scheme
#         #LC_wang = wang_mprage(T1[i],0,0.0071,0.4895,9*B1err[j],256,1,'linear') # Siemens
#         LC_wang = wang_mprage(T1[i], 0.450, 0.007964, 0, 12 * B1err[j], 180, 1, 'centric')  # GE
#         #LC_bock = bock_mprage(T1[i],0,0.0071,0.4895,9*B1err[j],256,1,'linear')
#         #HC_wang = wang_mprage(T1[i],0.37105,0.0061,1.82105,12*B1err[j],240,1,'linear') # Siemens
#         HC_wang = wang_mprage(T1[i], 1, 0.007916, 1.1, 12 * B1err[j], 100, 1, 'centric')  # GE
#         #new_flash = flash_steady_state(T1[i],0.0061,4*B1err[j],240) # Siemens
#         new_flash = flash_steady_state(T1[i], 0.007908, 4 * B1err[j], 189)  # GE
#         ratio_wang[i,j] = LC_wang/new_flash
#         #ratio_bock[i, j] = LC_bock/new_flash
#         ratio_wang_HC[i, j] = HC_wang / new_flash
#
# labels = []
# for i in range(0,len(T1)):
#     labels.append("T1={}".format(round(T1[i],2)))
#
# plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(T1)))))
# plt.plot(B1err,ratio_wang.T) #wang, bock & kim should give same thing
# plt.xlabel('B1 error(1 would be no deviation from prescribed flip angle)')
# plt.ylabel('Ratio(T1WHC/PDW)')
# plt.title("Effect of B1 Error on Low Contrast Ratio Image")
# plt.legend(labels, loc='upper left')
# plt.show() #You can see lower T1 (Higher R1) gives higher ratio
########################################################################################################################
def generate_lookup(TI, ES, TD, alpha, pe1, enc_scheme,ES_PD,alpha_PD,pe1_PD):
    minT1 = 300 # milliseconds
    maxT1 = 4010
    T1 = np.arange(minT1, maxT1, 10) / 1000 # Seconds
    B1err = np.arange(0.1, 2.05, 0.05) # B1+ deviations

    column_list = ["T1","B1err","Ratio"]
    T1_lookup = np.zeros(shape=(len(T1)*len(B1err),len(column_list)))
    index = 0
    for i in range(0, len(T1)):
        for j in range(0, len(B1err)):
            # T1, TI, ES, TD, flip angle, PE1, M0, phase encode scheme
            mprage = wang_mprage(T1[i],TI,ES,TD,alpha*B1err[j],pe1,1,enc_scheme) # bock_mprage gives same answer
            flash = flash_steady_state(T1[i],ES_PD,alpha_PD*B1err[j],pe1_PD)
            T1_lookup[index, 0] = T1[i]
            T1_lookup[index, 1] = B1err[j]
            T1_lookup[index, 2] = mprage/flash
            index = index + 1

    return T1_lookup