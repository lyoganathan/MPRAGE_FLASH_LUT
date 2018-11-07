import numpy as np

def wang_mprage(T1,TI,ES,TD,theta,N,M0,enc):
    # Based on Optimizing the Magnetization-Prepared Rapid Gradient-Echo (MP-RAGE) Sequence (2014) - Jinghua Wang
    # Also appears in Development and validation of a new MRI simulation technique that can reliably estimate optimal
    # in vivo scanning parameters in a glioblastoma murine model (2018) - Andrea Protti
    #MPRAGE T1
    theta = np.radians(theta)  # numpy trig functions expect radians
    TR = TI + N*ES + TD
    gamma = np.exp(-TI / T1)
    delta = np.exp(-ES / T1)
    phi = np.exp(-TD / T1)
    mu = delta * np.cos(theta)
    rho = np.exp(-TR/T1)
    #M0 = 1 # Shouldn't matter cause gets canceled when we divide
    #Meq is steady state magnetization after several TRs
    #Meq here is Meq/M0 in the paper
    Meq = (1 - phi +
           ((phi*np.cos(theta)*(1-delta)*(1-np.power(mu,N-1)))/(1-mu)) \
           + phi*np.cos(theta)*np.power(mu,N-1) \
           - rho * np.power(np.cos(theta),N)) / \
          (1 + rho * np.power(np.cos(theta),N))
    #M for ith read-out pulse
    M = np.zeros(N)
    for i in range(0,N):
        M[i] = M0 * ((1-delta) * (1-np.power(mu,i)) / (1-mu) \
              + np.power(mu,i) * (1-gamma) \
              - gamma*np.power(mu,i)*Meq)

    #Signal intensity at each phase encoding step:
    S = M*np.sin(theta)

    #Assume signal intensity is dominated by k-space centre
    if enc == 'centric':
        intensity = S[0]
    else:
        intensity = S[round(N/2)]

    return intensity

def bock_mprage(T1,TI,TR,TD,alpha,N,M0,enc): #TR is small TR
    #Based on Optimizing T1-weighted imaging of cortical myelin content at 3.0 T (2013) - Nicholas Bock
    alpha = np.radians(alpha)
    EI = np.exp(-TI/T1)
    ED = np.exp(-TD/T1)
    ETR = np.exp(-TR/T1)
    T1star = (T1 * TR) / (TR - T1 * np.log(np.cos(alpha)))
    ETRstar = np.exp(-TR/T1star)
    tau = N*TR
    Etaustar = np.exp(-tau/T1star)
    #M0 = 1 # Shouldn't matter, get's cancelled out in division
    M0star = M0 * (1 - ETR) / (1 - ETRstar)
    #Magnetization for tissue at steady-state at first echo in FLASH segment...
    M = M0 * (((1-EI) - EI*(1-ED) -
              ED*EI * ((1-ETR)/(1-ETRstar)) * (1-Etaustar))
              / (1+ED*EI*Etaustar))
    #Signal over duration of FLASH segment, t in seconds
    #Multiple t by TR, to get time of each phase encode step
    Ms = np.zeros(N)
    for t in range(0, N):
        Ms[t] = M0star * (1-np.exp(-t*TR/T1star)) + M * np.exp(-t*TR/T1star)

    S = Ms*np.sin(alpha)

    if enc == 'centric':
        intensity = S[0]
    else:
        intensity = S[round(N/2)]

    return intensity

def flash_steady_state(T1,TR,theta,N):
    theta = np.radians(theta)
    M = np.zeros(N)
    for i in range(0,N):
        M[i] = (1-np.exp(-TR/T1)) / (1-np.cos(theta)*np.exp(-TR/T1))

    S = M*np.sin(theta)

    #Linear - centre of k-space is acquired after N/2 phase encoding steps
    #Centric - centre of k-space is acquired first
    #Dosen't seem to matter for this sequence? S[1] vs S[N/2] is the same...

    return S[0]

# Kim's original equations based on Diechmann
# def flash2_kld(M0,T1,rho,TR,alpha,npe1,seg,enc):
#     # Equations from Diechmann et al.NeuroImage expect alpha is in rads
#     # Number of readout lines(i.e.PE Matrix). # pe1 = 256
#
#     pe1 = np.arange(1,npe1/seg +1) # This creates an array and the rest of code is now array.
#
#     #a = len(T1) is one
#     for i in range (0,1):
#
#         # Effective relaxation time.
#         T1star = (1 / T1 - 1. / TR * np.log(np.cos(alpha)))**(-1)
#
#         # Saturation magnetization
#         M0star = rho * (1 - np.exp(-TR / T1)) / (1 - np.exp(-TR / T1star))
#
#         # Equations for magnetization(see text).
#
#         A1 = M0star * (1 - np.exp(-pe1 * TR / T1star)) #pe*TR = duration of acquisition block (tau?)
#         B1 = np.exp(-pe1 * TR / T1star)
#
#         MZ = A1 + B1 * M0 # eq 9 from Diechmann, uses M0(M1) from mprage_kld
#
#         S = MZ * np.sin(alpha) # signal strength, eq 12
#     # S_temp = repmat(S, [seg, 1])
#     # S_temp = S_temp(:)'
#     S_temp = S # S is an array of length pe1
#
#     # frequency space with 0, the center of k-space, at left side.
#
#     if enc == 'centric':
#         S_new = np.concatenate([S_temp[0:len(S_temp):2], np.flipud(S_temp[1:len(S_temp):2])])  # k - space weighting, centric encoding
#     else:
#         S_new = np.fft.fftshift(S_temp)  # k - space weighting, linear encoding
#
#     peakCent = int(np.ceil((len(S_new) + 0.1) / 2))
#
#     IFTS = np.fft.fftshift(abs(np.fft.ifft(S_new))) # Diff from matlab ifft?? Shouldn't matter cause not used in result
#
#     peakheight = IFTS[peakCent]
#
#     #peakwidth = fwhm([1:length(IFTS)], IFTS); % FWHM
#     peakwidth=1
#
#     # peak = peakheight * sign(S_new(1)); # adjusted for k - space signal variations.% this is for a single voxel PSF.
#     # peak = abs(trapz(IFTS)) * sign(S_new(1)); # eq 14 Deichmann
#
#     # THIS LINE WAS USED TO MAKE LOOKUP TABLES
#     peak = S_new[0] # this is assuming the signal is dominated by the signal amplitude at the center of k-space.i.e.for an infinite object.
#     # offpeak = IFTS(peakCent-1); # or +1, shouldn't matter.
#     offpeak=1
#     return peak
#
#
# def flashSS_kld(T1,rho,TR,alpha,pe1):
#     # simulate signal and contrast for steady-state FLASH sequence
#     # Equations from Nishmura % alpha is in rads % MZ is longitudinal magnetization for each readout.
#
#     #a = len(T1) is 1
#
#     for i in range(0,1):
#         MZ = rho * (1 - np.exp(-TR / T1)) / (1 - np.exp(-TR / T1) * np.cos(alpha))
#         S = MZ * np.sin(alpha)
#     #Generate one line of pe1 k - space to calculate the PSF(which will trivally be a % delta function).
#     S1=np.zeros(pe1) # initialize
#     for j in range (0,pe1):
#         S1[j]= S
#     # IFT of signal - don't need for now... also doesn't work
#     #for i in range(0,a):
#         #IFTS[i]=np.fft.fftshift(abs(np.fft.ifft(S1[i])))
#
#     # display('Peak height of IFT of WM, GML, GMH');
#     # peak = max(IFTS(1,:)); orig
#     # peakheight = abs(trapz(IFTS))*sign(S1(1,1)); %2018-02-20 KLD  %eq 14 Deichmann.
#     # For a large tissue region, the observed signal would be the sum of the PSF because
#     # you are adding the contribution from the neighbouring pixels
#     # peakheight = max(IFTS(1,:)); % Nick orig.
#     peakheight = S1[round(pe1 / 2)] # assume linear encoding.
#     # peakwidth = nearestpoint(peakheight / 2, IFTS(1,:)); % FWHM
#     # peak = peakheight / peakwidth; % should help to minimize blurring.
#     peak = peakheight
#
#     # peak2 = max(IFTS(2,:));
#     # peak3 = max(IFTS(3,:));
#     return peak
#
# def mprage_kld(T1,rho,alpha,TR,tau,TI, TD, IE):
#     # 3D-inversion recovery α and N equally-spaced readout RF pulses of flip angle θ and echo spacing tau
#     # TR - Time b/w two inversion recovery pulses (TR = TI + N*ES +TD)
#     # tau - duration of acquisition block
#     # TI - time interval b/w inversion recovery pulse & first RF readout pulse
#     # TD - delay time
#
#     # Directly before an acquisition block, M has the value M1. During the acquisition, M approaches the saturation
#     # value M0star with the time constant T1*
#
#     # T1* is Estimate of T1, usually lower equilibrium magnetization and shorter T1 estimates
#     # T1 is recovery of longitudinal magnetization back to equilibrium
#
#     # The MP-RAGE signal recovery curve is in turn dependent on the apparent relaxation rate 1/T1*
#     # during the acquisition phase and is given by
#
#     T1star = ((1/T1 - 1/TR*np.log(np.cos(alpha)))**(-1))
#
#     # Saturation Magnetization: number of atoms per unit volume * magnetic moment of each atom
#     # Msat(M0star) is maximum volume magnetization total magnetic moment per unit volume
#     M0star = (rho*1-np.exp(-TR/T1)) / (1-np.exp(-TR/T1star))
#
#     # After acquisition block, M = A1 + B1*M1
#     A1 = M0star*(1-np.exp(-tau/T1star))
#     B1 = np.exp(-tau/T1star)
#
#     # There is a delay TD, after which M3=A2+B2*M2
#     A2 = rho * (1 - np.exp(-TD / T1))
#     B2 = np.exp(-TD / T1)
#
#     A3 = rho * (1 - IE * np.exp(-TI / T1))
#     B3 = -IE * np.exp(-TI/ T1)
#
#     A = A3 + A2 * B3 + A1 * B2 * B3
#     B = B1 * B2 * B3
#
#     M1 = A/(1-B) # Used in flash2_kld
#
#     # S = M1. * sin(alpha);
#
#     return M1
