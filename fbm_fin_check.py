from fbm import FBM
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.interpolation import shift
from scipy.linalg import toeplitz,inv,svd
from scipy.optimize import curve_fit
import pywt
import pdb
def Qf(x, A, B,C): # this is your 'straight line' y=f(x)
    return A*x**2 + B*x + C

def get_wavelet_var(wavelet,H=0.3,samp_size=1024,plot_sil=1):
    f = FBM(n=samp_size, hurst=H, length=0.5, method='daviesharte')
    fbm_sample = f.fbm()
    t_values = f.times()
    cA, cD = pywt.dwt(fbm_sample, wavelet)
    
    M=len(cD)
    R=np.zeros(M+1)

    for k in range(M+1):
        cD_shift=shift(cD, k, cval=0)
        R[k]=np.sum(cD*np.conjugate(cD_shift))

    print(M+1)
    print(np.argmin(R))
    plt.plot(R)
    plt.show()
    AutoCorr_mat=toeplitz(R[0:M],R[0:M])+2.5
    #pdb.set_trace()

    mean_corr=(np.sum(AutoCorr_mat)-M*AutoCorr_mat[0][0])/(M**2-M)
    #print(mean_corr)
    max_corr=np.max(R[1:])
    #pdb.set_trace()
    U,S,V=svd(AutoCorr_mat)
    if(plot_sil!=1):
    	plt.plot(np.log10(np.arange(M)+1),np.log10(S))

    return S[0],1+(M-1)*mean_corr
print(get_wavelet_var('haar',0.3,1024,0))
#plt.show()
