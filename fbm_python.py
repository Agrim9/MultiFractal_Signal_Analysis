from fbm import FBM
import matplotlib.pyplot as plt

f = FBM(n=100, hurst=0.9, length=1, method='daviesharte')

# Generate a fBm realization
fbm_sample = f.fbm()

# Generate a fGn realization
fgn_sample = f.fgn()

# Get the times associated with the fBm
t_values = f.times()
plt.plot(t_values,fbm_sample)
plt.show()