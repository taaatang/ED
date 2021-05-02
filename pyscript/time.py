#%% import
import numpy as np
import matplotlib.pyplot as plt
#%% data directory
dataDir = "/Volumes/Sandisk/ProjectData/photodoping/test/square/periodic/2x2/3u2d/k1/Udp/pump/pol_z/i_1.00/w_1.70/phi_0.25/"
orbs = ["d", "px", "py", "pzu", "pzd"]
orbNum = len(orbs)
num = []
for i in range(orbNum):
    num.append(np.fromfile(dataDir + "orb" + str(i), dtype=np.float64))
num = np.array(num)
t = 0.01 * np.array(list(range(len(num[0,:]))))
#%%
for i,orb in enumerate(orbs):
    plt.plot(t,num[i,:] - num[i,0], label = orb)
plt.legend()
# %%
