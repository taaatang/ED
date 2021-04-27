# %%
import numpy as np
import matplotlib.pyplot as plt
import ExactDiagUtils as ed
#%% eigs
dataDir = "/Volumes/Sandisk/ProjectData/QSL/triangular/periodic/36/18u18d/k0/J2_0.24"
eigDir = dataDir + "/wavefunc"
evals = np.real(np.fromfile(eigDir + "/eval",dtype=np.complex128))
stateNum = 2
print("evals:\n{}".format(evals))
#%% Raman
delta = 0.01
w = np.linspace(0.0, 1.0, 500)
for channel in ["A1", "A2", "E21", "E22"]:
    ramanDir = dataDir + "/raman" + "/" + channel
    chi = np.zeros(len(w))
    for state in range(stateNum):
        spec = ed.read_spec(ramanDir, w, state, delta)
        spec_inv = ed.read_spec(ramanDir, -w, state, delta)
        chi += spec - spec_inv
    chi /= stateNum

    plt.plot(w, chi)
    plt.xlabel(r"$\omega$")
    plt.ylabel(r"{}".format(ed.getRamanLabel(channel)))
    plt.show()
# %%
