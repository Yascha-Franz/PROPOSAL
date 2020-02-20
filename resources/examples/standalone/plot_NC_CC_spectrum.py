import pyPROPOSAL as pp
import numpy as np
import matplotlib.pyplot as plt

NuTauDef = pp.particle.NuTauDef.get()
TauDef = pp.particle.TauMinusDef.get()

Direction = pp.Vector3D(0,0,-1)
Position = pp.Vector3D(0,0,6.374134e12)
statistic = 10**(3)
Energy_bins = 10**(np.linspace(5, 12,100))

nutauprop = pp.Propagator(NuTauDef, "../../config_rock_neutrino.json")
#tauprop = pp.Propagator(TauDef, "../../config_rock_lepton.json")
print("Propagator erstellt")
p_conditon = pp.particle.DynamicData(0)
def propagate(prop, direction, pos, e):
    p_conditon.direction = direction
    p_conditon.position = pos
    p_conditon.energy = e
    p_conditon.propagated_distance = 0
    p_conditon.time = 0

    secondarys = prop.propagate(p_conditon)

    return secondarys


def E(e):
    return np.max(np.append(Energy_bins[Energy_bins<e],0))

y_CC = np.array([[0,0]])
y_NC = np.array([[0,0]])
for i in range(0,np.size(Energy_bins)):
    for j in range(0, statistic):
        print(i, '  ', np.log10(j+1))
        secondarys = propagate(nutauprop, Direction, Position, Energy_bins[i])
        for secondary in secondarys.particles:
            if secondary.id == int(pp.particle.Interaction_Id.WeakInt):
                y_CC = np.append(y_CC, [[1-secondary.energy/secondary.parent_particle_energy, E(secondary.parent_particle_energy)]], axis = 0)
            if secondary.id == int(pp.particle.Interaction_Id.WeakIntNC):
                y_NC = np.append(y_NC, [[1-secondary.energy/secondary.parent_particle_energy, E(secondary.parent_particle_energy)]], axis = 0)

def mean(Y):
    mean_ = np.array([])
    for e in Energy_bins:
        m = 0
        if np.size(Y[Y[:,1]==e][:,0])>0:
            m = np.mean(Y[Y[:,1]==e][:,0])
        mean_ = np.append(mean_, m)
    return mean_

mean_CC = mean(y_CC)
mean_NC = mean(y_NC)

plt.subplot(2,1,1)

plt.step(Energy_bins/1e3, mean_CC)
plt.title('CC')
plt.xlabel('E [GeV]')
plt.xlim(np.min(Energy_bins), np.max(Energy_bins))
plt.ylabel('y_expected')
plt.xscale('log')

plt.subplot(2,1,2)

plt.step(Energy_bins/1e3, mean_NC)
plt.title('NC')
plt.xlabel('E [GeV]')
plt.xlim(np.min(Energy_bins), np.max(Energy_bins))
plt.ylabel('y_expected')
plt.xscale('log')

plt.tight_layout()
plt.savefig('NC_CC_spectrum.pdf')
