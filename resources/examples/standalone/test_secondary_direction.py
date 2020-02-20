import pyPROPOSAL as pp
import numpy as np
import matplotlib.pyplot as plt

NuTauDef = pp.particle.NuTauDef.get()
TauDef = pp.particle.TauMinusDef.get()

Direction = pp.Vector3D(0,0,-1)
Position = pp.Vector3D(0,0,6.374134e12)
statistic = 10**(4)
Energy_bins = 10**(np.linspace(10, 11 ,2))

#nutauprop = pp.Propagator(NuTauDef, "../../config_rock_neutrino.json")
tauprop = pp.Propagator(TauDef, "../../config_rock_lepton.json")
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

direction = pp.Vector3D(0,0,-1)
count_wrong = 0
for i in range(0,np.size(Energy_bins)):
    for j in range(0, statistic):
        #print(i, '  ', np.log10(j+1))
        secondarys = propagate(tauprop, Direction, Position, Energy_bins[i])
        break_ = False
        direction = secondarys.particles[-1].direction
        #print(np.log10(j+1))
        if(direction.z>-0.9):
            count_wrong+=1
            for particle in secondarys.particles:
                print(particle.name)
            print(direction)
            print(count_wrong, "/", j)
            print()

print("Hallo")