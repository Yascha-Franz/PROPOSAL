import pyPROPOSAL as pp
import numpy as np
import matplotlib.pyplot as plt

NuTauDef = pp.particle.NuTauDef.get()
TauDef = pp.particle.TauMinusDef.get()

Direction = pp.Vector3D(0,0,-1)
Position = pp.Vector3D(0,0,6.374134e12)
Energy = 1e10
statistic = 1

nutauprop = pp.Propagator(NuTauDef, "../../config_rock_neutrino.json")
tauprop = pp.Propagator(TauDef, "../../config_rock_lepton.json")

def propagate(prop, direction, pos, e):
    particle = prop.particle
    particle.direction = direction
    particle.position = pos
    particle.energy = e
    particle.propagated_distance = 0
    particle.time = 0

    particlepath = np.array([[pos.x/100, pos.y/100, pos.z/100]])

    secondarys = prop.propagate()

    for secondary in secondarys:
        particlepath = np.append(particlepath, [[secondary.position.x/100, secondary.position.y/100, secondary.position.z/100]], axis=0)

    return particlepath, secondarys

def propagate_and_plot(prop, direction, pos, e, colour):
    particlepath, secondarys = propagate(prop, direction, pos, e)
    plt.plot(particlepath[:,0],particlepath[:,2], c = colour)
    for secondary in secondarys:
        if secondary.id == pp.particle.Data.Particle:
            particleId = secondary.particle_def.particle_id
            if (abs(particleId) < 20 and abs(particleId) > 10 and particleId%2==1):
                propagate_and_plot(tauprop, secondary.direction, secondary.position, secondary.energy, 'k')
            elif (abs(particleId) < 20 and abs(particleId) > 10 and particleId%2==0):
                propagate_and_plot(nutauprop, secondary.direction, secondary.position, secondary.energy, 'b')
            
for i in range(0,statistic):
    propagate_and_plot(nutauprop, Direction, Position, Energy, 'b')

plt.xlabel('X')
plt.ylabel('Z')
plt.savefig('tauregeneration.pdf')
