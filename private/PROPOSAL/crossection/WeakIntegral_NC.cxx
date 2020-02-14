
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/WeakIntegral_NC.h"
#include "PROPOSAL/crossection/parametrization/WeakInteraction_NC.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

WeakIntegral_NC::WeakIntegral_NC(const WeakInteraction_NC& param)
        : CrossSectionIntegral(DynamicData::WeakIntNC, param)
{
}

WeakIntegral_NC::WeakIntegral_NC(const WeakIntegral_NC& weak)
        : CrossSectionIntegral(weak)
{
}

WeakIntegral_NC::~WeakIntegral_NC() {}

std::pair<std::vector<Particle*>, bool> WeakIntegral_NC::CalculateProducedParticles(double energy, double energy_loss, const Vector3D initial_direction){
    // interaction energises a nukleus, which is not monitored in PROPOSAL

    return std::make_pair(std::vector<Particle*>{}, false);
}