
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/Secondaries.h"

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/decay/DecayChannel.h"

#include "PROPOSAL/medium/density_distr/density_distr.h"

#include "PROPOSAL/particle/Particle.h"

#include "PROPOSAL/geometry/Geometry.h"
#include "PROPOSAL/geometry/Sphere.h"

#include "PROPOSAL/SectorStochastic.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/methods.h"
#include "PROPOSAL/propagation_utility/ContinuousRandomizer.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

#include <algorithm>
#include <array>
#include <memory>
#include <utility>
using namespace PROPOSAL;

/******************************************************************************
 *                                 SectorStochastic                                 *
 ******************************************************************************/
// ------------------------------------------------------------------------- //
// Constructors
// ------------------------------------------------------------------------- //

SectorStochastic::SectorStochastic(const ParticleDef& particle_def, const Definition& sector_def)
    : Sector::Sector(particle_def, sector_def)
{
}

SectorStochastic::SectorStochastic(const ParticleDef& particle_def,
               const Definition& sector_def,
               const InterpolationDef& interpolation_def)
    : Sector::Sector(particle_def, sector_def, interpolation_def)
{
}

SectorStochastic::SectorStochastic(const ParticleDef& particle_def, const Sector& sector)
    : Sector::Sector(particle_def, sector)
{
}

SectorStochastic::SectorStochastic(const Sector& sector)
    : Sector::Sector(sector)
{
}

SectorStochastic::~SectorStochastic()
{
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %                          SectorStochastic Utilities                     %
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double SectorStochastic::CalculateTime(const DynamicData& p_condition,
    const double final_energy, const double displacement)
{
    if (exact_time_calculator_) {
        
        double square_momentum = (final_energy - utility_.GetParticleDef().mass) *
                                 (final_energy + utility_.GetParticleDef().mass);
        double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));
        double propagating_time = displacement * final_energy / (particle_momentum * SPEED);

        return p_condition.GetTime()
            + propagating_time
            / utility_.GetMedium().GetDensityDistribution().Evaluate(
                  p_condition.GetPosition());
    }

    return p_condition.GetTime() + displacement / SPEED;
}
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %                          LossEnergies                                   %
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double SectorStochastic::CalculateDecay(const double energy, const double rnd)
{
    double rndd = -std::log(rnd);

    // solving the tracking integral
    if (particle_def_.lifetime < 0) {
        return INFINITY;
    }

    double square_momentum =
        (energy - particle_def_.mass) * (energy + particle_def_.mass);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));

    double aux = 1.0 / std::max((particle_momentum / particle_def_.mass) *
                             particle_def_.lifetime * SPEED,
                         PARTICLE_POSITION_RESOLUTION);

    return aux/rndd;
}

double SectorStochastic::CalculateInteraction(const double energy, const double rnd)
{
    double rndi = -std::log(rnd);

    double total_rate = 0.0;
    const std::vector<CrossSection*>& crosssections = utility_.GetCrosssections();
    for (std::vector<CrossSection*>::const_iterator iter = crosssections.begin(); iter != crosssections.end(); ++iter) {
        total_rate += (*iter)->CalculatedNdx(energy);
    }

    return rndi / total_rate;
}

double SectorStochastic::CalculateMinimalEnergy(const double current_energy, const double cut)
{
    if(current_energy > cut){
        return INFINITY;
    }
    return 0;
}

double SectorStochastic::CalculateDistance(const double energy, const double distance) {
    (void) energy;
    return distance;
}

int SectorStochastic::ChooseLoss(const std::array<double, 4>& LossLengths)
{
    const auto minmax
        = std::minmax_element(begin(LossLengths), end(LossLengths));

    if (*minmax.first == LossLengths[LossType::Decay]) {
        return LossType::Decay;
    } else {
        return std::distance(LossLengths.begin(), minmax.first);
    }
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %                               Do Loss                                   %
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std::shared_ptr<DynamicData> SectorStochastic::DoDisplacement(
    const DynamicData& p_condition, double energy, double displacement)
{
    double dist = p_condition.GetPropagatedDistance() + displacement;
    double time = CalculateTime(p_condition, energy, displacement);

    Vector3D position{ p_condition.GetPosition() };
    Vector3D direction{ p_condition.GetDirection() };

    /*p_condition.SetPosition(p_condition.GetPosition() + displacement * p_condition.GetDirection());
    p_condition.SetTime(time);
    p_condition.SetPropagatedDistance(dist);*/
    return std::make_shared<DynamicData>(
        static_cast<int>(InteractionType::ContinuousEnergyLoss), position + direction * displacement,
        direction, energy, energy, time, dist);
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %                               Propagate                                 %
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Secondaries SectorStochastic::Propagate(
    const DynamicData& p_initial, double distance, const double minimal_energy)
{
    Secondaries secondaries(std::make_shared<ParticleDef>(particle_def_));

    auto p_condition = std::make_shared<DynamicData>(p_initial);
    double dist_limit{ p_initial.GetPropagatedDistance() + distance };
    double rnd;
    int Loss;
    std::array<double, 4> LossLengths;

    while (true) {
        rnd = RandomGenerator::Get().RandomDouble();
        LossLengths[LossType::Decay]
            = CalculateDecay(p_condition->GetEnergy(), rnd);

        rnd = RandomGenerator::Get().RandomDouble();
        LossLengths[LossType::Interaction]
            = CalculateInteraction(p_condition->GetEnergy(), rnd);

        distance = CalculateDistance(p_condition->GetEnergy(), dist_limit - p_condition->GetPropagatedDistance());
        LossLengths[LossType::Distance] = distance;

        LossLengths[LossType::MinimalE]
            = CalculateMinimalEnergy(p_condition->GetEnergy(), minimal_energy);

        Loss = ChooseLoss(LossLengths);

        p_condition
            = DoDisplacement(*p_condition, p_condition->GetEnergy(), LossLengths[Loss]);

        if (Loss == LossType::Interaction)
        {
            p_condition = Sector::DoInteraction(*p_condition);
            if(utility_.GetCrosssection(p_condition->GetTypeId())->GetParametrization().isFatal()) break;
            secondaries.push_back(*p_condition);
        }
        else
        {
            break;
        }
    };

    if (Loss == LossType::Decay)
    {
        p_condition = Sector::DoDecay(*p_condition);
    }

    secondaries.push_back(*p_condition);

    return secondaries;
}

/*bool SectorStochastic::only_stochastic_loss(){
    std::vector<CrossSection*> cross_sections = utility_.GetCrosssections();
    for (unsigned int i = 0; i < cross_sections.size(); i++) {
        if(cross_sections[i]->CalculatedEdx(particle_.GetEnergy())!=0){
            return false;
        }
    }
    return true;
}*/