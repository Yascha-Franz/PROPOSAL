
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

#include "PROPOSAL/SectorMixed.h"
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
 *                                 SectorMixed                                *
 ******************************************************************************/
// ------------------------------------------------------------------------- //
// Constructors
// ------------------------------------------------------------------------- //

SectorMixed::SectorMixed(const ParticleDef& particle_def, const Definition& sector_def)
    : Sector::Sector(particle_def, sector_def)
{
    displacement_calculator_ = new UtilityIntegralDisplacement(utility_);
    interaction_calculator_  = new UtilityIntegralInteraction(utility_);
}

SectorMixed::SectorMixed(const ParticleDef& particle_def,
               const Definition& sector_def,
               const InterpolationDef& interpolation_def)
    : Sector::Sector(particle_def, sector_def, interpolation_def)
{
    displacement_calculator_ = new UtilityInterpolantDisplacement(utility_, interpolation_def);
    interaction_calculator_  = new UtilityInterpolantInteraction(utility_, interpolation_def);
}

SectorMixed::SectorMixed(const ParticleDef& particle_def, const Sector& sector)
    : Sector::Sector(particle_def, sector)
{
}

SectorMixed::SectorMixed(const Sector& sector)
    : Sector::Sector(sector)
{
}

SectorMixed::~SectorMixed()
{
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %                          SectorMixed Utilities                          %
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double SectorMixed::CalculateTime(const DynamicData& p_condition,
    const double final_energy, const double displacement)
{
    if (exact_time_calculator_) {
        // DensityDistribution Approximation: Use the DensityDistribution at the
        // position of initial energy
        return p_condition.GetTime()
            + exact_time_calculator_->Calculate(
                  p_condition.GetEnergy(), final_energy, 0.0)
            / utility_.GetMedium().GetDensityDistribution().Evaluate(
                  p_condition.GetPosition());
    }

    return p_condition.GetTime() + displacement / SPEED;
}

void SectorMixed::Scatter(const double displacement, const double initial_energy,
    const double final_energy, Vector3D& position, Vector3D& direction)
{
    if (sector_def_.scattering_model != ScatteringFactory::Enum::NoScattering) {
        Directions directions = scattering_->Scatter(
            displacement, initial_energy, final_energy, position, direction);
        position = position + displacement * directions.u_;
        direction = directions.n_i_;
    } else {
        position = position + displacement * direction;
    }
}

double SectorMixed::ContinuousRandomize(
    const double initial_energy, const double final_energy)
{
    if (cont_rand_) {
        if (final_energy != particle_def_.low) {
            double rnd = RandomGenerator::Get().RandomDouble();
            return cont_rand_->Randomize(initial_energy, final_energy, rnd);
        }
    }
    return final_energy;
}

double SectorMixed::Displacement(const DynamicData& p_condition,
    const double final_energy, const double border_length)
{
    return displacement_calculator_->Calculate(p_condition.GetEnergy(),
        final_energy, border_length, p_condition.GetPosition(),
        p_condition.GetDirection());
}
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %                          LossEnergies                                   %
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double SectorMixed::CalculateDecay(const double initial_energy, const double rnd)
{
    double rndd = -std::log(rnd);
    double rnddMin = 0;

    // solving the tracking integral
    if (particle_def_.lifetime < 0) {
        return particle_def_.low;
    }

    rnddMin
        = decay_calculator_->Calculate(initial_energy, particle_def_.low, rndd);

    // evaluating the energy loss
    if (rndd >= rnddMin || rnddMin <= 0) {
        return particle_def_.low;
    }

    return decay_calculator_->GetUpperLimit(initial_energy, rndd);
}

double SectorMixed::CalculateInteraction(const double initial_energy, const double rnd)
{
    double rndi = -std::log(rnd);
    double rndiMin = 0;

    // solving the tracking integral
    rndiMin = interaction_calculator_->Calculate(
        initial_energy, particle_def_.low, rndi);

    if (rndi >= rndiMin || rndiMin <= 0) {
        return particle_def_.low;
    }

    return interaction_calculator_->GetUpperLimit(initial_energy, rndi);
}

double SectorMixed::CalculateMinimalEnergy(const double current_energy, const double cut)
{
    return std::min({ current_energy, cut });
}

double SectorMixed::CalculateDistance(
    const double initial_energy, const double distance)
{
    return displacement_calculator_->GetUpperLimit(initial_energy, distance);
}

int SectorMixed::ChooseLoss(const std::array<double, 4>& LossEnergies)
{
    const auto minmax
        = std::minmax_element(begin(LossEnergies), end(LossEnergies));

    if (*minmax.second == LossEnergies[LossType::Decay]) {
        return LossType::Decay;
    } else {
        return std::distance(LossEnergies.begin(), minmax.second);
    }
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %                               Do Loss                                   %
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std::shared_ptr<DynamicData> SectorMixed::DoContinuous(
    const DynamicData& p_condition, double final_energy, double sector_border)
{

    double initial_energy{ p_condition.GetEnergy() };
    double displacement
        = Displacement(p_condition, final_energy, sector_border);

    double dist = p_condition.GetPropagatedDistance() + displacement;
    double time = CalculateTime(p_condition, final_energy, displacement);

    Vector3D position{ p_condition.GetPosition() };
    Vector3D direction{ p_condition.GetDirection() };

    Scatter(displacement, p_condition.GetEnergy(), final_energy, position,
        direction);
    final_energy = ContinuousRandomize(p_condition.GetEnergy(), final_energy);

    return std::make_shared<DynamicData>(
        static_cast<int>(InteractionType::ContinuousEnergyLoss), position,
        direction, final_energy, initial_energy, time, dist);
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %                               Propagate                                 %
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Secondaries SectorMixed::Propagate(
    const DynamicData& p_initial, double distance, const double minimal_energy)
{
    Secondaries secondaries(std::make_shared<ParticleDef>(particle_def_));

    auto p_condition = std::make_shared<DynamicData>(p_initial);
    double dist_limit{ p_initial.GetPropagatedDistance() + distance };
    double rnd;
    int minimalLoss;
    std::array<double, 4> LossEnergies;

    while (true) {
        rnd = RandomGenerator::Get().RandomDouble();
        LossEnergies[LossType::Decay]
            = CalculateDecay(p_condition->GetEnergy(), rnd);

        rnd = RandomGenerator::Get().RandomDouble();
        LossEnergies[LossType::Interaction]
            = CalculateInteraction(p_condition->GetEnergy(), rnd);

        distance = dist_limit - p_condition->GetPropagatedDistance();
        LossEnergies[LossType::Distance]
            = CalculateDistance(p_condition->GetEnergy(), distance);

        LossEnergies[LossType::MinimalE]
            = CalculateMinimalEnergy(p_condition->GetEnergy(), minimal_energy);

        minimalLoss = ChooseLoss(LossEnergies);

        p_condition
            = DoContinuous(*p_condition, LossEnergies[minimalLoss], distance);
        if (sector_def_.do_continuous_energy_loss_output)
            secondaries.push_back(*p_condition);

        if (minimalLoss == LossType::Interaction)
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

    if (minimalLoss == LossType::Decay)
    {
        p_condition = Sector::DoDecay(*p_condition);
    }

    secondaries.push_back(*p_condition);

    return secondaries;
}