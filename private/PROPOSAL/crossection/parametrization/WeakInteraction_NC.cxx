
#include <cmath>

#include "PROPOSAL/crossection/parametrization/WeakInteraction_NC.h"

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/crossection/parametrization/ParamTables.h"
#include <iostream>


using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

WeakInteraction_NC::WeakInteraction_NC(const ParticleDef& particle_def,
                                 const Medium& medium,
                                 double multiplier)
        : Parametrization(particle_def, medium, EnergyCutSettings(), multiplier, false)
{
}

WeakInteraction_NC::WeakInteraction_NC(const WeakInteraction_NC& param)
        : Parametrization(param)
{
}

WeakInteraction_NC::~WeakInteraction_NC() {}

bool WeakInteraction_NC::compare(const Parametrization& parametrization) const
{
    return Parametrization::compare(parametrization);
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

Parametrization::IntegralLimits WeakInteraction_NC::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    double aux = (MP + MN) / 2; //for isoscalar targets
    aux = 2 * energy * aux + pow(aux, 2);

    limits.vMin = 1e6 / (aux ); //q^2_min = 1e6 MeV for experimental reasons
    limits.vUp = limits.vMin; //treat interaction as fully stochastic
    limits.vMax = 1;

    return limits;
}

size_t WeakInteraction_NC::GetHash() const
{
    size_t seed = Parametrization::GetHash();
    hash_combine(seed, particle_def_.charge);

    return seed;
}

// ------------------------------------------------------------------------- //
// Specific implementations
// ------------------------------------------------------------------------- //
double WeakInteraction_NC::GetWeakPartnerCharge(const ParticleDef& particle_def){
    if(particle_def.weak_partner == static_cast<int>(ParticleType::None)){
        log_fatal("WeakPartner not defined for particle %s.", particle_def.name.c_str());
        return 0; //To avoid warnings
    }
    else{
        std::map<const int,const ParticleDef&>::iterator it = Id_Particle_Map.find(particle_def.weak_partner);
        if(it->first == particle_def.weak_partner)
            return (it->second).charge;
        log_fatal("WeakPartner (Id: %s) is not defined in PROPOSAL", particle_def.weak_partner);
        return 0; //To avoid warnings
    }
}

WeakCooperSarkarMertsch_NC::WeakCooperSarkarMertsch_NC(const ParticleDef& particle_def,
                                                 const Medium& medium,
                                                 double multiplier)
        : WeakInteraction_NC(particle_def, medium, multiplier)
        , interpolant_(2, NULL)
{
    if(particle_def.charge < 0.|| GetWeakPartnerCharge(particle_def) > 0.)
    {
        // Initialize interpolant for particles (remember crossing symmetry rules)
        interpolant_[0] = new Interpolant(energies, y_nubar_p_NC, sigma_nubar_p_NC, IROMB, false, false, IROMB, false, false);
        interpolant_[1] = new Interpolant(energies, y_nubar_n_NC, sigma_nubar_n_NC, IROMB, false, false, IROMB, false, false);
    }
    else if(particle_def.charge > 0.|| GetWeakPartnerCharge(particle_def) < 0.){
        // Initialize interpolant for antiparticles (remember crossing symmetry rules)
        interpolant_[0] = new Interpolant(energies, y_nu_p_NC, sigma_nu_p_NC, IROMB, false, false, IROMB, false, false);
        interpolant_[1] = new Interpolant(energies, y_nu_n_NC, sigma_nu_n_NC, IROMB, false, false, IROMB, false, false);
    }else{
        log_fatal("Weak interaction: Particle to propagate is not a charged lepton");
    }

}

WeakCooperSarkarMertsch_NC::WeakCooperSarkarMertsch_NC(const WeakCooperSarkarMertsch_NC& param)
        : WeakInteraction_NC(param)
        , interpolant_()
{
    interpolant_.resize(param.interpolant_.size());

    for (unsigned int i = 0; i < param.interpolant_.size(); ++i)
    {
        interpolant_[i] = new Interpolant(*param.interpolant_[i]);
    }
}

WeakCooperSarkarMertsch_NC::~WeakCooperSarkarMertsch_NC()
{
    for (std::vector<Interpolant*>::const_iterator iter = interpolant_.begin(); iter != interpolant_.end(); ++iter)
    {
        delete *iter;
    }

    interpolant_.clear();
}

bool WeakCooperSarkarMertsch_NC::compare(const Parametrization& parametrization) const
{
    const WeakCooperSarkarMertsch_NC* weak = static_cast<const WeakCooperSarkarMertsch_NC*>(&parametrization);

    if (interpolant_.size() != weak->interpolant_.size())
        return false;

    for (unsigned int i = 0; i < interpolant_.size(); ++i)
    {
        if (*interpolant_[i] != *weak->interpolant_[i])
            return false;
    }

    return WeakInteraction_NC::compare(parametrization);
}

double WeakCooperSarkarMertsch_NC::DifferentialCrossSection(double energy, double v)
{
    double proton_contribution = components_[component_index_]->GetNucCharge() * interpolant_.at(0)->InterpolateArray(std::log10(energy), v);
    double neutron_contribution = (components_[component_index_]->GetAtomicNum() - components_[component_index_]->GetNucCharge()) * interpolant_.at(1)->InterpolateArray(std::log10(energy), v);
    double mean_contribution = (proton_contribution + neutron_contribution) / (components_[component_index_]->GetAtomicNum());
    
    return medium_->GetMolDensity() * components_[component_index_]->GetAtomInMolecule() * 1e-36 * std::max(0.0, mean_contribution); //factor 1e-36: conversion from pb to cm^2
}


const std::string WeakCooperSarkarMertsch_NC::name_ = "WeakCooperSarkarMertsch_NC";

