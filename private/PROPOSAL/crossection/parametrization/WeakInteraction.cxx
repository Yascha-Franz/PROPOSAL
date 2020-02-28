
#include <cmath>

#include "PROPOSAL/crossection/parametrization/WeakInteraction.h"

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

WeakInteraction::WeakInteraction(const ParticleDef& particle_def,
                                 const Medium& medium,
                                 double multiplier)
        : Parametrization(particle_def, medium, EnergyCutSettings(), multiplier, true)
{
}

WeakInteraction::WeakInteraction(const WeakInteraction& param)
        : Parametrization(param)
{
}

WeakInteraction::~WeakInteraction() {}

bool WeakInteraction::compare(const Parametrization& parametrization) const
{
    return Parametrization::compare(parametrization);
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

Parametrization::IntegralLimits WeakInteraction::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    double aux = (MP + MN) / 2; //for isoscalar targets
    aux = 2 * energy * aux + pow(aux, 2);

    limits.vMin = 1e6 / (aux ); //q^2_min = 1e6 MeV for experimental reasons
    limits.vUp = limits.vMin; //treat interaction as fully stochastic
    limits.vMax = 1;

    return limits;
}

size_t WeakInteraction::GetHash() const
{
    size_t seed = Parametrization::GetHash();
    hash_combine(seed, particle_def_.charge);

    return seed;
}

// ------------------------------------------------------------------------- //
// Specific implementations
// ------------------------------------------------------------------------- //
double WeakInteraction::GetWeakPartnerCharge(const ParticleDef& particle_def){
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

WeakCooperSarkarMertsch::WeakCooperSarkarMertsch(const ParticleDef& particle_def,
                                                 const Medium& medium,
                                                 double multiplier)
        : WeakInteraction(particle_def, medium, multiplier)
        , interpolant_(2, NULL)
{
    if(particle_def.charge < 0.|| GetWeakPartnerCharge(particle_def) > 0.)
    {
        // Initialize interpolant for particles (remember crossing symmetry rules)
        interpolant_[0] = new Interpolant(energies, y_nubar_p, sigma_nubar_p, IROMB, false, false, IROMB, false, false);
        interpolant_[1] = new Interpolant(energies, y_nubar_n, sigma_nubar_n, IROMB, false, false, IROMB, false, false);
    }
    else if(particle_def.charge > 0.|| GetWeakPartnerCharge(particle_def) < 0.){
        // Initialize interpolant for antiparticles (remember crossing symmetry rules)
        interpolant_[0] = new Interpolant(energies, y_nu_p, sigma_nu_p, IROMB, false, false, IROMB, false, false);
        interpolant_[1] = new Interpolant(energies, y_nu_n, sigma_nu_n, IROMB, false, false, IROMB, false, false);
    }else{
        log_fatal("Weak interaction: Particle to propagate is not a charged lepton");
    }

}

WeakCooperSarkarMertsch::WeakCooperSarkarMertsch(const WeakCooperSarkarMertsch& param)
        : WeakInteraction(param)
        , interpolant_()
{
    interpolant_.resize(param.interpolant_.size());

    for (unsigned int i = 0; i < param.interpolant_.size(); ++i)
    {
        interpolant_[i] = new Interpolant(*param.interpolant_[i]);
    }
}

WeakCooperSarkarMertsch::~WeakCooperSarkarMertsch()
{
    for (std::vector<Interpolant*>::const_iterator iter = interpolant_.begin(); iter != interpolant_.end(); ++iter)
    {
        delete *iter;
    }

    interpolant_.clear();
}

bool WeakCooperSarkarMertsch::compare(const Parametrization& parametrization) const
{
    const WeakCooperSarkarMertsch* weak = static_cast<const WeakCooperSarkarMertsch*>(&parametrization);

    if (interpolant_.size() != weak->interpolant_.size())
        return false;

    for (unsigned int i = 0; i < interpolant_.size(); ++i)
    {
        if (*interpolant_[i] != *weak->interpolant_[i])
            return false;
    }

    return WeakInteraction::compare(parametrization);
}

double WeakCooperSarkarMertsch::DifferentialCrossSection(double energy, double v)
{
    double proton_contribution = components_[component_index_]->GetNucCharge() * interpolant_.at(0)->InterpolateArray(std::log10(energy), v);
    double neutron_contribution = (components_[component_index_]->GetAtomicNum() - components_[component_index_]->GetNucCharge()) * interpolant_.at(1)->InterpolateArray(std::log10(energy), v);
    double total_contribution = (proton_contribution + neutron_contribution);

    return medium_->GetMassDensity() * medium_->GetMolDensity() * components_[component_index_]->GetAtomInMolecule() * 1e-36 * std::max(0.0, total_contribution); //factor 1e-36: conversion from pb to cm^2
}


const std::string WeakCooperSarkarMertsch::name_ = "WeakCooperSarkarMertsch";

