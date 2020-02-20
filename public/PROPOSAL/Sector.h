
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/

#pragma once

// #include <string>
// #include <vector>
#include <memory>
#include <tuple>

#include "PROPOSAL/Secondaries.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/scattering/ScatteringFactory.h"

#include "PROPOSAL/propagation_utility/PropagationUtility.h"

namespace PROPOSAL {

class ContinuousRandomizer;
// class CrossSection;
// class Medium;
// class EnergyCutSettings;
class Geometry;

enum LossType : int {
    Interaction = 0,
    Decay = 1,
    Distance = 2,
    MinimalE = 3,
};


/*! \class ProcessSector ProcessSector.h "CrossSections.h"
    \brief initializes all cross sections and keeps references to them
 */
class Sector {
public:
    struct ParticleLocation {
        enum Enum { InfrontDetector = 0, InsideDetector, BehindDetector };
    };

    class Definition {
    public:
        Definition();
        Definition(const Definition&);
        ~Definition();

        bool operator==(const Definition&) const;
        bool operator!=(const Definition&) const;
        Definition& operator=(const Definition&);
        void swap(Definition&);

        void SetMedium(const Medium&);
        const Medium& GetMedium() const { return *medium_; }

        void SetGeometry(const Geometry&);
        const Geometry& GetGeometry() const { return *geometry_; }

        bool do_stochastic_loss_weighting; //!< Do weigthing of stochastic
                                           //!< losses. Set to false in
                                           //!< constructor.
        double stochastic_loss_weighting;  //!< weigth of stochastic losses. Set
                                           //!< to 0 in constructor
        bool stopping_decay; //!< Let particle decay if elow is reached but no
                             //!< decay was sampled

        bool do_continuous_randomization;
        bool do_continuous_energy_loss_output;
        bool do_exact_time_calculation;
        bool only_loss_inside_detector;

        ScatteringFactory::Enum scattering_model;

        Sector::ParticleLocation::Enum location;

        Utility::Definition utility_def;

        EnergyCutSettings cut_settings;

    private:
        Medium* medium_;
        Geometry* geometry_;
    };

public:
    // Sector(Particle&);
    Sector(const ParticleDef&, const Definition&);
    Sector(const ParticleDef&, const Definition&, const InterpolationDef&);
    Sector(const ParticleDef&, const Sector&);
    // Sector(Particle&, const Geometry&, const Utility&, const Scattering&,
    // bool do_interpolation, const Definition& def = Definition());
    Sector(const Sector&);
    virtual ~Sector();

    bool operator==(const Sector&) const;
    bool operator!=(const Sector&) const;

    // Sector& operator=(const Sector& collection);
    // friend std::ostream& operator<<(std::ostream& os, Sector const&
    // collection);

    // --------------------------------------------------------------------- //
    // Member functions
    // --------------------------------------------------------------------- //

    /**
     * Propagates the particle of initial energy e to the distance r.
     * Returns the final energy if the
     * particle has survived or the track length to the
     * point of disappearance with a minus sign otherwise.
     *
     *  \param  distance   maximum track length
     *  \param  energy   initial energy
     *  \return energy at distance OR -(track length)
     */

    // Utilites
    virtual double CalculateTime(const DynamicData& p_condition,
        const double final_energy, const double displacement) = 0;
    double BorderLength(const Vector3D&, const Vector3D& );

    virtual double CalculateMinimalEnergy(const double inital_energy, const double cut) = 0;
    virtual double CalculateDecay(const double initial_energy, const double rnd) = 0;
    virtual double CalculateInteraction(const double initial_energy, const double rnd) = 0;
    virtual double CalculateDistance(const double initial_energy, const double distance) = 0;
    virtual int ChooseLoss(const std::array<double, 4>&) = 0;


    std::shared_ptr<DynamicData> DoInteraction(const DynamicData&);
    std::shared_ptr<DynamicData> DoDecay(const DynamicData&);
    /* std::shared_ptr<DynamicData> DoBorder(const DynamicData& ); */

    virtual Secondaries Propagate(const DynamicData& particle_condition,
        double max_distance=1e20, double minimal_energy=0.) = 0;

    /**
     *  Makes Stochastic Energyloss
     *
     *  \return tuple of energy loss [MeV], kind of interaction and list of
     * produced particles
     */
    std::pair<double, int> MakeStochasticLoss(double particle_energy);

    // --------------------------------------------------------------------- //
    // Enable options & Setter
    // --------------------------------------------------------------------- //

    void SetLocation(ParticleLocation::Enum location)
    {
        sector_def_.location = location;
    }

    // --------------------------------------------------------------------- //
    // Getter
    // --------------------------------------------------------------------- //

    ParticleLocation::Enum GetLocation() const { return sector_def_.location; }

    Scattering* GetScattering() const { return scattering_; }
    const ParticleDef GetParticleDef() const { return particle_def_; }
    Geometry* GetGeometry() const { return geometry_; }
    const Utility& GetUtility() const { return utility_; }
    const Medium* GetMedium() const { return &utility_.GetMedium(); }
    const Definition& GetSectorDef() const { return sector_def_; }
    Definition& GetSectorDef() { return sector_def_; }

protected:
    Sector& operator=(const Sector&); // Undefined & not allowed

    // --------------------------------------------------------------------- //
    // Protected members
    // --------------------------------------------------------------------- //

    Definition sector_def_;

    ParticleDef particle_def_;
    Geometry* geometry_;

    Utility utility_;
    UtilityDecorator* displacement_calculator_;
    UtilityDecorator* interaction_calculator_;
    UtilityDecorator* decay_calculator_;
    UtilityDecorator* exact_time_calculator_;

    ContinuousRandomizer* cont_rand_;
    Scattering* scattering_;

    std::pair<double, double> produced_particle_moments_{ 100., 10000. };
    unsigned int n_th_call_{ 1 };
};
} // namespace PROPOSAL
