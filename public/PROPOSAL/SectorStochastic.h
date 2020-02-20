
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
#include "PROPOSAL/Sector.h"

namespace PROPOSAL {
/*! \class ProcessSector ProcessSector.h "CrossSections.h"
    \brief initializes all cross sections and keeps references to them
 */
class SectorStochastic : public Sector{
public:
    // Sector(Particle&);
    SectorStochastic(const ParticleDef&, const Sector::Definition&);
    SectorStochastic(const ParticleDef&, const Sector::Definition&, const InterpolationDef&);
    SectorStochastic(const ParticleDef&, const Sector&);
    // Sector(Particle&, const Geometry&, const Utility&, const Scattering&,
    // bool do_interpolation, const Sector::Definition& def = Sector::Definition());
    SectorStochastic(const Sector&);
    ~SectorStochastic();
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
    double CalculateTime(const DynamicData& p_condition,
        const double final_energy, const double displacement);

    // Loss Lengths
    double CalculateMinimalEnergy(const double inital_energy, const double cut);
    double CalculateDecay(const double initial_energy, const double rnd);
    double CalculateInteraction(const double initial_energy, const double rnd);
    double CalculateDistance(const double initial_energy, const double distance);
    int ChooseLoss(const std::array<double, 4>& LossLengths);

    std::shared_ptr<DynamicData> DoDisplacement(
    const DynamicData& p_condition, double energy, double displacement);


    std::shared_ptr<DynamicData> DoContinuous(
        const DynamicData&, double, double);
    /* std::shared_ptr<DynamicData> DoBorder(const DynamicData& ); */

    Secondaries Propagate(const DynamicData& particle_condition,
        double max_distance=1e20, double minimal_energy=0.);

protected:
    Sector& operator=(const Sector&); // Undefined & not allowed
};
} // namespace PROPOSAL
