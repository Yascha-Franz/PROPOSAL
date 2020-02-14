
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

#include "PROPOSAL/crossection/CrossSectionInterpolant.h"

namespace PROPOSAL {

    class WeakInteraction_NC;

    class WeakInterpolant_NC : public CrossSectionInterpolant
    {
    public:
        WeakInterpolant_NC(const WeakInteraction_NC&, InterpolationDef);
        WeakInterpolant_NC(const WeakInterpolant_NC&);
        virtual ~WeakInterpolant_NC();

        CrossSection* clone() const { return new WeakInterpolant_NC(*this); }

        // ----------------------------------------------------------------- //
        // Public methods
        // ----------------------------------------------------------------- //

        //these methods return zero because the weak interaction contribution is stochastic only
        double CalculatedEdx(double energy){ (void)energy; return 0; }
        double CalculatedEdxWithoutMultiplier(double energy){ (void)energy; return 0; }
        double CalculatedE2dx(double energy){ (void)energy; return 0; }
        std::pair<std::vector<Particle*>, bool> CalculateProducedParticles(double energy, double energy_loss, const Vector3D initial_direction);
    protected:
        virtual bool compare(const CrossSection&) const;
        virtual void InitdNdxInterpolation(const InterpolationDef& def);
    };

} // namespace PROPOSAL
