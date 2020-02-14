
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

#include <functional>

#include "PROPOSAL/methods.h"

namespace PROPOSAL {

    class CrossSection;
    class WeakInteraction_NC;

    struct ParticleDef;
    class Medium;
    class EnergyCutSettings;

    class WeakInteractionFactory_NC
    {
    public:
        enum Enum
        {
            Fail = 0,
            None,
            CooperSarkarMertsch_NC,
        };

        struct Definition
        {
            Definition()
                    : parametrization(None)
                    , multiplier(1.0)
            {
            }

            bool operator==(const WeakInteractionFactory_NC::Definition& def) const
            {
                if (parametrization != def.parametrization)
                    return false;
                else if (multiplier != def.multiplier)
                    return false;

                return true;
            }

            bool operator!=(const WeakInteractionFactory_NC::Definition& def) const
            {
                return !(*this == def);
            }

            Enum parametrization;
            double multiplier;
        };

        // --------------------------------------------------------------------- //
        // Typedefs for readablitiy
        // --------------------------------------------------------------------- //

        typedef std::function<WeakInteraction_NC*(const ParticleDef&,
                                               const Medium&,
                                               double multiplier)>
                RegisterFunction;


        typedef std::map<std::string, RegisterFunction > Weak_NCMapString;
        typedef std::map<Enum, RegisterFunction > Weak_NCMapEnum;

        typedef Helper::Bimap<std::string, Enum> BimapStringEnum;

        // --------------------------------------------------------------------- //
        // Most general creation
        // --------------------------------------------------------------------- //

        CrossSection* CreateWeakInteraction_NC(const ParticleDef&,
                                            const Medium&,
                                            const Definition&) const;

        CrossSection* CreateWeakInteraction_NC(const ParticleDef&,
                                            const Medium&,
                                            const Definition&,
                                            InterpolationDef) const;


        // ----------------------------------------------------------------------------
        /// @brief string to enum conversation for photo parametrizations
        // ----------------------------------------------------------------------------
        Enum GetEnumFromString(const std::string&);

        // ----------------------------------------------------------------------------
        /// @brief enum to string conversation for photo parametrizations
        // ----------------------------------------------------------------------------
        std::string GetStringFromEnum(const Enum&);

        // --------------------------------------------------------------------- //
        // Singleton pattern
        // --------------------------------------------------------------------- //

        static WeakInteractionFactory_NC& Get()
        {
            static WeakInteractionFactory_NC instance;
            return instance;
        }

    private:
        WeakInteractionFactory_NC();
        ~WeakInteractionFactory_NC();

        // ----------------------------------------------------------------------------
        /// @brief Register Weak Interaction parametrizations
        ///
        /// @param name
        /// @param Enum
        /// @param RegisterFunction
        // ----------------------------------------------------------------------------
        void Register(const std::string& name, Enum, RegisterFunction);

        Weak_NCMapString weak_nc_map_str_;
        Weak_NCMapEnum weak_nc_map_enum_;
        BimapStringEnum string_enum_;
    };

} // namespace PROPOSAL
