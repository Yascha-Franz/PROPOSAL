
#include <functional>

#include "PROPOSAL/crossection/WeakIntegral_NC.h"
#include "PROPOSAL/crossection/WeakInterpolant_NC.h"
#include "PROPOSAL/crossection/parametrization/WeakInteraction_NC.h"

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

WeakInterpolant_NC::WeakInterpolant_NC(const WeakInteraction_NC& param, InterpolationDef def)
        : CrossSectionInterpolant(DynamicData::WeakIntNC, param) {
    // Use own dNdx interpolation
    WeakInterpolant_NC::InitdNdxInterpolation(def);
}

WeakInterpolant_NC::WeakInterpolant_NC(const WeakInterpolant_NC& param)
        : CrossSectionInterpolant(param)
{
}

WeakInterpolant_NC::~WeakInterpolant_NC() {}

bool WeakInterpolant_NC::compare(const CrossSection& cross_section) const
{
    const WeakInterpolant_NC* cross_section_interpolant =
            static_cast<const WeakInterpolant_NC*>(&cross_section);

    if (dndx_interpolant_1d_.size() != cross_section_interpolant->dndx_interpolant_1d_.size())
        return false;
    else if (dndx_interpolant_2d_.size() != cross_section_interpolant->dndx_interpolant_2d_.size())
        return false;

    for (unsigned int i = 0; i < dndx_interpolant_1d_.size(); ++i)
    {
        if (*dndx_interpolant_1d_[i] != *cross_section_interpolant->dndx_interpolant_1d_[i])
            return false;
    }
    for (unsigned int i = 0; i < dndx_interpolant_2d_.size(); ++i)
    {
        if (*dndx_interpolant_2d_[i] != *cross_section_interpolant->dndx_interpolant_2d_[i])
            return false;
    }

    return true;
}

// ------------------------------------------------------------------------- //
void WeakInterpolant_NC::InitdNdxInterpolation(const InterpolationDef& def)
{
    // --------------------------------------------------------------------- //
    // Builder for dNdx
    // --------------------------------------------------------------------- //

    std::vector<Interpolant1DBuilder> builder1d(components_.size());
    std::vector<Interpolant2DBuilder> builder2d(components_.size());

    Helper::InterpolantBuilderContainer builder_container1d(components_.size());
    Helper::InterpolantBuilderContainer builder_container2d(components_.size());
    Helper::InterpolantBuilderContainer builder_return;

    Integral integral(IROMB, IMAXS, IPREC);

    for (unsigned int i = 0; i < components_.size(); ++i)
    {
        // !!! IMPORTANT !!!
        // Order of builder matter because the functions needed for 1d interpolation
        // needs the already intitialized 2d interpolants.
        builder2d[i]
            .SetMax1(def.nodes_cross_section)
            .SetX1Min(parametrization_->GetParticleDef().low)
            .SetX1Max(def.max_node_energy)
            .SetMax2(def.nodes_cross_section)
            .SetX2Min(0.0)
            .SetX2Max(1.0)
            .SetRomberg1(def.order_of_interpolation)
            .SetRational1(false)
            .SetRelative1(false)
            .SetIsLog1(true)
            .SetRomberg2(def.order_of_interpolation)
            .SetRational2(false)
            .SetRelative2(false)
            .SetIsLog2(false)
            .SetRombergY(def.order_of_interpolation)
            .SetRationalY(true)
            .SetRelativeY(false)
            .SetLogSubst(false)
            .SetFunction2D(std::bind(
                &CrossSectionInterpolant::FunctionToBuildDNdxInterpolant2D,
                this,
                std::placeholders::_1,
                std::placeholders::_2,
                std::ref(integral),
                i));

        builder_container2d[i].first  = &builder2d[i];
        builder_container2d[i].second = &dndx_interpolant_2d_[i];

        builder1d[i]
            .SetMax(def.nodes_cross_section)
            .SetXMin(parametrization_->GetParticleDef().low)
            .SetXMax(def.max_node_energy)
            .SetRomberg(def.order_of_interpolation)
            .SetRational(false)
            .SetRelative(false)
            .SetIsLog(true)
            .SetRombergY(def.order_of_interpolation)
            .SetRationalY(true)
            .SetRelativeY(false)
            .SetLogSubst(false)
            .SetFunction1D(std::bind(&CrossSectionInterpolant::FunctionToBuildDNdxInterpolant, this, std::placeholders::_1, i));

        builder_container1d[i].first  = &builder1d[i];
        builder_container1d[i].second = &dndx_interpolant_1d_[i];
    }

    builder_return.insert(builder_return.end(), builder_container2d.begin(), builder_container2d.end());
    builder_return.insert(builder_return.end(), builder_container1d.begin(), builder_container1d.end());
    // builder2d.insert(builder2d.end(), builder1d.begin(), builder1d.end());

    Helper::InitializeInterpolation("dNdx", builder_return, std::vector<Parametrization*>(1, parametrization_), def);
}

std::pair<std::vector<Particle*>, bool> WeakInterpolant_NC::CalculateProducedParticles(double energy, double energy_loss, const Vector3D initial_direction){
    // interaction energises a nukleus, which is not monitored in PROPOSAL

    return std::make_pair(std::vector<Particle*>{}, false);
}
