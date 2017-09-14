
#include <boost/bind.hpp>

#include <cmath>

#include "PROPOSAL/crossection/IonizInterpolant.h"
#include "PROPOSAL/crossection/IonizIntegral.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"

#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Output.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

IonizInterpolant::IonizInterpolant(const Parametrization& param)
    : CrossSectionInterpolant(DynamicData::DeltaE, param)
{
    Parametrization::Definition param_def = parametrization_->GetDefinition();
    Interpolant1DBuilder builder1d;
    Helper::InterpolantBuilderContainer builder_container;

    IonizIntegral ioniz(param);

    builder1d.SetMax(NUM1)
        .SetXMin(param.GetParticleDef().low)
        .SetXMax(BIGENERGY)
        .SetRomberg(param_def.order_of_interpolation)
        .SetRational(true)
        .SetRelative(false)
        .SetIsLog(true)
        .SetRombergY(param_def.order_of_interpolation)
        .SetRationalY(false)
        .SetRelativeY(false)
        .SetLogSubst(true)
        .SetFunction1D(boost::bind(&CrossSection::CalculatedEdx, &ioniz, _1));

    builder_container.push_back(std::make_pair(&builder1d, &dedx_interpolant_));

    log_info("Initialize dEdx for %s", typeid(parametrization_).name());
    Helper::InitializeInterpolation("dEdx", builder_container, std::vector<Parametrization*>(1, parametrization_));
    log_info("Initialization dEdx for %s done!", typeid(parametrization_).name());
}

IonizInterpolant::IonizInterpolant(const IonizInterpolant& brems): CrossSectionInterpolant(brems)
{
}

IonizInterpolant::~IonizInterpolant()
{
}

// ----------------------------------------------------------------- //
// Public methods
// ----------------------------------------------------------------- //

double IonizInterpolant::CalculatedEdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    return std::max(dedx_interpolant_->Interpolate(energy), 0.);
}

// ------------------------------------------------------------------------- //
double IonizInterpolant::CalculatedNdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    return std::max(dndx_interpolant_1d_[0]->Interpolate(energy), 0.);
}

// ------------------------------------------------------------------------- //
double IonizInterpolant::CalculatedNdx(double energy, double rnd)
{
    (void) rnd;

    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    return std::max(dndx_interpolant_1d_[0]->Interpolate(energy), 0.);
}

// ------------------------------------------------------------------------- //
double IonizInterpolant::FunctionToBuildDNdxInterpolant(double energy, int component)
{
    (void) component;
    return dndx_interpolant_2d_[0]->Interpolate(energy, 1.0);
}

// ------------------------------------------------------------------------- //
double IonizInterpolant::FunctionToBuildDNdxInterpolant2D(double energy, double v, Integral& integral, int component)
{
    (void)component;
    std::cout << "Ioniz FunctionToBuildDNdxInterpolant2D called" << std::endl;

    Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);
    ;

    if (limits.vUp == limits.vMax)
    {
        return 0;
    }

    v = limits.vUp * exp(v * log(limits.vMax / limits.vUp));

    return integral.Integrate(
        limits.vUp, v, boost::bind(&Parametrization::FunctionToDNdxIntegral, parametrization_, energy, _1), 3, 1);
}

// ------------------------------------------------------------------------- //
double IonizInterpolant::CalculateStochasticLoss(double energy, double rnd1)
{
    double rnd, rsum;
    const Medium& medium = parametrization_->GetMedium();

    rnd  = medium.GetSumCharge() * rnd1;
    rsum = 0;

    for (unsigned int i = 0; i < components_.size(); i++)
    {
        rsum += components_[i]->GetAtomInMolecule() * components_[i]->GetNucCharge();

        if (rsum > rnd)
        {
            Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

            if (limits.vUp == limits.vMax)
            {
                return energy * limits.vUp;
            }
            return energy * (limits.vUp * exp(dndx_interpolant_2d_[0]->FindLimit(energy, rnd * sum_of_rates_) *
                                              log(limits.vMax / limits.vUp)));
        }
    }

    log_fatal("m.totZ was not initialized correctly");

    return 0;
}