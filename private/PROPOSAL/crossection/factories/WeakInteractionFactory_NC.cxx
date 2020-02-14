
#include <algorithm>

#include "PROPOSAL/crossection/WeakIntegral_NC.h"
#include "PROPOSAL/crossection/WeakInterpolant_NC.h"
#include "PROPOSAL/crossection/factories/WeakInteractionFactory_NC.h"
#include "PROPOSAL/crossection/parametrization/WeakInteraction_NC.h"

#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

WeakInteractionFactory_NC::WeakInteractionFactory_NC()
        : weak_nc_map_str_()
        , weak_nc_map_enum_()
        , string_enum_()
{
    Register("weakcoopersarkarmertsch_nc", CooperSarkarMertsch_NC, &WeakCooperSarkarMertsch_NC::create);
    Register("none", None, nullptr);
}

WeakInteractionFactory_NC::~WeakInteractionFactory_NC()
{
    string_enum_.clear();
    weak_nc_map_str_.clear();
    weak_nc_map_enum_.clear();
}

// ------------------------------------------------------------------------- //
// Most general creator
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
CrossSection* WeakInteractionFactory_NC::CreateWeakInteraction_NC(const ParticleDef& particle_def,
                                                            const Medium& medium,
                                                            const Definition& def) const
{
    if(def.parametrization == WeakInteractionFactory_NC::Enum::None){
        log_fatal("Can't return Weakinteraction Crosssection if parametrization is None");
        return NULL;
    }

    Weak_NCMapEnum::const_iterator it = weak_nc_map_enum_.find(def.parametrization);

    if (it != weak_nc_map_enum_.end())
    {
        return new WeakIntegral_NC(*it->second(particle_def, medium, def.multiplier));
    } else
    {
        log_fatal("WeakInteraction %s not registered!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
CrossSection* WeakInteractionFactory_NC::CreateWeakInteraction_NC(const ParticleDef& particle_def,
                                                            const Medium& medium,
                                                            const Definition& def,
                                                            InterpolationDef interpolation_def) const
{
    if(def.parametrization == WeakInteractionFactory_NC::Enum::None){
        log_fatal("Can't return Weakinteraction Crosssection if parametrization is None");
        return NULL;
    }

    Weak_NCMapEnum::const_iterator it = weak_nc_map_enum_.find(def.parametrization);

    if (it != weak_nc_map_enum_.end())
    {
        return new WeakInterpolant_NC(*it->second(particle_def, medium, def.multiplier), interpolation_def);
    } else
    {
        log_fatal("WeakInteraction %s not registered!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
void WeakInteractionFactory_NC::Register(const std::string& name,
                                      Enum enum_t,
                                      RegisterFunction create)
{
    weak_nc_map_str_[name]    = create;
    weak_nc_map_enum_[enum_t] = create;
    string_enum_.insert(name, enum_t);
}

// ------------------------------------------------------------------------- //
WeakInteractionFactory_NC::Enum WeakInteractionFactory_NC::GetEnumFromString(const std::string& name)
{
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    auto& left = string_enum_.GetLeft();
    auto it = left.find(name_lower);
    if (it != left.end())
    {
        return it->second;
    } else
    {
        log_fatal("WeakInteraction %s not registered!", name.c_str());
        return Fail; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
std::string WeakInteractionFactory_NC::GetStringFromEnum(const WeakInteractionFactory_NC::Enum& enum_t)
{
    auto& right = string_enum_.GetRight();
    auto it = right.find(enum_t);
    if (it != right.end())
    {
        return it->second;
    } else
    {
        log_fatal("WeakInteraction %s not registered!", typeid(enum_t).name());
        return ""; // Just to prevent warnings
    }
}
