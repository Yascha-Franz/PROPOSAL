
#include <PROPOSAL/crossection/factories/PhotoPairFactory.h>
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/propagation_utility/PropagationUtility.h"

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"

using namespace PROPOSAL;

/******************************************************************************
 *                            Propagation utility                              *
 ******************************************************************************/

Utility::Definition::Definition()
    // : do_interpolation(true)
    // , interpolation_def()
    : brems_def()
    , compton_def()
    , photo_def()
    , epair_def()
    , ioniz_def()
    , mupair_def()
    , weak_def()
    , weaknc_def()
    , photopair_def()
    , annihilation_def()
{
}


bool Utility::Definition::operator==(
    const Utility::Definition& utility_def) const {
    if (brems_def != utility_def.brems_def)
        return false;
    else if (compton_def != utility_def.compton_def)
        return false;
    else if (photo_def != utility_def.photo_def)
        return false;
    else if (epair_def != utility_def.epair_def)
        return false;
    else if (ioniz_def != utility_def.ioniz_def)
        return false;
    else if (mupair_def != utility_def.mupair_def)
        return false;
    else if (weak_def != utility_def.weak_def)
        return false;
    else if (weaknc_def != utility_def.weaknc_def)
        return false;
    else if (photopair_def != utility_def.photopair_def)
        return false;
    else if (annihilation_def != utility_def.annihilation_def)
        return false;

    return true;
}

bool Utility::Definition::operator!=(
    const Utility::Definition& utility_def) const {
    return !(*this == utility_def);
}

Utility::Definition::~Definition() {}

// ------------------------------------------------------------------------- //
// Constructors
// ------------------------------------------------------------------------- //

Utility::Utility(const ParticleDef& particle_def,
                 const Medium& medium,
                 const EnergyCutSettings& cut_settings,
                 Definition utility_def)
    : particle_def_(particle_def)
    , medium_(medium.clone())
    , cut_settings_(cut_settings)
    , crosssections_()
{
    if(utility_def.brems_def.parametrization!=BremsstrahlungFactory::Enum::None) {
        crosssections_.push_back(BremsstrahlungFactory::Get().CreateBremsstrahlung(
                particle_def_, *medium_, cut_settings_, utility_def.brems_def));
    }

    if(utility_def.photo_def.parametrization!=PhotonuclearFactory::Enum::None) {
        crosssections_.push_back(PhotonuclearFactory::Get().CreatePhotonuclear(
                particle_def_, *medium_, cut_settings_,utility_def.photo_def));
    }

    if(utility_def.epair_def.parametrization!=EpairProductionFactory::Enum::None) {
        crosssections_.push_back(EpairProductionFactory::Get().CreateEpairProduction(
                particle_def_, *medium_, cut_settings_, utility_def.epair_def));
    }

    if(utility_def.ioniz_def.parametrization!=IonizationFactory::Enum::None) {
        crosssections_.push_back(IonizationFactory::Get().CreateIonization(
                particle_def_, *medium_, cut_settings_, utility_def.ioniz_def));
    }
    else{
        log_debug("No Ionization cross section chosen. For lepton propagation,Initialization may fail because no cross"
                  "section for small energies are available. You may have to enable Ionization or set a higher e_low"
                  "parameter for the particle.");
    }

    if(utility_def.annihilation_def.parametrization!=AnnihilationFactory::Enum::None) {
        crosssections_.push_back(AnnihilationFactory::Get().CreateAnnihilation(
                particle_def_, *medium_, utility_def.annihilation_def));
        log_debug("Annihilation enabled");
    }

    if(utility_def.mupair_def.parametrization!=MupairProductionFactory::Enum::None) {
        crosssections_.push_back(MupairProductionFactory::Get().CreateMupairProduction(
            particle_def_, *medium_, cut_settings_, utility_def.mupair_def));
        log_debug("Mupair Production enabled");
    }

    if(utility_def.weak_def.parametrization!=WeakInteractionFactory::Enum::None) {
        crosssections_.push_back(WeakInteractionFactory::Get().CreateWeakInteraction(
                    particle_def_, *medium_, utility_def.weak_def));
        log_debug("Weak Interaction enabled");
    }

    if(utility_def.weaknc_def.parametrization!=WeakInteractionFactory_NC::Enum::None) {
        crosssections_.push_back(WeakInteractionFactory_NC::Get().CreateWeakInteraction_NC(
                    particle_def_, *medium_, utility_def.weaknc_def));
        log_debug("Weak NC Interaction enabled");
    }

    // Photon interactions

    if(utility_def.compton_def.parametrization!=ComptonFactory::Enum::None) {
        crosssections_.push_back(ComptonFactory::Get().CreateCompton(
                particle_def_, *medium_, cut_settings_, utility_def.compton_def));
        log_debug("Compton enabled");
    }

    if(utility_def.photopair_def.parametrization!=PhotoPairFactory::Enum::None) {
        crosssections_.push_back(PhotoPairFactory::Get().CreatePhotoPair(
                particle_def_, *medium_, utility_def.photopair_def));
        log_debug("PhotoPairProduction enabled");
    }
}

Utility::Utility(const ParticleDef& particle_def,
                 const Medium& medium,
                 const EnergyCutSettings& cut_settings,
                 Definition utility_def,
                 const InterpolationDef& interpolation_def)
    : particle_def_(particle_def)
    , medium_(medium.clone())
    , cut_settings_(cut_settings)
    , crosssections_()
{
    if(utility_def.brems_def.parametrization!=BremsstrahlungFactory::Enum::None) {
        crosssections_.push_back(BremsstrahlungFactory::Get().CreateBremsstrahlung(
                particle_def_, *medium_, cut_settings_, utility_def.brems_def, interpolation_def));
    }

    if(utility_def.photo_def.parametrization!=PhotonuclearFactory::Enum::None) {
        crosssections_.push_back(PhotonuclearFactory::Get().CreatePhotonuclear(
                particle_def_, *medium_, cut_settings_, utility_def.photo_def, interpolation_def));
    }

    if(utility_def.epair_def.parametrization!=EpairProductionFactory::Enum::None) {
        crosssections_.push_back(EpairProductionFactory::Get().CreateEpairProduction(
                particle_def_, *medium_, cut_settings_, utility_def.epair_def, interpolation_def));
    }

    if(utility_def.ioniz_def.parametrization!=IonizationFactory::Enum::None) {
        crosssections_.push_back(IonizationFactory::Get().CreateIonization(
                particle_def_, *medium_, cut_settings_, utility_def.ioniz_def, interpolation_def));
    }else{
        log_debug("No Ionization cross section chosen. For lepton propagation,Initialization may fail because no cross"
                  "section for small energies are available. You may have to enable Ionization or set a higher e_low"
                  "parameter for the particle.");
    }

    if(utility_def.annihilation_def.parametrization!=AnnihilationFactory::Enum::None) {
        crosssections_.push_back(AnnihilationFactory::Get().CreateAnnihilation(
                particle_def_, *medium_, utility_def.annihilation_def, interpolation_def));
        log_debug("Annihilation enabled");
    }

    if(utility_def.mupair_def.parametrization!=MupairProductionFactory::Enum::None) {
        crosssections_.push_back(MupairProductionFactory::Get().CreateMupairProduction(
                    particle_def_, *medium_, cut_settings_, utility_def.mupair_def, interpolation_def));
        log_debug("Mupair Production enabled");
    }

    if(utility_def.weak_def.parametrization!=WeakInteractionFactory::Enum::None) {
        crosssections_.push_back(WeakInteractionFactory::Get().CreateWeakInteraction(
                    particle_def_, *medium_, utility_def.weak_def, interpolation_def));
        log_debug("Weak Interaction enabled");
    }

    if(utility_def.weaknc_def.parametrization!=WeakInteractionFactory_NC::Enum::None) {
        crosssections_.push_back(WeakInteractionFactory_NC::Get().CreateWeakInteraction_NC(
                    particle_def_, *medium_, utility_def.weaknc_def, interpolation_def));
        log_debug("Weak NC Interaction enabled");
    }

    // Photon interactions

    if(utility_def.compton_def.parametrization!=ComptonFactory::Enum::None) {
        crosssections_.push_back(ComptonFactory::Get().CreateCompton(
                particle_def_, *medium_, cut_settings_, utility_def.compton_def, interpolation_def));
        log_debug("Compton enabled");
    }

    if(utility_def.photopair_def.parametrization!=PhotoPairFactory::Enum::None) {
        crosssections_.push_back(PhotoPairFactory::Get().CreatePhotoPair(
                particle_def_, *medium_, utility_def.photopair_def, interpolation_def));
        log_debug("PhotoPairProduction enabled");
    }
}

Utility::Utility(const std::vector<CrossSection*>& crosssections) try
    : particle_def_(crosssections.at(0)->GetParametrization().GetParticleDef()),
      medium_(crosssections.at(0)->GetParametrization().GetMedium().clone()),
      cut_settings_(crosssections.at(0)->GetParametrization().GetEnergyCuts()) {
    for (std::vector<CrossSection*>::const_iterator it = crosssections.begin();
         it != crosssections.end(); ++it) {
        if ((*it)->GetParametrization().GetParticleDef() != particle_def_) {
            log_fatal(
                "Particle definition of the cross section must be equal!");
        }

        if ((*it)->GetParametrization().GetMedium() != *medium_) {
            log_fatal("Medium of the cross section must be equal!");
        }

        if ((*it)->GetParametrization().GetEnergyCuts() != cut_settings_) {
            log_fatal("Energy cuts of the cross section must be equal!");
        }

        crosssections_.push_back((*it)->clone());
    }
} catch (const std::out_of_range& e) {
    log_fatal("At least one cross section is needed for initializing utility.");
}

Utility::Utility(const Utility& collection)
    : particle_def_(collection.particle_def_),
      medium_(collection.medium_->clone()),
      cut_settings_(collection.cut_settings_),
      crosssections_(collection.crosssections_.size(), NULL) {
    for (unsigned int i = 0; i < crosssections_.size(); ++i) {
        crosssections_[i] = collection.crosssections_[i]->clone();
    }
}

Utility::~Utility() {
    delete medium_;

    for (std::vector<CrossSection*>::const_iterator iter =
             crosssections_.begin();
         iter != crosssections_.end(); ++iter) {
        delete *iter;
    }

    crosssections_.clear();
}

bool Utility::operator==(const Utility& utility) const {
    if (particle_def_ != utility.particle_def_)
        return false;
    else if (*medium_ != *utility.medium_)
        return false;
    else if (cut_settings_ != utility.cut_settings_)
        return false;
    else if (crosssections_.size() != utility.crosssections_.size())
        return false;

    for (unsigned int i = 0; i < crosssections_.size(); ++i) {
        if (*crosssections_[i] != *utility.crosssections_[i])
            return false;
    }

    return true;
}

bool Utility::operator!=(const Utility& utility) const {
    return !(*this == utility);
}

/******************************************************************************
 *                            Utility Decorator                            *
 ******************************************************************************/

UtilityDecorator::UtilityDecorator(const Utility& utility)
    : utility_(*utility.clone()) {}

UtilityDecorator::UtilityDecorator(const UtilityDecorator& decorator)
    : utility_(*decorator.utility_.clone()) {}

UtilityDecorator::~UtilityDecorator() {}

bool UtilityDecorator::operator==(
    const UtilityDecorator& utility_decorator) const {
    if (typeid(*this) != typeid(utility_decorator))
        return false;
    if (utility_ != utility_decorator.utility_)
        return false;
    else
        return this->compare(utility_decorator);
}

bool UtilityDecorator::operator!=(
    const UtilityDecorator& utility_decorator) const {
    return !(*this == utility_decorator);
}

// ------------------------------------------------------------------------- //
double UtilityDecorator::FunctionToIntegral(double energy) {
    double result = 0.0;

    const std::vector<CrossSection*> crosssections =
        utility_.GetCrosssections();

    for (std::vector<CrossSection*>::const_iterator iter =
             crosssections.begin();
         iter != crosssections.end(); ++iter) {
        result += (*iter)->CalculatedEdx(energy);
    }

    return -1.0 / result;
}

double UtilityDecorator::Calculate(double ei,
                                   double ef,
                                   double rnd,
                                   Vector3D xi,
                                   Vector3D direction) {
    (void)xi;
    (void)direction;

    return this->Calculate(ei, ef, rnd);
}
