/*! \file   ParticleDef.cxx
 *   \brief  Source file for definition of the ParticleDef class object.
 *
 *   For more details see the class documentation.
 *
 *   \date   Sat Aug  5 14:47:16 CEST 2017
 *   \author Mario Dunsch
 */

#include <iostream>
#include <sstream>
#include <string>
#include <map>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/decay/LeptonicDecayChannel.h"
#include "PROPOSAL/decay/ManyBodyPhaseSpace.h"
#include "PROPOSAL/decay/StableChannel.h"
#include "PROPOSAL/decay/TwoBodyPhaseSpace.h"
#include "PROPOSAL/particle/ParticleDef.h"

#include "PROPOSAL/methods.h"
#include "PROPOSAL/Logging.h"

#define PARTICLE_IMP(cls, MASS, LIFETIME, CHARGE, ID)                                                                  \
    cls##Def::cls##Def()                                                                                               \
        : ParticleDef(#cls,                                                                                            \
                      MASS,                                                                                            \
                      MASS,                                                                                            \
                      LIFETIME,                                                                                        \
                      CHARGE,                                                                                          \
                      HardComponentTables::EmptyTable,                                                                 \
                      DecayTable().addChannel(1.1, StableChannel()),                                                   \
                      ID)                                                                                              \
    {                                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    cls##Def::~cls##Def() {}

using namespace PROPOSAL;

/******************************************************************************
 *                            HardComponentTables                              *
 ******************************************************************************/

const HardComponentTables::VecType HardComponentTables::MuonTable = {
    { 7.174409e-4, 1.7132e-3, 4.082304e-3, 8.628455e-3, 0.01244159, 0.02204591, 0.03228755 },
    { -0.2436045, -0.5756682, -1.553973, -3.251305, -5.976818, -9.495636, -13.92918 },
    { -0.2942209, -0.68615, -2.004218, -3.999623, -6.855045, -10.05705, -14.37232 },
    { -0.1658391, -0.3825223, -1.207777, -2.33175, -3.88775, -5.636636, -8.418409 },
    { -0.05227727, -0.1196482, -0.4033373, -0.7614046, -1.270677, -1.883845, -2.948277 },
    { -9.328318e-3, -0.02124577, -0.07555636, -0.1402496, -0.2370768, -0.3614146, -0.5819409 },
    { -8.751909e-4, -1.987841e-3, -7.399682e-3, -0.01354059, -0.02325118, -0.03629659, -0.059275 },
    { -3.343145e-5, -7.584046e-5, -2.943396e-4, -5.3155e-4, -9.265136e-4, -1.473118e-3, -2.419946e-3 }
};

const HardComponentTables::VecType HardComponentTables::TauTable = {
    { -1.269205e-4, -2.843877e-4, -5.761546e-4, -1.195445e-3, -1.317386e-3, -9.689228e-15, -6.4595e-15 },
    { -0.01563032, -0.03589573, -0.07768545, -0.157375, -0.2720009, -0.4186136, -0.8045046 },
    { 0.04693954, 0.1162945, 0.3064255, 0.7041273, 1.440518, 2.533355, 3.217832 },
    { 0.05338546, 0.130975, 0.3410341, 0.7529364, 1.425927, 2.284968, 2.5487 },
    { 0.02240132, 0.05496, 0.144945, 0.3119032, 0.5576727, 0.8360727, 0.8085682 },
    { 4.658909e-3, 0.01146659, 0.03090286, 0.06514455, 0.1109868, 0.1589677, 0.1344223 },
    { 4.822364e-4, 1.193018e-3, 3.302773e-3, 6.843364e-3, 0.011191, 0.015614, 0.01173827 },
    { 1.9837e-5, 4.940182e-5, 1.409573e-4, 2.877909e-4, 4.544877e-4, 6.280818e-4, 4.281932e-4 }
};
const HardComponentTables::VecType HardComponentTables::EmptyTable;

/******************************************************************************
 *                                ParticleMap                                 *
 ******************************************************************************/

static std::map<const int,const ParticleDef&> ParticleDefs{
    {11,  EMinusDef::Get()},
    {-11, EPlusDef::Get()},
    {12,  NuEDef::Get()},
    {-12, NuEBarDef::Get()},
    {13,  MuMinusDef::Get()},
    {-13, MuPlusDef::Get()},
    {14,  NuMuDef::Get()},
    {-14, NuMuBarDef::Get()},
    {15,  TauMinusDef::Get()},
    {-15, TauPlusDef::Get()},
    {16,  NuTauDef::Get()},
    {-16, NuTauBarDef::Get()},
    {21,  GammaDef::Get()}
};

//ParticleDef* e = ParticleDefs.find(11)->second;

/******************************************************************************
 *                                ParticleDef                                 *
 ******************************************************************************/

std::ostream& PROPOSAL::operator<<(std::ostream& os, ParticleDef const& def)
{
    std::stringstream ss;
    ss << " ParticleDef (" << &def << ") ";

    os << Helper::Centered(60, ss.str()) << '\n';

    os << def.name << '\n';
    os << "Mass:"
       << "\t\t" << def.mass << '\n';
    os << "Low:"
       << "\t\t" << def.low << '\n';
    os << "Lifetime:"
       << "\t" << def.lifetime << '\n';
    os << "Charge:"
       << "\t\t" << def.charge << '\n';
    os << "HardComponentTables:" << '\n';

    for (unsigned int i = 0; i < def.hard_component_table.size(); ++i)
    {
        for (unsigned int j = 0; j < def.hard_component_table[i].size(); ++j)
        {
            os << def.hard_component_table[i][j] << "\t";
        }
        os << '\n';
    }

    os << Helper::Centered(60, "");
    return os;
}

/*bool IsStandard(std::string name){
    if(name.compare("EMinus")   == 0 || name.compare("EPlus")    == 0||
       name.compare("MuMinus")  == 0 || name.compare("MuPlus")   == 0||
       name.compare("TauMinus") == 0 || name.compare("TauPlus")  == 0||
       name.compare("NuE")      == 0 || name.compare("NuEBar")   == 0||
       name.compare("NuMu")     == 0 || name.compare("NuMuBar")  == 0||
       name.compare("NuTau")    == 0 || name.compare("NuTauBar") == 0||
       name.compare("Gamma")    == 0)
        return true;
    return false;
}*/

ParticleDef::ParticleDef()
    : name("")
    , mass(0.0)
    , low(0.0)
    , lifetime(0.0)
    , charge(0.0)
    , hard_component_table(HardComponentTables::EmptyTable)
    , decay_table()
    , particle_id(0)
    , weak_partner(0)
{
}

ParticleDef::ParticleDef(std::string name,
                         double mass,
                         double low,
                         double lifetime,
                         double charge,
                         const HardComponentTables::VecType& table,
                         const DecayTable& decay_table,
                         int ParticleId)
    : name(name)
    , mass(mass)
    , low(low)
    , lifetime(lifetime)
    , charge(charge)
    , hard_component_table(table)
    , decay_table(decay_table)
    , particle_id(ParticleId)
    , weak_partner(0)
{
    //ParticleDefs.insert(std::make_pair<const int,const ParticleDef&>(particle_id,*this));
    //if(!IsStandard(name))
    //ParticleDefs.insert({particle_id,*this});
}

ParticleDef::ParticleDef(std::string name,
                         double mass,
                         double low,
                         double lifetime,
                         double charge,
                         const HardComponentTables::VecType& table,
                         const DecayTable& decay_table,
                         int ParticleId,
                         int partner)
        : name(name)
        , mass(mass)
        , low(low)
        , lifetime(lifetime)
        , charge(charge)
        , hard_component_table(table)
        , decay_table(decay_table)
        , particle_id(ParticleId)
        , weak_partner(partner)
{
    //if(!IsStandard(name))
    //ParticleDefs.insert({particle_id,*this});
}
ParticleDef::~ParticleDef() {}

ParticleDef::ParticleDef(const ParticleDef& def)
    : name(def.name)
    , mass(def.mass)
    , low(def.low)
    , lifetime(def.lifetime)
    , charge(def.charge)
    , hard_component_table(def.hard_component_table)
    , decay_table(def.decay_table)
    , particle_id(def.particle_id)
    , weak_partner(def.weak_partner)
{
    //if(!IsStandard(name))
    //ParticleDefs.insert({particle_id,*this});
}

const ParticleDef* ParticleDef::GetWeakPartner() const {
    if(weak_partner == 0){
        log_fatal("WeakPartner not defined for particle %s.", name.c_str());
        return nullptr;
    }
    else{
        std::map<const int,const ParticleDef&>::iterator it = ParticleDefs.find(weak_partner);
        if(it->first==weak_partner)
            return &(it->second);
        log_fatal("WeakPartner (Id: %s) is not defined in PROPOSAL",weak_partner);
        return nullptr;
    }
}

// ParticleDef& ParticleDef::operator=(const ParticleDef& def)
// {
//     if (this != &def)
//     {
//       ParticleDef tmp(def);
//       swap(tmp);
//     }
//     return *this;
// }
//
// void ParticleDef::swap(ParticleDef& def)
// {
//     using std::swap;
//
//     swap(name, def.name);
//     swap(mass, def.mass);
//     swap(low, def.low);
//     swap(lifetime, def.lifetime);
//     swap(charge, def.charge);
//     swap(hard_component_table, def.hard_component_table);
//     swap(decay_table, def.decay_table);
// }

bool ParticleDef::operator==(const ParticleDef& def) const
{
    if (name != def.name)
    {
        return false;
    } else if (mass != def.mass)
    {
        return false;
    } else if (low != def.low)
    {
        return false;
    } else if (lifetime != def.lifetime)
    {
        return false;
    } else if (charge != def.charge)
    {
        return false;
    } else if (hard_component_table != def.hard_component_table)
    {
        return false;
    } else if (decay_table != def.decay_table)
    {
        return false;
    } else if (particle_id != def.particle_id)
    {
        return false;
    }else if (weak_partner != def.weak_partner)
    {
        return false;
    } else
    {
        return true;
    }
}

bool ParticleDef::operator!=(const ParticleDef& def) const
{
    return !(*this == def);
}

/******************************************************************************
 *                                  Builder                                    *
 ******************************************************************************/

ParticleDef::Builder::Builder()
    : name("")
    , mass(0)
    , low(0)
    , lifetime(-1)
    , charge(-1)
    , hard_component_table(&HardComponentTables::EmptyTable)
    , decay_table()
    , particle_id(0)
    , weak_partner(0)
{
}

// ------------------------------------------------------------------------- //
// Special Particle definitions
// ------------------------------------------------------------------------- //

EMinusDef::EMinusDef()
        : ParticleDef(
            "EMinus",
            ME,
            ME,
            STABLE_PARTICLE,
            -1.0,
            HardComponentTables::EmptyTable,
            DecayTable().addChannel(1.1, StableChannel()),
            11,
            12)
{
}

EMinusDef::~EMinusDef() {}

EPlusDef::EPlusDef()
        : ParticleDef(
            "EPlus",
            ME,
            ME,
            STABLE_PARTICLE,
            1.0,
            HardComponentTables::EmptyTable,
            DecayTable().addChannel(1.1, StableChannel()),
            -11,
            -12)
{
}

EPlusDef::~EPlusDef() {}

MuMinusDef::MuMinusDef()
    : ParticleDef(
            "MuMinus",
            MMU,
            MMU,
            LMU,
            -1.0,
            HardComponentTables::MuonTable,
            DecayTable().addChannel(1.0, LeptonicDecayChannelApprox(EMinusDef::Get(), NuMuDef::Get(), NuEBarDef::Get())),
            13,
            14)
{
}

MuMinusDef::~MuMinusDef() {}

MuPlusDef::MuPlusDef()
    : ParticleDef("MuPlus",
                  MMU,
                  MMU,
                  LMU,
                  1.0,
                  HardComponentTables::MuonTable,
                  DecayTable().addChannel(1.0, LeptonicDecayChannelApprox(EPlusDef::Get(), NuEDef::Get(), NuMuBarDef::Get())),
                  -13,
                  -14)
{
}

MuPlusDef::~MuPlusDef() {}

TauMinusDef::TauMinusDef()
    : ParticleDef("TauMinus",
                  MTAU,
                  MTAU,
                  LTAU,
                  -1.0,
                  HardComponentTables::TauTable,
                  DecayTable()
                      .addChannel(0.1737, LeptonicDecayChannel(MuMinusDef::Get(), NuTauDef::Get(), NuMuBarDef::Get()))
                      .addChannel(0.1783, LeptonicDecayChannelApprox(EMinusDef::Get(), NuTauDef::Get(), NuEBarDef::Get()))
                      .addChannel(0.1153, TwoBodyPhaseSpace(PiMinusDef::Get(), NuTauDef::Get()))
                      .addChannel(0.2595,
                                  ManyBodyPhaseSpace::Builder()
                                      .addDaughter(PiMinusDef::Get())
                                      .addDaughter(Pi0Def::Get())
                                      .addDaughter(NuTauDef::Get())
                                      .build())
                      .addChannel(0.0952,
                                  ManyBodyPhaseSpace::Builder()
                                      .addDaughter(PiMinusDef::Get())
                                      .addDaughter(Pi0Def::Get())
                                      .addDaughter(Pi0Def::Get())
                                      .addDaughter(NuTauDef::Get())
                                      .build())
                      .addChannel(0.098,
                                  ManyBodyPhaseSpace::Builder()
                                      .addDaughter(PiMinusDef::Get())
                                      .addDaughter(PiMinusDef::Get())
                                      .addDaughter(PiPlusDef::Get())
                                      .addDaughter(NuTauDef::Get())
                                      .build())
                      .addChannel(0.0457,
                                  ManyBodyPhaseSpace::Builder()
                                      .addDaughter(PiMinusDef::Get())
                                      .addDaughter(PiMinusDef::Get())
                                      .addDaughter(PiPlusDef::Get())
                                      .addDaughter(Pi0Def::Get())
                                      .addDaughter(NuTauDef::Get())
                                      .build())
                      .addChannel(0.0119,
                                  ManyBodyPhaseSpace::Builder()
                                      .addDaughter(PiMinusDef::Get())
                                      .addDaughter(Pi0Def::Get())
                                      .addDaughter(Pi0Def::Get())
                                      .addDaughter(Pi0Def::Get())
                                      .addDaughter(NuTauDef::Get())
                                      .build())
                      .addChannel(0.01,
                                  ManyBodyPhaseSpace::Builder()
                                      .addDaughter(PiMinusDef::Get())
                                      .addDaughter(K0Def::Get())
                                      .addDaughter(NuTauDef::Get())
                                      .build())
                      .addChannel(0.00349,
                                  ManyBodyPhaseSpace::Builder()
                                      .addDaughter(KMinusDef::Get())
                                      .addDaughter(PiPlusDef::Get())
                                      .addDaughter(PiMinusDef::Get())
                                      .addDaughter(NuTauDef::Get())
                                      .build()),
                  15,
                  16)
{
}

TauMinusDef::~TauMinusDef() {}

TauPlusDef::TauPlusDef()
    : ParticleDef("TauPlus",
                  MTAU,
                  MTAU,
                  LTAU,
                  1.0,
                  HardComponentTables::TauTable,
                  DecayTable()
                      .addChannel(0.1737, LeptonicDecayChannel(MuPlusDef::Get(), NuTauBarDef::Get(), NuMuDef::Get()))
                      .addChannel(0.1783, LeptonicDecayChannelApprox(EPlusDef::Get(), NuTauBarDef::Get(), NuEDef::Get()))
                      .addChannel(0.1153, TwoBodyPhaseSpace(PiPlusDef::Get(), NuTauBarDef::Get()))
                      .addChannel(0.2595,
                                  ManyBodyPhaseSpace::Builder()
                                      .addDaughter(PiPlusDef::Get())
                                      .addDaughter(Pi0Def::Get())
                                      .addDaughter(NuTauBarDef::Get())
                                      .build())
                      .addChannel(0.0952,
                                  ManyBodyPhaseSpace::Builder()
                                      .addDaughter(PiPlusDef::Get())
                                      .addDaughter(Pi0Def::Get())
                                      .addDaughter(Pi0Def::Get())
                                      .addDaughter(NuTauBarDef::Get())
                                      .build())
                      .addChannel(0.098,
                                  ManyBodyPhaseSpace::Builder()
                                      .addDaughter(PiPlusDef::Get())
                                      .addDaughter(PiPlusDef::Get())
                                      .addDaughter(PiMinusDef::Get())
                                      .addDaughter(NuTauBarDef::Get())
                                      .build())
                      .addChannel(0.0457,
                                  ManyBodyPhaseSpace::Builder()
                                      .addDaughter(PiPlusDef::Get())
                                      .addDaughter(PiPlusDef::Get())
                                      .addDaughter(PiMinusDef::Get())
                                      .addDaughter(Pi0Def::Get())
                                      .addDaughter(NuTauBarDef::Get())
                                      .build())
                      .addChannel(0.0119,
                                  ManyBodyPhaseSpace::Builder()
                                      .addDaughter(PiPlusDef::Get())
                                      .addDaughter(Pi0Def::Get())
                                      .addDaughter(Pi0Def::Get())
                                      .addDaughter(Pi0Def::Get())
                                      .addDaughter(NuTauBarDef::Get())
                                      .build())
                      .addChannel(0.01,
                                  ManyBodyPhaseSpace::Builder()
                                      .addDaughter(PiPlusDef::Get())
                                      .addDaughter(K0Def::Get())
                                      .addDaughter(NuTauBarDef::Get())
                                      .build())
                      .addChannel(0.00349,
                                  ManyBodyPhaseSpace::Builder()
                                      .addDaughter(KPlusDef::Get())
                                      .addDaughter(PiPlusDef::Get())
                                      .addDaughter(PiMinusDef::Get())
                                      .addDaughter(NuTauBarDef::Get())
                                      .build()),
                  -15,
                  -16)
{
}

TauPlusDef::~TauPlusDef() {}

NuEDef::NuEDef()
        : ParticleDef(
            "NuE",
            0,
            1,
            STABLE_PARTICLE,
            0.0,
            HardComponentTables::EmptyTable,
            DecayTable().addChannel(1.1, StableChannel()),
            12,
            11)
{
}

NuEDef::~NuEDef() {}

NuEBarDef::NuEBarDef()
        : ParticleDef(
            "NuEBar",
            0,
            1,
            STABLE_PARTICLE,
            0.0,
            HardComponentTables::EmptyTable,
            DecayTable().addChannel(1.1, StableChannel()),
            -12,
            -11)
{
}

NuEBarDef::~NuEBarDef() {}

NuMuDef::NuMuDef()
        : ParticleDef(
            "NuMu",
            0,
            1,
            STABLE_PARTICLE,
            0.0,
            HardComponentTables::EmptyTable,
            DecayTable().addChannel(1.1, StableChannel()),
            14,
            13)
{
}

NuMuDef::~NuMuDef() {}

NuMuBarDef::NuMuBarDef()
        : ParticleDef(
            "NuMuBar",
            0,
            1,
            STABLE_PARTICLE,
            0.0,
            HardComponentTables::EmptyTable,
            DecayTable().addChannel(1.1, StableChannel()),
            -14,
            -13)
{
}

NuMuBarDef::~NuMuBarDef() {}

NuTauDef::NuTauDef()
        : ParticleDef(
            "NuTau",
            0,
            1,
            STABLE_PARTICLE,
            0.0,
            HardComponentTables::EmptyTable,
            DecayTable().addChannel(1.1, StableChannel()),
            16,
            15)
{
}

NuTauDef::~NuTauDef() {}

NuTauBarDef::NuTauBarDef()
        : ParticleDef(
            "NuTauBar",
            0,
            1,
            STABLE_PARTICLE,
            0.0,
            HardComponentTables::EmptyTable,
            DecayTable().addChannel(1.1, StableChannel()),
            -16,
            -15)
{
}

NuTauBarDef::~NuTauBarDef() {}

// ------------------------------------------------------------------------- //
// Photon definition:
// e_low is set to twice the electron mass, since photons with a lower energy
// can not produce any leptons and are therefore irrelevant for us...
// ------------------------------------------------------------------------- //

GammaDef::GammaDef()
        : ParticleDef(
        "Gamma",
        0.,
        2*ME,
        STABLE_PARTICLE,
        0.0,
        HardComponentTables::EmptyTable,
        DecayTable().addChannel(1.1, StableChannel()),
        21)
{
}

GammaDef::~GammaDef() {}

// ------------------------------------------------------------------------- //
// Signature for following macro definitions:
//
// class name, mass, lifetime, charge
//
// ------------------------------------------------------------------------- //


PARTICLE_IMP(StauMinus, MSTAU, STABLE_PARTICLE, -1.0, 1000015)
PARTICLE_IMP(StauPlus, MSTAU, STABLE_PARTICLE, 1.0, -1000015)

PARTICLE_IMP(Pi0, MPI0, LPI0, 0.0, 111)
PARTICLE_IMP(PiMinus, MPI, LPI, -1.0, -211)
PARTICLE_IMP(PiPlus, MPI, LPI, 1.0, 211)

PARTICLE_IMP(K0, MKAON, -1.0, 0.0, 311)
PARTICLE_IMP(KMinus, MKAON, LKAON, -1.0, -321)
PARTICLE_IMP(KPlus, MKAON, LKAON, 1.0, 321)

PARTICLE_IMP(PMinus, MP, STABLE_PARTICLE, -1.0, -2212)
PARTICLE_IMP(PPlus, MP, STABLE_PARTICLE, 1.0, 2212)

//PARTICLE_IMP(NuE, 0.0, STABLE_PARTICLE, 0.0)
//PARTICLE_IMP(NuEBar, 0.0, STABLE_PARTICLE, 0.0)

//PARTICLE_IMP(NuMu, 0.0, STABLE_PARTICLE, 0.0)
//PARTICLE_IMP(NuMuBar, 0.0, STABLE_PARTICLE, 0.0)

//PARTICLE_IMP(NuTau, 0.0, STABLE_PARTICLE, 0.0)
//PARTICLE_IMP(NuTauBar, 0.0, STABLE_PARTICLE, 0.0)

PARTICLE_IMP(Monopole, MMON, STABLE_PARTICLE, CMON, -1)

PARTICLE_IMP(SMPMinus, MSMP, STABLE_PARTICLE, -1.0, -2)
PARTICLE_IMP(SMPPlus, MSMP, STABLE_PARTICLE, 1.0, -3)

#undef PARTICLE_IMP
