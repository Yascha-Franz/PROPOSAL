
#include "gtest/gtest.h"
#include <algorithm>

#include "PROPOSAL/PROPOSAL.h"
#include <string>
using namespace PROPOSAL;

std::string PATH_TO_TABLES = "~/.local/share/PROPOSAL/tables";
/* std::string PATH_TO_TABLES = ""; */

ParticleDef getParticleDef(const std::string& name)
{
    if (name == "MuMinus") {
        return MuMinusDef::Get();
    } else if (name == "TauMinus") {
        return TauMinusDef::Get();
    } else {
        return EMinusDef::Get();
    }
}

TEST(MixedComparison, Comparison_equal)
{
    ParticleDef mu_def = MuMinusDef::Get();
    Water medium;
    Sphere sphere;
    EnergyCutSettings ecuts;

    Sector::Definition sector_def;
    sector_def.location = Sector::ParticleLocation::InsideDetector;
    sector_def.SetMedium(medium);
    sector_def.SetGeometry(sphere);
    sector_def.scattering_model = ScatteringFactory::Moliere;
    sector_def.cut_settings = ecuts;
    sector_def.do_continuous_randomization = true;

    SectorMixed sector(mu_def, sector_def);
    SectorMixed sector_2 = SectorMixed(mu_def, sector_def);

    EXPECT_TRUE(sector == sector_2);
}

TEST(MixedComparison, Comparison_not_equal)
{
    ParticleDef mu_def = MuMinusDef::Get();
    ParticleDef tau_def = TauMinusDef::Get();
    Water medium1;
    Ice medium2;
    Sphere geometry1;
    Cylinder geometry2;
    EnergyCutSettings ecuts1(500, 0.05);
    EnergyCutSettings ecuts2(400, 0.005);

    Sector::Definition sector_def1;
    sector_def1.location = Sector::ParticleLocation::InsideDetector;
    sector_def1.SetMedium(medium1);
    sector_def1.SetGeometry(geometry1);
    sector_def1.scattering_model = ScatteringFactory::Moliere;
    sector_def1.cut_settings = ecuts1;
    sector_def1.do_continuous_randomization = true;
    sector_def1.do_exact_time_calculation = true;
    sector_def1.stopping_decay = true;

    Sector::Definition sector_def2 = sector_def1;
    sector_def2.location = Sector::ParticleLocation::InfrontDetector;

    Sector::Definition sector_def3 = sector_def1;
    sector_def3.SetMedium(medium2);

    Sector::Definition sector_def4 = sector_def1;
    sector_def4.SetGeometry(geometry2);

    Sector::Definition sector_def5 = sector_def1;
    sector_def5.scattering_model = ScatteringFactory::Highland;

    Sector::Definition sector_def6 = sector_def1;
    sector_def6.cut_settings = ecuts2;

    Sector::Definition sector_def7 = sector_def1;
    sector_def7.do_continuous_randomization = false;

    Sector::Definition sector_def8 = sector_def1;
    sector_def8.do_exact_time_calculation = false;

    Sector::Definition sector_def9 = sector_def1;
    sector_def9.stopping_decay = false;

    SectorMixed sector(mu_def, sector_def1);
    SectorMixed sector_1(tau_def, sector_def1);
    SectorMixed sector_2(mu_def, sector_def2);
    SectorMixed sector_3(mu_def, sector_def3);
    SectorMixed sector_4(mu_def, sector_def4);
    SectorMixed sector_5(mu_def, sector_def5);
    SectorMixed sector_6(mu_def, sector_def6);
    SectorMixed sector_7(mu_def, sector_def7);
    SectorMixed sector_8(mu_def, sector_def8);
    SectorMixed sector_9(mu_def, sector_def9);

    EXPECT_TRUE(sector != sector_1);
    EXPECT_TRUE(sector != sector_2);
    EXPECT_TRUE(sector != sector_3);
    EXPECT_TRUE(sector != sector_4);
    EXPECT_TRUE(sector != sector_5);
    EXPECT_TRUE(sector != sector_6);
    EXPECT_TRUE(sector != sector_7);
    EXPECT_TRUE(sector != sector_8);
    EXPECT_TRUE(sector != sector_9);
}

TEST(MixedAssignment, Copyconstructor)
{
    ParticleDef mu_def = MuMinusDef::Get();
    Water water(1.0);
    Sphere geometry(Vector3D(0, 0, 0), 1000, 0);
    EnergyCutSettings ecuts(500, 0.05);

    Sector::Definition sector_def;
    sector_def.location = Sector::ParticleLocation::InsideDetector;
    sector_def.SetMedium(water);
    sector_def.SetGeometry(geometry);
    sector_def.scattering_model = ScatteringFactory::Moliere;
    sector_def.cut_settings = ecuts;

    SectorMixed sector_1(mu_def, sector_def);
    SectorMixed sector_2 = sector_1;
    EXPECT_TRUE(sector_1 == sector_2);
}

TEST(MixedAssignment, Copyconstructor2)
{
    ParticleDef mu = MuMinusDef::Get();
    Water water(1.0);
    Sphere geometry(Vector3D(), 1000, 0);
    EnergyCutSettings ecuts(500, 0.05);

    Sector::Definition sector_def;
    sector_def.location = Sector::ParticleLocation::InsideDetector;
    sector_def.SetMedium(water);
    sector_def.SetGeometry(geometry);
    sector_def.scattering_model = ScatteringFactory::Moliere;
    sector_def.cut_settings = ecuts;

    SectorMixed sector_1(mu, sector_def);
    SectorMixed sector_2(sector_1);
    EXPECT_TRUE(sector_1 == sector_2);
}

TEST(MixedSector, Continuous)
{
    std::ifstream in;
    std::string filename = "bin/TestFiles/Sector_ContinousLoss.txt";
    in.open(filename.c_str());

    if (!in.good()) {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(1234);
    std::string particleName;
    std::string mediumName;
    double ecut, vcut;
    double energy, initial_energy;
    InterpolationDef inter_def;
    inter_def.path_to_tables = PATH_TO_TABLES;
    inter_def.path_to_tables_readonly = PATH_TO_TABLES;

    Sector::Definition sector_def;
    ParticleDef* particle = new ParticleDef(MuMinusDef::Get());
    Medium* medium = new Medium(sector_def.GetMedium());
    EnergyCutSettings* cuts = new EnergyCutSettings(sector_def.cut_settings);
    Sector* sector = new SectorMixed(*particle, sector_def, inter_def);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        if (particle->name != particleName || cuts->GetEcut() != ecut
            || cuts->GetVcut() != vcut || medium->GetName() != mediumName) {
            delete particle, sector, medium, cuts;

            cuts = new EnergyCutSettings(ecut, vcut);
            medium = MediumFactory::Get().CreateMedium(mediumName);
            particle = new ParticleDef(getParticleDef(particleName));

            sector_def.SetMedium(*medium);
            sector_def.cut_settings = *cuts;

            sector = new SectorMixed(*particle, sector_def, inter_def);
        }

        while (ss >> energy) {
            double rndd = RandomGenerator::Get().RandomDouble();
            double decay_energy = sector->CalculateDecay(initial_energy, rndd);
            double rndi = RandomGenerator::Get().RandomDouble();
            double inter_energy
                = sector->CalculateInteraction(initial_energy, rndi);
            double energy_calc = std::max(decay_energy, inter_energy);
            ASSERT_NEAR(energy_calc, energy, std::abs(1e-3 * energy_calc));
            initial_energy = energy;
        }
    }
}

TEST(MixedSector, Stochastic)
{
    std::ifstream in;
    std::string filename = "bin/TestFiles/Sector_StochasticLoss.txt";
    in.open(filename.c_str());

    if (!in.good()) {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(1234);
    std::string particleName;
    std::string mediumName;
    double ecut, vcut;
    double energy, initial_energy, rnd;
    int interaction_type;
    InterpolationDef inter_def;
    inter_def.path_to_tables = PATH_TO_TABLES;
    inter_def.path_to_tables_readonly = PATH_TO_TABLES;

    Sector::Definition sector_def;
    ParticleDef* particle = new ParticleDef(MuMinusDef::Get());
    Medium* medium = new Medium(sector_def.GetMedium());
    EnergyCutSettings* cuts = new EnergyCutSettings(sector_def.cut_settings);
    Sector* sector = new SectorMixed(*particle, sector_def, inter_def);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        if (particle->name != particleName || cuts->GetEcut() != ecut
            || cuts->GetVcut() != vcut || medium->GetName() != mediumName) {
            delete particle, sector, medium, cuts;

            cuts = new EnergyCutSettings(ecut, vcut);
            medium = MediumFactory::Get().CreateMedium(mediumName);
            particle = new ParticleDef(getParticleDef(particleName));

            sector_def.SetMedium(*medium);
            sector_def.cut_settings = *cuts;

            sector = new SectorMixed(*particle, sector_def, inter_def);
        }

        while (ss >> energy >> interaction_type >> rnd) {
            std::pair<double, int> loss = sector->MakeStochasticLoss(initial_energy);
            double energy_calc = initial_energy - loss.first;
            double random = RandomGenerator::Get().RandomDouble();
            ASSERT_NEAR(random, rnd, std::abs(1e-3 * energy_calc));
            ASSERT_NEAR(loss.second, interaction_type, std::abs(1e-3 * energy_calc));
            ASSERT_NEAR(energy_calc, energy, std::abs(1e-3 * energy_calc));
            initial_energy = energy;
        }
    }
}

TEST(MixedSector, EnergyDisplacement)
{
    std::ifstream in;
    std::string filename = "bin/TestFiles/Sector_Energy_Distance.txt";
    in.open(filename.c_str());

    if (!in.good()) {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(1234);
    std::string particleName;
    std::string mediumName;
    double ecut, vcut;
    double displacement, energy, initial_energy;
    InterpolationDef inter_def;
    inter_def.path_to_tables = PATH_TO_TABLES;
    inter_def.path_to_tables_readonly = PATH_TO_TABLES;

    Sector::Definition sector_def;
    ParticleDef* particle = new ParticleDef(MuMinusDef::Get());
    Medium* medium = new Medium(sector_def.GetMedium());
    EnergyCutSettings* cuts = new EnergyCutSettings(sector_def.cut_settings);
    Sector* sector = new SectorMixed(*particle, sector_def, inter_def);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        if (particle->name != particleName || cuts->GetEcut() != ecut
            || cuts->GetVcut() != vcut || medium->GetName() != mediumName) {
            delete particle, sector, medium, cuts;

            cuts = new EnergyCutSettings(ecut, vcut);
            medium = MediumFactory::Get().CreateMedium(mediumName);
            particle = new ParticleDef(getParticleDef(particleName));

            sector_def.SetMedium(*medium);
            sector_def.cut_settings = *cuts;

            sector = new SectorMixed(*particle, sector_def, inter_def);
        }

        while (ss >> displacement >> energy) {
            double energy_calc
                = sector->CalculateDistance(initial_energy, displacement);
            ASSERT_NEAR(energy_calc, energy, std::abs(1e-3 * energy_calc));
        }
    }
}

TEST(MixedSector, Propagate)
{
    std::ifstream in;
    std::string filename = "bin/TestFiles/Sector_Propagate.txt";
    in.open(filename.c_str());

    if (!in.good()) {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(1234);
    std::string particleName;
    std::string mediumName;
    double ecut, vcut;
    double energy, initial_energy;
    int produced_particle;
    InterpolationDef inter_def;
    inter_def.path_to_tables = PATH_TO_TABLES;
    inter_def.path_to_tables_readonly = PATH_TO_TABLES;

    Sector::Definition sector_def;
    ParticleDef* particle = new ParticleDef(MuMinusDef::Get());
    Medium* medium = new Medium(sector_def.GetMedium());
    EnergyCutSettings* cuts = new EnergyCutSettings(sector_def.cut_settings);
    Sector* sector = new SectorMixed(*particle, sector_def, inter_def);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        if (particle->name != particleName || cuts->GetEcut() != ecut
            || cuts->GetVcut() != vcut || medium->GetName() != mediumName) {
            delete particle, sector, medium, cuts;

            cuts = new EnergyCutSettings(ecut, vcut);
            medium = MediumFactory::Get().CreateMedium(mediumName);
            particle = new ParticleDef(getParticleDef(particleName));

            sector_def.SetMedium(*medium);
            sector_def.cut_settings = *cuts;

            sector = new SectorMixed(*particle, sector_def, inter_def);
        }

        while (ss >> produced_particle >> energy) {
            DynamicData p_condition;
            p_condition.SetDirection(Vector3D(0, 0, -1));
            p_condition.SetPosition(Vector3D(0, 0, 0));
            p_condition.SetEnergy(initial_energy);

            Secondaries secondaries = sector->Propagate(p_condition, 1000, 0);
            int produced_particle = secondaries.GetNumberOfParticles();

            double energy_calc = secondaries.GetSecondaries().back().GetEnergy();
            DynamicData last_condition = secondaries.GetSecondaries().back();
            if (last_condition.GetTypeId() != static_cast<int>(InteractionType::Decay)) {
                ASSERT_NEAR(last_condition.GetEnergy(), energy, std::abs(1e-3 * energy_calc));
            } else {
                ASSERT_NEAR(-last_condition.GetPropagatedDistance(), energy, std::abs(1e-3 * energy_calc));
            }
            /* double energy_calc = secondaries.GetSecondaries().back().GetEnergy(); */
            /* ASSERT_NEAR(energy_calc, energy, std::abs(1e-3 * energy_calc)); */
        }
    }
    delete medium;
}

TEST(StochasticComparison, Comparison_equal)
{
    ParticleDef mu_def = MuMinusDef::Get();
    Water medium;
    Sphere sphere;
    EnergyCutSettings ecuts;

    Sector::Definition sector_def;
    sector_def.location = Sector::ParticleLocation::InsideDetector;
    sector_def.SetMedium(medium);
    sector_def.SetGeometry(sphere);
    sector_def.scattering_model = ScatteringFactory::Moliere;
    sector_def.cut_settings = ecuts;
    sector_def.do_continuous_randomization = true;

    SectorMixed sector(mu_def, sector_def);
    SectorMixed sector_2 = SectorMixed(mu_def, sector_def);

    EXPECT_TRUE(sector == sector_2);
}

TEST(StochasticComparison, Comparison_not_equal)
{
    ParticleDef mu_def = MuMinusDef::Get();
    ParticleDef tau_def = TauMinusDef::Get();
    Water medium1;
    Ice medium2;
    Sphere geometry1;
    Cylinder geometry2;
    EnergyCutSettings ecuts1(500, 0.05);
    EnergyCutSettings ecuts2(400, 0.005);

    Sector::Definition sector_def1;
    sector_def1.location = Sector::ParticleLocation::InsideDetector;
    sector_def1.SetMedium(medium1);
    sector_def1.SetGeometry(geometry1);
    sector_def1.scattering_model = ScatteringFactory::Moliere;
    sector_def1.cut_settings = ecuts1;
    sector_def1.do_continuous_randomization = true;
    sector_def1.do_exact_time_calculation = true;
    sector_def1.stopping_decay = true;

    Sector::Definition sector_def2 = sector_def1;
    sector_def2.location = Sector::ParticleLocation::InfrontDetector;

    Sector::Definition sector_def3 = sector_def1;
    sector_def3.SetMedium(medium2);

    Sector::Definition sector_def4 = sector_def1;
    sector_def4.SetGeometry(geometry2);

    Sector::Definition sector_def5 = sector_def1;
    sector_def5.scattering_model = ScatteringFactory::Highland;

    Sector::Definition sector_def6 = sector_def1;
    sector_def6.cut_settings = ecuts2;

    Sector::Definition sector_def7 = sector_def1;
    sector_def7.do_continuous_randomization = false;

    Sector::Definition sector_def8 = sector_def1;
    sector_def8.do_exact_time_calculation = false;

    Sector::Definition sector_def9 = sector_def1;
    sector_def9.stopping_decay = false;

    SectorStochastic sector(mu_def, sector_def1);
    SectorStochastic sector_1(tau_def, sector_def1);
    SectorStochastic sector_2(mu_def, sector_def2);
    SectorStochastic sector_3(mu_def, sector_def3);
    SectorStochastic sector_4(mu_def, sector_def4);
    SectorStochastic sector_5(mu_def, sector_def5);
    SectorStochastic sector_6(mu_def, sector_def6);
    SectorStochastic sector_7(mu_def, sector_def7);
    SectorStochastic sector_8(mu_def, sector_def8);
    SectorStochastic sector_9(mu_def, sector_def9);

    EXPECT_TRUE(sector != sector_1);
    EXPECT_TRUE(sector != sector_2);
    EXPECT_TRUE(sector != sector_3);
    EXPECT_TRUE(sector != sector_4);
    EXPECT_TRUE(sector != sector_5);
    EXPECT_TRUE(sector != sector_6);
    EXPECT_TRUE(sector != sector_7);
    EXPECT_TRUE(sector != sector_8);
    EXPECT_TRUE(sector != sector_9);
}

TEST(StochasticAssignment, Copyconstructor)
{
    ParticleDef mu_def = MuMinusDef::Get();
    Water water(1.0);
    Sphere geometry(Vector3D(0, 0, 0), 1000, 0);
    EnergyCutSettings ecuts(500, 0.05);

    Sector::Definition sector_def;
    sector_def.location = Sector::ParticleLocation::InsideDetector;
    sector_def.SetMedium(water);
    sector_def.SetGeometry(geometry);
    sector_def.scattering_model = ScatteringFactory::Moliere;
    sector_def.cut_settings = ecuts;

    SectorStochastic sector_1(mu_def, sector_def);
    SectorStochastic sector_2 = sector_1;
    EXPECT_TRUE(sector_1 == sector_2);
}

TEST(StochasticAssignment, Copyconstructor2)
{
    ParticleDef mu = MuMinusDef::Get();
    Water water(1.0);
    Sphere geometry(Vector3D(), 1000, 0);
    EnergyCutSettings ecuts(500, 0.05);

    Sector::Definition sector_def;
    sector_def.location = Sector::ParticleLocation::InsideDetector;
    sector_def.SetMedium(water);
    sector_def.SetGeometry(geometry);
    sector_def.scattering_model = ScatteringFactory::Moliere;
    sector_def.cut_settings = ecuts;

    SectorStochastic sector_1(mu, sector_def);
    SectorStochastic sector_2(sector_1);
    EXPECT_TRUE(sector_1 == sector_2);
}

/*TEST(StochasticSector, Continuous)
{
    std::ifstream in;
    std::string filename = "bin/TestFiles/Sector_ContinousLoss.txt";
    in.open(filename.c_str());

    if (!in.good()) {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(1234);
    std::string particleName;
    std::string mediumName;
    double ecut, vcut;
    double energy, initial_energy;
    InterpolationDef inter_def;
    inter_def.path_to_tables = PATH_TO_TABLES;
    inter_def.path_to_tables_readonly = PATH_TO_TABLES;

    Sector::Definition sector_def;
    ParticleDef* particle = new ParticleDef(MuMinusDef::Get());
    Medium* medium = new Medium(sector_def.GetMedium());
    EnergyCutSettings* cuts = new EnergyCutSettings(sector_def.cut_settings);
    Sector* sector = new SectorStochastic(*particle, sector_def, inter_def);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        if (particle->name != particleName || cuts->GetEcut() != ecut
            || cuts->GetVcut() != vcut || medium->GetName() != mediumName) {
            delete particle, sector, medium, cuts;

            cuts = new EnergyCutSettings(ecut, vcut);
            medium = MediumFactory::Get().CreateMedium(mediumName);
            particle = new ParticleDef(getParticleDef(particleName));

            sector_def.SetMedium(*medium);
            sector_def.cut_settings = *cuts;

            sector = new SectorStochastic(*particle, sector_def, inter_def);
        }

        while (ss >> energy) {
            double rndd = RandomGenerator::Get().RandomDouble();
            double decay_energy = sector->CalculateDecay(initial_energy, rndd);
            double rndi = RandomGenerator::Get().RandomDouble();
            double inter_energy
                = sector->CalculateInteraction(initial_energy, rndi);
            double energy_calc = std::max(decay_energy, inter_energy);
            ASSERT_NEAR(energy_calc, energy, std::abs(1e-3 * energy_calc));
            initial_energy = energy;
        }
    }
}*/

TEST(StochasticSector, Stochastic)
{
    std::ifstream in;
    std::string filename = "bin/TestFiles/Sector_StochasticLoss.txt";
    in.open(filename.c_str());

    if (!in.good()) {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(1234);
    std::string particleName;
    std::string mediumName;
    double ecut, vcut;
    double energy, initial_energy, rnd;
    int interaction_type;
    InterpolationDef inter_def;
    inter_def.path_to_tables = PATH_TO_TABLES;
    inter_def.path_to_tables_readonly = PATH_TO_TABLES;

    Sector::Definition sector_def;
    ParticleDef* particle = new ParticleDef(MuMinusDef::Get());
    Medium* medium = new Medium(sector_def.GetMedium());
    EnergyCutSettings* cuts = new EnergyCutSettings(sector_def.cut_settings);
    Sector* sector = new SectorStochastic(*particle, sector_def, inter_def);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        if (particle->name != particleName || cuts->GetEcut() != ecut
            || cuts->GetVcut() != vcut || medium->GetName() != mediumName) {
            delete particle, sector, medium, cuts;

            cuts = new EnergyCutSettings(ecut, vcut);
            medium = MediumFactory::Get().CreateMedium(mediumName);
            particle = new ParticleDef(getParticleDef(particleName));

            sector_def.SetMedium(*medium);
            sector_def.cut_settings = *cuts;

            sector = new SectorStochastic(*particle, sector_def, inter_def);
        }

        while (ss >> energy >> interaction_type >> rnd) {
            std::pair<double, int> loss =sector->MakeStochasticLoss(initial_energy);
            double energy_calc = initial_energy - loss.first;
            double random = RandomGenerator::Get().RandomDouble();
            ASSERT_NEAR(random, rnd, std::abs(1e-3 * energy_calc));
            ASSERT_NEAR(loss.second, interaction_type, std::abs(1e-3 * energy_calc));
            ASSERT_NEAR(energy_calc, energy, std::abs(1e-3 * energy_calc));
            initial_energy = energy;
        }
    }
}

/*TEST(StochasticSector, EnergyDisplacement)
{
    std::ifstream in;
    std::string filename = "bin/TestFiles/Sector_Energy_Distance.txt";
    in.open(filename.c_str());

    if (!in.good()) {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(1234);
    std::string particleName;
    std::string mediumName;
    double ecut, vcut;
    double displacement, energy, initial_energy;
    InterpolationDef inter_def;
    inter_def.path_to_tables = PATH_TO_TABLES;
    inter_def.path_to_tables_readonly = PATH_TO_TABLES;

    Sector::Definition sector_def;
    ParticleDef* particle = new ParticleDef(MuMinusDef::Get());
    Medium* medium = new Medium(sector_def.GetMedium());
    EnergyCutSettings* cuts = new EnergyCutSettings(sector_def.cut_settings);
    Sector* sector = new SectorStochastic(*particle, sector_def, inter_def);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        if (particle->name != particleName || cuts->GetEcut() != ecut
            || cuts->GetVcut() != vcut || medium->GetName() != mediumName) {
            delete particle, sector, medium, cuts;

            cuts = new EnergyCutSettings(ecut, vcut);
            medium = MediumFactory::Get().CreateMedium(mediumName);
            particle = new ParticleDef(getParticleDef(particleName));

            sector_def.SetMedium(*medium);
            sector_def.cut_settings = *cuts;

            sector = new SectorStochastic(*particle, sector_def, inter_def);
        }

        while (ss >> displacement >> energy) {
            double energy_calc
                = sector->CalculateDistance(initial_energy, displacement);
            ASSERT_NEAR(energy_calc, energy, std::abs(1e-3 * energy_calc));
        }
    }
}*/

TEST(StochasticSector, Propagate)
{
    std::ifstream in;
    std::string filename = "bin/TestFiles/Sector_Propagate.txt";
    in.open(filename.c_str());

    if (!in.good()) {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(1234);
    std::string particleName;
    std::string mediumName;
    double ecut, vcut;
    double energy, initial_energy;
    int produced_particle;
    InterpolationDef inter_def;
    inter_def.path_to_tables = PATH_TO_TABLES;
    inter_def.path_to_tables_readonly = PATH_TO_TABLES;

    Sector::Definition sector_def;
    ParticleDef* particle = new ParticleDef(MuMinusDef::Get());
    Medium* medium = new Medium(sector_def.GetMedium());
    EnergyCutSettings* cuts = new EnergyCutSettings(sector_def.cut_settings);
    Sector* sector = new SectorStochastic(*particle, sector_def, inter_def);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        if (particle->name != particleName || cuts->GetEcut() != ecut
            || cuts->GetVcut() != vcut || medium->GetName() != mediumName) {
            delete particle, sector, medium, cuts;

            cuts = new EnergyCutSettings(ecut, vcut);
            medium = MediumFactory::Get().CreateMedium(mediumName);
            particle = new ParticleDef(getParticleDef(particleName));

            sector_def.SetMedium(*medium);
            sector_def.cut_settings = *cuts;

            sector = new SectorStochastic(*particle, sector_def, inter_def);
        }

        while (ss >> produced_particle >> energy) {
            DynamicData p_condition;
            p_condition.SetDirection(Vector3D(0, 0, -1));
            p_condition.SetPosition(Vector3D(0, 0, 0));
            p_condition.SetEnergy(initial_energy);

            Secondaries secondaries = sector->Propagate(p_condition, 1000, 0);
            int produced_particle = secondaries.GetNumberOfParticles();

            double energy_calc = secondaries.GetSecondaries().back().GetEnergy();
            DynamicData last_condition = secondaries.GetSecondaries().back();
            if (last_condition.GetTypeId() != static_cast<int>(InteractionType::Decay)) {
                ASSERT_NEAR(last_condition.GetEnergy(), energy, std::abs(1e-3 * energy_calc));
            } else {
                ASSERT_NEAR(-last_condition.GetPropagatedDistance(), energy, std::abs(1e-3 * energy_calc));
            }
            /* double energy_calc = secondaries.GetSecondaries().back().GetEnergy(); */
            /* ASSERT_NEAR(energy_calc, energy, std::abs(1e-3 * energy_calc)); */
        }
    }
    delete medium;
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
