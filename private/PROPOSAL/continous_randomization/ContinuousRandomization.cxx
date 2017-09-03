
// #include <algorithm>
#include <cmath>

#include <boost/bind.hpp>

#include "PROPOSAL/crossection/Bremsstrahlung.h"
#include "PROPOSAL/crossection/Epairproduction.h"
#include "PROPOSAL/crossection/Ionization.h"
#include "PROPOSAL/crossection/Photonuclear.h"
#include "PROPOSAL/continous_randomization/ContinuousRandomization.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"

using namespace std;
using namespace PROPOSAL;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double ContinuousRandomization::Randomize( const std::vector<CrossSections*>& cross_sections, double initial_energy, double final_energy, double rnd)
{
    double sigma,xhi,xlo,rndtmp;

//this happens if small distances are propagated and the
//energy loss is so small that it is smaller than the precision
//which is checked for in during the calculation.
    if(initial_energy == final_energy)
    {
        return final_energy;
    }

    sigma   =   sqrt( DE2de( cross_sections, initial_energy, final_energy) );

//It is not drawn from the real gaus distribution but rather from the
//area which is possible due to the limits of the initial energy and the
//particle mass. Another possibility would be to draw again but that would be
//more expensive.
//
//calculate the allowed region
    xhi     =   0.5+boost::math::erf((initial_energy        -final_energy)  /(SQRT2*sigma))/2;
    xlo     =   0.5+boost::math::erf((particle_def_.low   -final_energy)  /(SQRT2*sigma))/2;

//draw random number from the allowed region.
    rndtmp =  xlo + (xhi-xlo)*rnd;

//Calculate and return the needed value.
    return SQRT2*sigma*erfInv( 2*(rndtmp-0.5) )+final_energy;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-----------------------Enable and disable interpolation---------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ContinuousRandomization::EnableDE2dxInterpolation( const std::vector<CrossSections*>& cross_sections, std::string path, bool raw)
{
    if(do_dE2dx_Interpolation_)return;

    bool reading_worked =   true;
    bool storing_failed =   false;

    // charged anti leptons have the same cross sections like charged leptons
    // (except of diffractive Bremsstrahlung, where one can analyse the interference term if implemented)
    // so they use the same interpolation tables
    string particle_name = particle_def_.name;
    // string particle_name;
    // switch (particle_.GetType())
    // {
    //     case ParticleType::MuPlus:
    //         particle_name = PROPOSALParticle::GetName(ParticleType::MuMinus);
    //         break;
    //     case ParticleType::TauPlus:
    //         particle_name = PROPOSALParticle::GetName(ParticleType::TauMinus);
    //         break;
    //     case ParticleType::EPlus:
    //         particle_name = PROPOSALParticle::GetName(ParticleType::EMinus);
    //         break;
    //     case ParticleType::STauPlus:
    //         particle_name = PROPOSALParticle::GetName(ParticleType::STauMinus);
    //         break;
    //     default:
    //         particle_name = particle_def_.name;
    //         break;
    // }

    if(!path.empty())
    {
        stringstream filename;
        filename<<path<<"/Cont_dE2dx"
                <<"_"<<particle_name
                <<"_mass_"<<particle_def_.mass
                <<"_charge_"<<particle_def_.charge
                <<"_lifetime_"<<particle_def_.lifetime
                <<"_charge_"<<particle_def_.charge
                <<"_"<<cross_sections.at(0)->GetMedium()->GetName()
                <<"_"<<cross_sections.at(0)->GetMedium()->GetMassDensity();


        for(unsigned int i =0; i<cross_sections.size(); i++)
        {
            switch (cross_sections.at(i)->GetType())
            {
                case ParticleType::Brems:
                    filename << "_b"
                        << "_" << cross_sections.at(i)->GetParametrization()
                        << "_" << cross_sections.at(i)->GetLpmEffectEnabled();
                    break;
                case ParticleType::DeltaE:
                    filename << "_i";
                    break;
                case ParticleType::EPair:
                    filename << "_e"
                        << "_" << cross_sections.at(i)->GetLpmEffectEnabled();
                    break;
                case ParticleType::NuclInt:
                    filename << "_p"
                        << "_" << cross_sections.at(i)->GetParametrization();
                    break;
                default:
                    log_fatal("Unknown cross section");
                    exit(1);
            }
            filename<< "_" << cross_sections.at(i)->GetMultiplier()
                    << "_" << cross_sections.at(i)->GetEnergyCutSettings()->GetEcut()
                    << "_" << cross_sections.at(i)->GetEnergyCutSettings()->GetVcut();
        }

        if(!raw)
            filename<<".txt";

        if( FileExist(filename.str()) )
        {
            log_debug("Continuous Randomization parametrisation tables (dE2dx) will be read from file:\t%s",filename.str().c_str());
            ifstream input;

            if(raw)
            {
                input.open(filename.str().c_str(), ios::binary);
            }
            else
            {
                input.open(filename.str().c_str());
            }

            dE2dx_interpolant_ = new Interpolant();
            reading_worked = dE2dx_interpolant_->Load(input,raw);

            input.close();
        }
        if(!FileExist(filename.str()) || !reading_worked )
        {
            if(!reading_worked)
            {
                log_info("File %s is corrupted! Write it again!",filename.str().c_str());
            }

            log_info("Continuous Randomization parametrisation tables (dE2dx) will be saved to file:\t%s",filename.str().c_str());

            ofstream output;

            if(raw)
            {
                output.open(filename.str().c_str(), ios::binary);
            }
            else
            {
                output.open(filename.str().c_str());
            }
            if(output.good())
            {
                output.precision(16);

                dE2dx_interpolant_    =   new Interpolant(NUM2, particle_def_.low, BIGENERGY, boost::bind(&ContinuousRandomization::FunctionToBuildDE2dxInterplant, this,  boost::cref(cross_sections), _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);

                dE2dx_interpolant_->Save(output,raw);

            }
            else
            {
                storing_failed  =   true;
                log_warn("Can not open file %s for writing! Table will not be stored!",filename.str().c_str());
            }

            output.close();
        }
    }
    if(path.empty() || storing_failed)
    {
        dE2dx_interpolant_    =   new Interpolant(NUM2, particle_def_.low, BIGENERGY, boost::bind(&ContinuousRandomization::FunctionToBuildDE2dxInterplant, this,  boost::cref(cross_sections), _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);

    }

    do_dE2dx_Interpolation_=true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ContinuousRandomization::EnableDE2deInterpolation( const std::vector<CrossSections*>& cross_sections, std::string path, bool raw)
{
    if(do_dE2de_Interpolation_)return;

    bool reading_worked =   true;
    bool storing_failed =   false;

    // charged anti leptons have the same cross sections like charged leptons
    // (except of diffractive Bremsstrahlung, where one can analyse the interference term if implemented)
    // so they use the same interpolation tables
    string particle_name = particle_def_.name;
    // string particle_name;
    // switch (particle_.GetType())
    // {
    //     case ParticleType::MuPlus:
    //         particle_name = PROPOSALParticle::GetName(ParticleType::MuMinus);
    //         break;
    //     case ParticleType::TauPlus:
    //         particle_name = PROPOSALParticle::GetName(ParticleType::TauMinus);
    //         break;
    //     case ParticleType::EPlus:
    //         particle_name = PROPOSALParticle::GetName(ParticleType::EMinus);
    //         break;
    //     case ParticleType::STauPlus:
    //         particle_name = PROPOSALParticle::GetName(ParticleType::STauMinus);
    //         break;
    //     default:
    //         particle_name = particle_.GetName();
    //         break;
    // }

    if(!path.empty())
    {
        stringstream filename;
        filename<<path<<"/Cont_dE2de"
                <<"_"<<particle_name
                <<"_mass_"<<particle_def_.mass
                <<"_charge_"<<particle_def_.charge
                <<"_lifetime_"<<particle_def_.lifetime
                <<"_"<<cross_sections.at(0)->GetMedium()->GetName()
                <<"_"<<cross_sections.at(0)->GetMedium()->GetMassDensity();


        for(unsigned int i =0; i<cross_sections.size(); i++)
        {
            switch (cross_sections.at(i)->GetType())
            {
                case ParticleType::Brems:
                    filename << "_b"
                        << "_" << cross_sections.at(i)->GetParametrization()
                        << "_" << cross_sections.at(i)->GetLpmEffectEnabled();
                    break;
                case ParticleType::DeltaE:
                    filename << "_i";
                    break;
                case ParticleType::EPair:
                    filename << "_e"
                        << "_" << cross_sections.at(i)->GetLpmEffectEnabled();
                    break;
                case ParticleType::NuclInt:
                    filename << "_p"
                        << "_" << cross_sections.at(i)->GetParametrization();
                    break;
                default:
                    log_fatal("Unknown cross section");
                    exit(1);
            }
            filename<< "_" << cross_sections.at(i)->GetMultiplier()
                    << "_" << cross_sections.at(i)->GetEnergyCutSettings()->GetEcut()
                    << "_" << cross_sections.at(i)->GetEnergyCutSettings()->GetVcut();
        }

        if(!raw)
            filename<<".txt";

        if( FileExist(filename.str()) )
        {
            log_debug("Continuous Randomization parametrisation tables (dE2de) will be read from file:\t%s",filename.str().c_str());
            ifstream input;

            if(raw)
            {
                input.open(filename.str().c_str(), ios::binary);
            }
            else
            {
                input.open(filename.str().c_str());
            }
            dE2de_interpolant_ = new Interpolant();
            reading_worked = dE2de_interpolant_->Load(input,raw);

            dE2de_interpolant_diff_ = new Interpolant();
            reading_worked = dE2de_interpolant_diff_->Load(input,raw);

            input.close();
        }
        if(!FileExist(filename.str()) || !reading_worked )
        {
            if(!reading_worked)
            {
                log_info("File %s is corrupted! Write it again!",filename.str().c_str());
            }

            log_info("Continuous Randomization parametrisation tables (dE2de) will be saved to file:\t%s",filename.str().c_str());

            ofstream output;
            if(raw)
            {
                output.open(filename.str().c_str(), ios::binary);
            }
            else
            {
                output.open(filename.str().c_str());
            }
            if(output.good())
            {
                output.precision(16);

                dE2de_interpolant_       =   new Interpolant(NUM2, particle_def_.low, BIGENERGY, boost::bind(&ContinuousRandomization::FunctionToBuildDE2deInterplant, this,  boost::cref(cross_sections), _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);
                dE2de_interpolant_diff_  =   new Interpolant(NUM2, particle_def_.low, BIGENERGY, boost::bind(&ContinuousRandomization::FunctionToBuildDE2deInterplantDiff, this,  boost::cref(cross_sections), _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);

                dE2de_interpolant_->Save(output,raw);
                dE2de_interpolant_diff_->Save(output,raw);
            }
            else
            {
                storing_failed  =   true;
                log_warn("Can not open file %s for writing! Table will not be stored!",filename.str().c_str());
            }

            output.close();
        }
    }
    if(path.empty() || storing_failed)
    {
        dE2de_interpolant_       =   new Interpolant(NUM2, particle_def_.low, BIGENERGY, boost::bind(&ContinuousRandomization::FunctionToBuildDE2deInterplant, this,  boost::cref(cross_sections), _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);
        dE2de_interpolant_diff_  =   new Interpolant(NUM2, particle_def_.low, BIGENERGY, boost::bind(&ContinuousRandomization::FunctionToBuildDE2deInterplantDiff, this,  boost::cref(cross_sections), _1), order_of_interpolation_, false, false, true, order_of_interpolation_, false, false, false);
    }

    do_dE2de_Interpolation_=true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

void ContinuousRandomization::DisableDE2dxInterpolation()
{
    delete dE2dx_interpolant_;
    do_dE2dx_Interpolation_=false;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ContinuousRandomization::DisableDE2deInterpolation()
{
    delete dE2de_interpolant_;
    delete dE2de_interpolant_diff_;
    do_dE2de_Interpolation_=false;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


ContinuousRandomization::ContinuousRandomization(PROPOSALParticle& particle)
    : particle_def_(particle.GetParticleDef())
    , do_dE2dx_Interpolation_    ( false )
    ,do_dE2de_Interpolation_    ( false )
    ,which_cross_               ( 0 )
    ,order_of_interpolation_    ( 5 )
{
    dE2dx_integral_     =   new Integral(IROMB, IMAXS, IPREC);
    dE2de_integral_     =   new Integral(IROMB, IMAXS, IPREC2);

    dE2dx_interpolant_      =   NULL;
    dE2de_interpolant_      =   NULL;
    dE2de_interpolant_diff_ =   NULL;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// ContinuousRandomization::ContinuousRandomization(PROPOSALParticle*  Medium* medium, std::vector<CrossSections*> cross_sections)
//     :particle_                  ( particle )
//     ,medium_                    ( medium )
//     ,cross_sections_            ( cross_sections )
//     ,do_dE2dx_Interpolation_    ( false )
//     ,do_dE2de_Interpolation_    ( false )
//     ,which_cross_               ( 0 )
//     ,order_of_interpolation_    ( 5 )
// {
//     dE2dx_integral_     =   new Integral(IROMB, IMAXS, IPREC);
//     dE2de_integral_     =   new Integral(IROMB, IMAXS, IPREC2);
//
//     dE2dx_interpolant_      =   NULL;
//     dE2de_interpolant_      =   NULL;
//     dE2de_interpolant_diff_ =   NULL;
//
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

// ContinuousRandomization::ContinuousRandomization(const ContinuousRandomization &continuous_randomization)
//     :do_dE2dx_Interpolation_    ( continuous_randomization.do_dE2dx_Interpolation_)
//     ,do_dE2de_Interpolation_    ( continuous_randomization.do_dE2de_Interpolation_)
//     ,dE2dx_integral_            ( new Integral(*continuous_randomization.dE2dx_integral_) )
//     ,dE2de_integral_            ( new Integral(*continuous_randomization.dE2de_integral_) )
//     ,which_cross_               ( continuous_randomization.which_cross_ )
//     ,order_of_interpolation_    ( continuous_randomization.order_of_interpolation_ )
// {
//     cross_sections_.resize(continuous_randomization.cross_sections_.size());
//
//     for(unsigned int i =0; i<continuous_randomization.cross_sections_.size(); i++)
//     {
//         switch (continuous_randomization.cross_sections_.at(i)->GetType())
//         {
//             case ParticleType::Brems:
//                 cross_sections_.at(i) = new Bremsstrahlung( *(Bremsstrahlung*)continuous_randomization.cross_sections_.at(i) );
//                 break;
//             case ParticleType::DeltaE:
//                 cross_sections_.at(i) = new Ionization( *(Ionization*)continuous_randomization.cross_sections_.at(i) );
//                 break;
//             case ParticleType::EPair:
//                 cross_sections_.at(i) = new Epairproduction( *(Epairproduction*)continuous_randomization.cross_sections_.at(i) );
//                 break;
//             case ParticleType::NuclInt:
//                 cross_sections_.at(i) = new Photonuclear( *(Photonuclear*)continuous_randomization.cross_sections_.at(i) );
//                 break;
//             default:
//                 log_fatal("Unknown cross section of type '%i' ", continuous_randomization.cross_sections_.at(i)->GetType());
//         }
//     }
//
//     if(continuous_randomization.dE2dx_interpolant_ != NULL)
//     {
//         dE2dx_interpolant_ = new Interpolant(*continuous_randomization.dE2dx_interpolant_) ;
//     }
//     else
//     {
//         dE2dx_interpolant_ = NULL;
//     }
//
//     if(continuous_randomization.dE2de_interpolant_ != NULL)
//     {
//         dE2de_interpolant_ = new Interpolant(*continuous_randomization.dE2de_interpolant_) ;
//     }
//     else
//     {
//         dE2de_interpolant_ = NULL;
//     }
//
//     if(continuous_randomization.dE2de_interpolant_diff_ != NULL)
//     {
//         dE2de_interpolant_diff_ = new Interpolant(*continuous_randomization.dE2de_interpolant_diff_) ;
//     }
//     else
//     {
//         dE2de_interpolant_diff_ = NULL;
//     }
//
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// ContinuousRandomization& ContinuousRandomization::operator=(const ContinuousRandomization &continuous_randomization)
// {
//     if (this != &continuous_randomization)
//     {
//       ContinuousRandomization tmp(continuous_randomization);
//       swap(tmp);
//     }
//     return *this;
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool ContinuousRandomization::operator==(const ContinuousRandomization &continuous_randomization) const
{
    if( *dE2dx_integral_        != *continuous_randomization.dE2dx_integral_ )          return false;
    if( *dE2de_integral_        != *continuous_randomization.dE2de_integral_ )          return false;
    if( do_dE2dx_Interpolation_ != continuous_randomization.do_dE2dx_Interpolation_ )   return false;
    if( do_dE2de_Interpolation_ != continuous_randomization.do_dE2de_Interpolation_ )   return false;
    if( which_cross_            != continuous_randomization.which_cross_)               return false;
    if( order_of_interpolation_ != continuous_randomization.order_of_interpolation_ )   return false;


    // for(unsigned int i =0; i<continuous_randomization.cross_sections_.size(); i++)
    // {
    //     switch (continuous_randomization.cross_sections_.at(i)->GetType())
    //     {
    //         case ParticleType::Brems:
    //             if( *(Bremsstrahlung*)cross_sections_.at(i) !=  *(Bremsstrahlung*)continuous_randomization.cross_sections_.at(i) ) return false;
    //             break;
    //         case ParticleType::DeltaE:
    //             if( *(Ionization*)cross_sections_.at(i) != *(Ionization*)continuous_randomization.cross_sections_.at(i) ) return false;
    //             break;
    //         case ParticleType::EPair:
    //             if( *(Epairproduction*)cross_sections_.at(i) !=  *(Epairproduction*)continuous_randomization.cross_sections_.at(i) ) return false;
    //             break;
    //         case ParticleType::NuclInt:
    //             if( *(Photonuclear*)cross_sections_.at(i) !=  *(Photonuclear*)continuous_randomization.cross_sections_.at(i) )  return false;
    //             break;
    //         default:
    //             log_fatal("Unknown cross section of type '%i' ", continuous_randomization.cross_sections_.at(i)->GetType());
    //             return false;
    //     }
    // }

    if( dE2dx_interpolant_ != NULL && continuous_randomization.dE2dx_interpolant_ != NULL)
    {
        if( *dE2dx_interpolant_   != *continuous_randomization.dE2dx_interpolant_)                      return false;
    }
    else if( dE2dx_interpolant_ != continuous_randomization.dE2dx_interpolant_)                         return false;

    if( dE2de_interpolant_ != NULL && continuous_randomization.dE2de_interpolant_ != NULL)
    {
        if( *dE2de_interpolant_   != *continuous_randomization.dE2de_interpolant_)                      return false;
    }
    else if( dE2de_interpolant_ != continuous_randomization.dE2de_interpolant_)                         return false;

    if( dE2de_interpolant_diff_ != NULL && continuous_randomization.dE2de_interpolant_diff_ != NULL)
    {
        if( *dE2de_interpolant_diff_   != *continuous_randomization.dE2de_interpolant_diff_)            return false;
    }
    else if( dE2de_interpolant_diff_ != continuous_randomization.dE2de_interpolant_diff_)               return false;

    //else
    return true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool ContinuousRandomization::operator!=(const ContinuousRandomization &continuous_randomization) const
{
  return !(*this == continuous_randomization);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ContinuousRandomization::swap(ContinuousRandomization &continuous_randomization)
{
    using std::swap;

    dE2dx_integral_->swap( *continuous_randomization.dE2dx_integral_ );
    dE2de_integral_->swap( *continuous_randomization.dE2de_integral_ );

    swap( which_cross_ , continuous_randomization.which_cross_ );
    swap( do_dE2dx_Interpolation_, continuous_randomization.do_dE2dx_Interpolation_);
    swap( do_dE2de_Interpolation_, continuous_randomization.do_dE2de_Interpolation_);
    swap( order_of_interpolation_, continuous_randomization.order_of_interpolation_);

    // Set pointers again (to many swapping above....)
    //TODO(mario): hack Thu 2017/08/24
    // SetParticle( new PROPOSALParticle(tmp_particle1) );
    // continuous_randomization.SetParticle( new PROPOSALParticle(tmp_particle2) );

    if( dE2dx_interpolant_ != NULL && continuous_randomization.dE2dx_interpolant_ != NULL)
    {
        dE2dx_interpolant_->swap(*continuous_randomization.dE2dx_interpolant_);
    }
    else if( dE2dx_interpolant_ == NULL && continuous_randomization.dE2dx_interpolant_ != NULL)
    {
        dE2dx_interpolant_ = new Interpolant(*continuous_randomization.dE2dx_interpolant_);
        continuous_randomization.dE2dx_interpolant_ = NULL;
    }
    else if( dE2dx_interpolant_ != NULL && continuous_randomization.dE2dx_interpolant_ == NULL)
    {
        continuous_randomization.dE2dx_interpolant_ = new Interpolant(*dE2dx_interpolant_);
        dE2dx_interpolant_ = NULL;
    }

    if( dE2de_interpolant_ != NULL && continuous_randomization.dE2de_interpolant_ != NULL)
    {
        dE2de_interpolant_->swap(*continuous_randomization.dE2de_interpolant_);
    }
    else if( dE2de_interpolant_ == NULL && continuous_randomization.dE2de_interpolant_ != NULL)
    {
        dE2de_interpolant_ = new Interpolant(*continuous_randomization.dE2de_interpolant_);
        continuous_randomization.dE2de_interpolant_ = NULL;
    }
    else if( dE2de_interpolant_ != NULL && continuous_randomization.dE2de_interpolant_ == NULL)
    {
        continuous_randomization.dE2de_interpolant_ = new Interpolant(*dE2de_interpolant_);
        dE2de_interpolant_ = NULL;
    }

    if( dE2de_interpolant_diff_ != NULL && continuous_randomization.dE2de_interpolant_diff_ != NULL)
    {
        dE2de_interpolant_diff_->swap(*continuous_randomization.dE2de_interpolant_diff_);
    }
    else if( dE2de_interpolant_diff_ == NULL && continuous_randomization.dE2de_interpolant_diff_ != NULL)
    {
        dE2de_interpolant_diff_ = new Interpolant(*continuous_randomization.dE2de_interpolant_diff_);
        continuous_randomization.dE2de_interpolant_diff_ = NULL;
    }
    else if( dE2de_interpolant_diff_ != NULL && continuous_randomization.dE2de_interpolant_diff_ == NULL)
    {
        continuous_randomization.dE2de_interpolant_diff_ = new Interpolant(*dE2de_interpolant_diff_);
        dE2de_interpolant_diff_ = NULL;
    }
}

namespace PROPOSAL
{

ostream& operator<<(ostream& os, ContinuousRandomization const& continuous_randomization)
{
    os<<"---------------------------ContinuousRandomization( "<<&continuous_randomization<<" )---------------------------"<<endl;
    os<<fixed<<setprecision(5);
    os<<"\tdE2dx_integral:\t"<<continuous_randomization.dE2dx_integral_<<endl;
    os<<"\tdE2de_integral:\t"<<continuous_randomization.dE2de_integral_<<endl;
    os<<endl;
    os<<"\tdo_dE2dx_Interpolation:\t"<<continuous_randomization.do_dE2dx_Interpolation_<<endl;
    os<<"\tdE2dx_interpolant_:\t\t"<<continuous_randomization.dE2dx_interpolant_<<endl;
    os<<endl;
    os<<"\tdo_dE2de_Interpolation:\t"<<continuous_randomization.do_dE2de_Interpolation_<<endl;
    os<<"\tdE2de_interpolant_:\t\t"<<continuous_randomization.dE2de_interpolant_<<endl;
    os<<"\tdE2de_interpolant_diff:\t"<<continuous_randomization.dE2de_interpolant_diff_<<endl;
    os<<endl<<endl;
    os<<"\tTemp. variables: " << endl;
    os<<"\t\twhich_cros:\t\t\t" << continuous_randomization.which_cross_ << endl;
    os<<"\t\torder_of_interpolation:\t" << continuous_randomization.order_of_interpolation_ << endl;
    os<<"-----------------------------------------------------------------------------------------------";
    return os;
}

}  // namespace PROPOSAL

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------private Memberfunctions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ContinuousRandomization::DE2dx(const std::vector<CrossSections*>& cross_sections, double energy)
{
    if(do_dE2dx_Interpolation_)
    {
        return max( dE2dx_interpolant_->Interpolate(energy) , 0.0 );
    }

    double sum = 0;

    for(unsigned int i = 0 ; i < cross_sections.size() ; i++)
    {
        which_cross_ = i;

        for(int j=0 ; j < cross_sections.at(i)->GetMedium()->GetNumComponents() ; j++)
        {
            CrossSections::IntegralLimits limits = cross_sections.at(i)->SetIntegralLimits(energy, j);

            sum +=  dE2dx_integral_->Integrate (limits.vMin, limits.vUp, boost::bind(&ContinuousRandomization::FunctionToDE2dxIntegral, this,  boost::cref(cross_sections), energy, _1) ,2);

            if(cross_sections.at(i)->GetType() == ParticleType::DeltaE)
            {
                break;
            }
        }
    }

    return energy*energy*sum;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ContinuousRandomization::DE2de( const std::vector<CrossSections*>& cross_sections, double ei, double ef )
{
    if( do_dE2de_Interpolation_ )
    {
        if( abs( ei-ef ) > abs(ei)*HALF_PRECISION )
        {
            double aux  =   dE2de_interpolant_->Interpolate( ei );
            double aux2 =   aux - dE2de_interpolant_->Interpolate( ef );

            if( abs(aux2) > abs(aux)*HALF_PRECISION )
            {
                return max(aux2, 0.0);
            }
        }

        else
        {
            return max( dE2de_interpolant_diff_->Interpolate( (ei+ef)/2 )*(ef-ei) , 0.0 );
        }

    }

    return dE2de_integral_->Integrate(ei, ef, boost::bind(&ContinuousRandomization::FunctionToDE2deIntegral, this,  boost::cref(cross_sections), _1), 4);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------Functions to interpolate---------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ContinuousRandomization::FunctionToBuildDE2dxInterplant( const std::vector<CrossSections*>& cross_sections, double energy)
{
    return DE2dx(cross_sections, energy);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ContinuousRandomization::FunctionToBuildDE2deInterplant( const std::vector<CrossSections*>& cross_sections, double energy)
{
    return dE2de_integral_->Integrate(energy, particle_def_.low, boost::bind(&ContinuousRandomization::FunctionToDE2deIntegral, this,  boost::cref(cross_sections), _1) , 4);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ContinuousRandomization::FunctionToBuildDE2deInterplantDiff( const std::vector<CrossSections*>& cross_sections, double energy)
{
    return FunctionToDE2deIntegral( cross_sections, energy);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Functions to integrate----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ContinuousRandomization::FunctionToDE2dxIntegral(const std::vector<CrossSections*>& cross_sections, double energy, double v)
{
    return v*v*cross_sections.at(which_cross_)->FunctionToDNdxIntegral(energy, v);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ContinuousRandomization::FunctionToDE2deIntegral( const std::vector<CrossSections*>& cross_sections, double energy)
{

    double result = 0.0;
    double aux;

    for(unsigned int i = 0 ; i<cross_sections.size() ; i++)
    {
        aux     =   cross_sections.at(i)->CalculatedEdx(energy);
        result  +=  aux;
    }
    return -1/result*DE2dx( cross_sections, energy);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------------Setter---------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// void ContinuousRandomization::SetMedium(Medium* medium)
// {
//     medium_ = medium;
//     for(unsigned int i = 0 ; i < cross_sections_.size() ; i++)
//     {
//         cross_sections_.at(i)->SetMedium(medium_);
//     }
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// void ContinuousRandomization::SetParticle(PROPOSALParticle* particle)
// {
//     particle_ = particle;
//     for(unsigned int i = 0 ; i < cross_sections_.size() ; i++)
//     {
//         cross_sections_.at(i)->SetParticle(particle);
//     }
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// void ContinuousRandomization::SetCrosssections(
//         std::vector<CrossSections*> crosssections) {
//     cross_sections_ = crosssections;
// }
//
// PROPOSALParticle *ContinuousRandomization::GetBackup_particle() const
// {
//     return backup_particle_;
// }
//
// void ContinuousRandomization::SetBackup_particle(PROPOSALParticle *backup_particle)
// {
//     backup_particle_ = backup_particle;
// }
//
// void ContinuousRandomization::RestoreBackup_particle()
// {
//     particle_ = backup_particle_;
// }