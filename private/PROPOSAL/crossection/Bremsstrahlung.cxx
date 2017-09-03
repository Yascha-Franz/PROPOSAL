
// #include <iostream>
// #include <cmath>
// #include <algorithm>
// #include <sstream>

#include <boost/bind.hpp>

#include "PROPOSAL/crossection/Bremsstrahlung.h"
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

double Bremsstrahlung::CalculatedEdx(double energy)
{
    if (multiplier_ <= 0)
    {
        return 0;
    }

    if (do_dedx_Interpolation_)
    {
        return max(dedx_interpolant_->Interpolate(energy), 0.0);
    }

    double sum = 0;

    for (int i = 0; i < (medium_->GetNumComponents()); i++)
    {
        IntegralLimits limits = SetIntegralLimits(energy, i);
        sum += dedx_integral_->Integrate(
            limits.vMin, limits.vUp, boost::bind(&Bremsstrahlung::FunctionToDEdxIntegral, this, energy, _1), 2);
    }

    return multiplier_ * energy * sum;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Bremsstrahlung::CalculatedNdx(double energy)
{
    if(multiplier_<=0)
    {
        return 0;
    }

    sum_of_rates_ = 0;

    for(int i=0; i<medium_->GetNumComponents(); i++)
    {
        if(do_dndx_Interpolation_)
        {
            prob_for_component_.at(i) = max(dndx_interpolant_1d_.at(i)->Interpolate(energy), 0.);
        }
        else
        {
            IntegralLimits limits = SetIntegralLimits(energy, i);
            prob_for_component_.at(i) = dndx_integral_.at(i)->Integrate(limits.vUp, limits.vMax, boost::bind(&Bremsstrahlung::FunctionToDNdxIntegral, this, energy,  _1),4);
        }
        sum_of_rates_ += prob_for_component_.at(i);
    }
    return sum_of_rates_;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Bremsstrahlung::CalculatedNdx(double energy, double rnd)
{
    if(multiplier_<=0)
    {
        return 0.;
    }

    // The random number will be stored to be able
    // to check if dNdx is already calculated for this random number.
    // This avoids a second calculation in CalculateStochaticLoss

    rnd_ = rnd;
    sum_of_rates_ = 0;

    for(int i=0; i<medium_->GetNumComponents(); i++)
    {
        if(do_dndx_Interpolation_)
        {
            prob_for_component_.at(i) = max(dndx_interpolant_1d_.at(i)->Interpolate(energy), 0.);
        }
        else
        {
            IntegralLimits limits = SetIntegralLimits(energy, i);
            prob_for_component_.at(i) = dndx_integral_.at(i)->IntegrateWithRandomRatio(limits.vUp, limits.vMax, boost::bind(&Bremsstrahlung::FunctionToDNdxIntegral, this, energy,  _1), 4, rnd);
        }
        sum_of_rates_ += prob_for_component_.at(i);
    }

    return sum_of_rates_;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Bremsstrahlung::CalculateStochasticLoss(double energy,  double rnd1, double rnd2)
{
    if(rnd1 != rnd_ )
    {
        CalculatedNdx(energy, rnd1);
    }

    return CalculateStochasticLoss(energy, rnd2);

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

// TODO: Compare to RadiationLength in Medium; just one (this or the Medium-thing) is needed
//formerlly in setlpm Calculates the radiation length of electrons
double Bremsstrahlung::CalculateScatteringX0(double energy)
{
    double X0;
    Integral* integral_temp = new Integral(IROMB,IMAXS,IPREC);

    bool store_lpm_effect_enabled = lpm_effect_enabled_;
    lpm_effect_enabled_ = false;

    double sum      =   0;

    bool store_init_lpm_effect_ = init_lpm_effect_;
    init_lpm_effect_    =   false;

    // PROPOSALParticle tmp_particle = PROPOSALParticle(particle_);

    for(int i=0; i < medium_->GetNumComponents(); i++)
    {

       IntegralLimits limits = SetIntegralLimits(energy, i);

       sum +=  integral_temp->Integrate(limits.vMin, limits.vUp, boost::bind(&Bremsstrahlung::FunctionToDEdxIntegral, this, energy,  _1),2);
       sum +=  integral_temp->Integrate(limits.vUp, limits.vMax, boost::bind(&Bremsstrahlung::FunctionToDEdxIntegral, this, energy,  _1),4);
    }

//    eLpm_        =   ALPHA*(particle->GetMass());
//    eLpm_        *=  eLpm_/(4*PI*ME*RE*sum);

    //TODO(mario): Check! Thu 2017/08/24
    //restore everything
    // SetIntegralLimits(0);
    lpm_effect_enabled_ = store_lpm_effect_enabled;
    init_lpm_effect_ = store_init_lpm_effect_;
    delete integral_temp;

    X0  =   pow(particle_def_.charge , 2);
    X0  =   X0*X0/sum;
    return X0;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-----------------------Enable and disable interpolation---------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Bremsstrahlung::EnableDNdxInterpolation( std::string path ,bool raw)
{
    if(do_dndx_Interpolation_)return;

    bool storing_failed =   false;
    bool reading_worked =   true;

    // charged anti leptons have the same cross sections like charged leptons
    // (except of diffractive Bremsstrahlung, where one can analyse the interference term if implemented)
    // so they use the same interpolation tables

    string particlename = particle_def_.name;

    if(!path.empty())
    {
        stringstream filename;
        filename<<path<<"/Brems_dNdx"
                <<"_particle"<<particlename
                <<"_mass_"<<particle_def_.mass
                <<"_charge_"<<particle_def_.charge
                <<"_lifetime_"<<particle_def_.lifetime
                <<"_para_"<<parametrization_
                <<"_med_"<<medium_->GetName()
                <<"_"<<medium_->GetMassDensity()
                <<"_ecut_"<<cut_settings_->GetEcut()
                <<"_vcut_"<<cut_settings_->GetVcut()
                <<"_lpm_"<<lpm_effect_enabled_
                <<"_multiplier_"<<multiplier_;

        if(!raw)
            filename<<".txt";

        dndx_interpolant_1d_.resize( medium_->GetNumComponents() );
        dndx_interpolant_2d_.resize( medium_->GetNumComponents() );

        if( FileExist(filename.str()) )
        {
            log_debug("Bremsstrahlungs parametrisation tables (dNdx) will be read from file:\t%s",filename.str().c_str());

            ifstream input;
            if(raw)
            {
                input.open(filename.str().c_str(), ios::binary);
            }
            else
            {
                input.open(filename.str().c_str());
            }
            for(int i=0; i<(medium_->GetNumComponents()); i++)
            {
                component_ = i;
                dndx_interpolant_2d_.at(i) = new Interpolant();
                dndx_interpolant_1d_.at(i) = new Interpolant();
                reading_worked = dndx_interpolant_2d_.at(i)->Load(input,raw);
                reading_worked = dndx_interpolant_1d_.at(i)->Load(input,raw);

            }
            input.close();
        }
        if(!FileExist(filename.str()) || !reading_worked )
        {
            if(!reading_worked)
            {
                log_info("file %s is corrupted! Write it again!",filename.str().c_str());
            }

            log_info("Info: Bremsstrahlungs parametrisation tables (dNdx) will be saved to file:\t%s",filename.str().c_str());

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

                for(int i=0; i<(medium_->GetNumComponents()); i++)
                {
                    component_ = i;

                    dndx_interpolant_2d_.at(i) =
                        new Interpolant(NUM1,
                                        particle_def_.low,
                                        BIGENERGY,
                                        NUM1,
                                        0,
                                        1,
                                        boost::bind(&Bremsstrahlung::FunctionToBuildDNdxInterpolant2D, this, _1, _2),
                                        order_of_interpolation_,
                                        false,
                                        false,
                                        true,
                                        order_of_interpolation_,
                                        false,
                                        false,
                                        false,
                                        order_of_interpolation_,
                                        true,
                                        false,
                                        false);
                    dndx_interpolant_1d_.at(i) =
                        new Interpolant(NUM1,
                                        particle_def_.low,
                                        BIGENERGY,
                                        boost::bind(&Bremsstrahlung::FunctionToBuildDNdxInterpolant, this, _1),
                                        order_of_interpolation_,
                                        false,
                                        false,
                                        true,
                                        order_of_interpolation_,
                                        true,
                                        false,
                                        false);

                    dndx_interpolant_2d_.at(i)->Save(output,raw);
                    dndx_interpolant_1d_.at(i)->Save(output,raw);

                }
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
        dndx_interpolant_1d_.resize( medium_->GetNumComponents() );
        dndx_interpolant_2d_.resize( medium_->GetNumComponents() );
        for(int i=0; i<(medium_->GetNumComponents()); i++)
        {
            component_ = i;
            dndx_interpolant_2d_.at(i) =
                new Interpolant(NUM1,
                                particle_def_.low,
                                BIGENERGY,
                                NUM1,
                                0,
                                1,
                                boost::bind(&Bremsstrahlung::FunctionToBuildDNdxInterpolant2D, this, _1, _2),
                                order_of_interpolation_,
                                false,
                                false,
                                true,
                                order_of_interpolation_,
                                false,
                                false,
                                false,
                                order_of_interpolation_,
                                true,
                                false,
                                false);
            dndx_interpolant_1d_.at(i) =
                new Interpolant(NUM1,
                                particle_def_.low,
                                BIGENERGY,
                                boost::bind(&Bremsstrahlung::FunctionToBuildDNdxInterpolant, this, _1),
                                order_of_interpolation_,
                                false,
                                false,
                                true,
                                order_of_interpolation_,
                                true,
                                false,
                                false);
        }
    }

    do_dndx_Interpolation_ = true;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

void Bremsstrahlung::EnableDEdxInterpolation( std::string path, bool raw)
{
    if(do_dedx_Interpolation_)return;

    bool reading_worked =   true;
    bool storing_failed =   false;

    // charged anti leptons have the same cross sections like charged leptons
    // (except of diffractive Bremsstrahlung, where one can analyse the interference term if implemented)
    // so they use the same interpolation tables
    string particlename = particle_def_.name;

    if(!path.empty())
    {
        stringstream filename;
        filename<<path<<"/Brems_dEdx"
                <<"_particle"<<particlename
                <<"_mass_"<<particle_def_.mass
                <<"_charge_"<<particle_def_.charge
                <<"_lifetime_"<<particle_def_.lifetime
                <<"_para_"<<parametrization_
                <<"_med_"<<medium_->GetName()
                <<"_"<<medium_->GetMassDensity()
                <<"_ecut_"<<cut_settings_->GetEcut()
                <<"_vcut_"<<cut_settings_->GetVcut()
                <<"_lpm_"<<lpm_effect_enabled_
                <<"_multiplier_"<<multiplier_;
        if(!raw)
            filename<<".txt";

        if( FileExist(filename.str()) )
        {
            log_debug("Bremsstrahlungs parametrisation tables (dEdx) will be read from file:\t%s",filename.str().c_str());
            ifstream input;

            if(raw)
            {
                input.open(filename.str().c_str(), ios::binary);
            }
            else
            {
                input.open(filename.str().c_str());
            }

            dedx_interpolant_ = new Interpolant();
            reading_worked = dedx_interpolant_->Load(input,raw);

            input.close();
        }
        if(!FileExist(filename.str()) || !reading_worked )
        {
            if(!reading_worked)
            {
                log_info("File %s is corrupted! Write it again!",filename.str().c_str());
            }

            log_info("Bremsstrahlungs parametrisation tables (dEdx) will be saved to file:%s\t",filename.str().c_str());

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

                dedx_interpolant_ = new Interpolant(NUM1, particle_def_.low, BIGENERGY, boost::bind(&Bremsstrahlung::FunctionToBuildDEdxInterpolant, this, _1),
                                                    order_of_interpolation_, true, false, true, order_of_interpolation_, false, false, true); //changed from ...,false,false,false)
                dedx_interpolant_->Save(output,raw);
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
        dedx_interpolant_ = new Interpolant(NUM1, particle_def_.low, BIGENERGY, boost::bind(&Bremsstrahlung::FunctionToBuildDEdxInterpolant, this, _1),
                                            order_of_interpolation_, true, false, true, order_of_interpolation_, false, false, true); //changed from ...,false,false,false)
    }

    do_dedx_Interpolation_=true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Bremsstrahlung::DisableDNdxInterpolation()
{

    for(unsigned int i = 0 ; i < dndx_interpolant_1d_.size() ; i++ )
    {
        delete dndx_interpolant_1d_.at(i);
    }

    for(unsigned int i = 0 ; i < dndx_interpolant_2d_.size() ; i++ )
    {
        delete dndx_interpolant_2d_.at(i);
    }

    dndx_interpolant_1d_.clear();
    dndx_interpolant_2d_.clear();

    do_dndx_Interpolation_  =   false;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Bremsstrahlung::DisableDEdxInterpolation()
{
    delete dedx_interpolant_;
    dedx_interpolant_   = NULL;
    do_dedx_Interpolation_  =   false;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Set and validate options--------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Bremsstrahlung::ValidateOptions()
{
    switch (parametrization_)
    {
        case ParametrizationType::BremsKelnerKokoulinPetrukhin:
        case ParametrizationType::BremsAndreevBezrukovBugaev:
        case ParametrizationType::BremsPetrukhinShestakov:
        case ParametrizationType::BremsCompleteScreeningCase:
            break;
        default:
            cerr<<"Bremsstrahlung: Parametrization type number is not vaild. \
                Set to default parametrization of KelnerKokoulinPetrukhin"<<endl;
            parametrization_ = ParametrizationType::BremsKelnerKokoulinPetrukhin;
            break;
    }
    if(order_of_interpolation_ < 2)
    {
        order_of_interpolation_ = 5;
        cerr<<"Bremsstrahlung: Order of Interpolation is not a vaild number\t"<<"Set to 5"<<endl;
    }
    if(order_of_interpolation_ > 6)
    {
        cerr<<"Bremsstrahlung: Order of Interpolation is set to "<<order_of_interpolation_
            <<".\t Note a order of interpolation > 6 will slow down the program"<<endl;
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Bremsstrahlung::Bremsstrahlung( )
    :CrossSections          ( )
    ,lorenz_                ( false )
    ,lorenz_cut_            ( 1e6 )
    ,component_             ( 0 )
    ,dndx_integral_         ( )
    ,dndx_interpolant_1d_   ( )
    ,dndx_interpolant_2d_   ( )
    ,eLpm_                  ( 0 )
    ,prob_for_component_    ( )
{
    dedx_integral_   =  new Integral(IROMB, IMAXS, IPREC);
    dedx_interpolant_= NULL;

    dndx_integral_.resize(medium_->GetNumComponents());

    for(int i =0; i<medium_->GetNumComponents(); i++)
    {
        dndx_integral_.at(i) = new Integral(IROMB, IMAXS, IPREC);
    }

    do_dedx_Interpolation_  = false;
    do_dndx_Interpolation_  = false;
    name_                   = "Bremsstrahlung";
    type_                   = ParticleType::Brems;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//Copyconstructor

Bremsstrahlung::Bremsstrahlung(const Bremsstrahlung &brems)
    :CrossSections                      ( brems )
    ,lorenz_                            ( brems.lorenz_ )
    ,lorenz_cut_                        ( brems.lorenz_cut_ )
    ,component_                         ( brems.component_ )
    ,dedx_integral_                     ( new Integral(*brems.dedx_integral_) )
    ,eLpm_                              ( brems.eLpm_)
    ,prob_for_component_                ( brems.prob_for_component_)
{
    if(brems.dedx_interpolant_ != NULL)
    {
        dedx_interpolant_ = new Interpolant(*brems.dedx_interpolant_) ;
    }
    else
    {
        dedx_interpolant_ = NULL;
    }


    dndx_integral_.resize( brems.dndx_integral_.size() );
    dndx_interpolant_1d_.resize( brems.dndx_interpolant_1d_.size() );
    dndx_interpolant_2d_.resize( brems.dndx_interpolant_2d_.size() );

    for(unsigned int i =0; i<brems.dndx_integral_.size(); i++)
    {
        dndx_integral_.at(i) = new Integral( *brems.dndx_integral_.at(i) );
    }
    for(unsigned int i =0; i<brems.dndx_interpolant_1d_.size(); i++)
    {
        dndx_interpolant_1d_.at(i) = new Interpolant( *brems.dndx_interpolant_1d_.at(i) );
    }
    for(unsigned int i =0; i<brems.dndx_interpolant_2d_.size(); i++)
    {
        dndx_interpolant_2d_.at(i) = new Interpolant( *brems.dndx_interpolant_2d_.at(i) );
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Bremsstrahlung::Bremsstrahlung(PROPOSALParticle& particle, Medium* medium,
                             EnergyCutSettings* cut_settings)
    :CrossSections          (particle, medium, cut_settings)
    ,lorenz_                ( false )
    ,lorenz_cut_            ( 1.e6 )
    ,dndx_integral_         ( )
    ,dndx_interpolant_1d_   ( )
    ,dndx_interpolant_2d_   ( )
    ,eLpm_                  ( 0 )
    ,prob_for_component_    ( )
{
    name_                       = "Bremsstrahlung";
    type_                       = ParticleType::Brems;
    ebig_                       = BIGENERGY;
    do_dedx_Interpolation_      = false;
    do_dndx_Interpolation_      = false;
    multiplier_                 = 1.;
    parametrization_            = ParametrizationType::BremsKelnerKokoulinPetrukhin;
    lpm_effect_enabled_         = false;
    init_lpm_effect_            = true;
    component_                  = 0;

    dedx_integral_   =  new Integral(IROMB, IMAXS, IPREC);
    dedx_interpolant_=  NULL;

    dndx_integral_.resize( medium_->GetNumComponents() );

    for(int i =0; i<(medium_->GetNumComponents()); i++)
    {
            dndx_integral_.at(i) =   new Integral(IROMB, IMAXS, IPREC);
    }

    prob_for_component_.resize(medium_->GetNumComponents());

    do_dedx_Interpolation_  = false;
    do_dndx_Interpolation_  = false;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Bremsstrahlung& Bremsstrahlung::operator=(const Bremsstrahlung &brems)
{

    if (this != &brems)
    {
      Bremsstrahlung tmp(brems);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



bool Bremsstrahlung::operator==(const Bremsstrahlung &brems) const
{
    if( this->CrossSections::operator !=(brems) )                           return false;

    if( lorenz_                     !=  brems.lorenz_)                      return false;
    if( lorenz_cut_                 !=  brems.lorenz_cut_)                  return false;
    if( component_                  !=  brems.component_)                   return false;
    if( eLpm_                       !=  brems.eLpm_)                        return false;
    if( *dedx_integral_             != *brems.dedx_integral_)               return false;
    if( prob_for_component_.size()  !=  brems.prob_for_component_.size())   return false;
    if( dndx_integral_.size()       !=  brems.dndx_integral_.size())        return false;
    if( dndx_interpolant_1d_.size() !=  brems.dndx_interpolant_1d_.size())  return false;
    if( dndx_interpolant_2d_.size() !=  brems.dndx_interpolant_2d_.size())  return false;

    for(unsigned int i =0; i<brems.dndx_integral_.size(); i++)
    {
        if( *dndx_integral_.at(i)       != *brems.dndx_integral_.at(i) )        return false;
    }
    for(unsigned int i =0; i<brems.dndx_interpolant_1d_.size(); i++)
    {
        if( *dndx_interpolant_1d_.at(i) != *brems.dndx_interpolant_1d_.at(i) )  return false;
    }
    for(unsigned int i =0; i<brems.dndx_interpolant_2d_.size(); i++)
    {
        if( *dndx_interpolant_2d_.at(i) != *brems.dndx_interpolant_2d_.at(i) )  return false;
    }
    for(unsigned int i =0; i<brems.prob_for_component_.size(); i++)
    {
        if( prob_for_component_.at(i) != brems.prob_for_component_.at(i) )      return false;
    }

    if( dedx_interpolant_ != NULL && brems.dedx_interpolant_ != NULL)
    {
        if( *dedx_interpolant_          != *brems.dedx_interpolant_)            return false;
    }
    else if( dedx_interpolant_ != brems.dedx_interpolant_)                      return false;

    //else
    return true;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Bremsstrahlung::operator!=(const Bremsstrahlung &bremsstrahlung) const
{
    return !(*this == bremsstrahlung);

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Bremsstrahlung::swap(Bremsstrahlung &brems)
{
    using std::swap;

    this->CrossSections::swap(brems);

    swap(lorenz_,brems.lorenz_);
    swap(lorenz_cut_,brems.lorenz_cut_);
    swap(component_,brems.component_);
    dedx_integral_->swap(*brems.dedx_integral_);

    if( dedx_interpolant_ != NULL && brems.dedx_interpolant_ != NULL)
    {
        dedx_interpolant_->swap(*brems.dedx_interpolant_);
    }
    else if( dedx_interpolant_ == NULL && brems.dedx_interpolant_ != NULL)
    {
        dedx_interpolant_ = new Interpolant(*brems.dedx_interpolant_);
        brems.dedx_interpolant_ = NULL;
    }
    else if( dedx_interpolant_ != NULL && brems.dedx_interpolant_ == NULL)
    {
        brems.dedx_interpolant_ = new Interpolant(*dedx_interpolant_);
        dedx_interpolant_ = NULL;
    }

    swap(eLpm_,brems.eLpm_);

    prob_for_component_.swap(brems.prob_for_component_);
    dndx_integral_.swap(brems.dndx_integral_);
    dndx_interpolant_1d_.swap(brems.dndx_interpolant_1d_);
    dndx_interpolant_2d_.swap(brems.dndx_interpolant_2d_);

}

namespace PROPOSAL
{

ostream& operator<<(std::ostream& os, Bremsstrahlung const &brems)
{
    os<<"---------------------------Bremsstrahlung( "<<&brems<<" )---------------------------"<<std::endl;
    os<< static_cast <const CrossSections &>( brems ) << endl;
    os<< "------- Class Specific: " << endl;
    os<< "\tlorenz:\t\t" << brems.lorenz_ << endl;
    os<< "\tlorenz_cut:\t\t" << brems.lorenz_cut_ << endl;
    os<<endl;
    os<<"\tdedx_integral:\t"<<brems.dedx_integral_ << endl;
    os<<"\tdndx_integral:\t"<<brems.dndx_integral_.size()<<endl;
    for(unsigned int i=0;i<brems.dndx_integral_.size();i++)
    {
        os<<"\t\tadress:\t\t"<<brems.dndx_integral_.at(i)<<endl;
    }
    os<<endl;
    os<<"\tdo_dedx_Interpolation:\t\t"<<brems.do_dedx_Interpolation_<<endl;
    os<<"\tdedx_interpolant:\t\t"<<brems.dedx_interpolant_<<endl;
    os<<endl;
    os<<"\tdo_dndx_Interpolation:\t\t"<<brems.do_dndx_Interpolation_<<endl;
    os<<"\tdndx_interpolant_1d:\t\t"<<brems.dndx_interpolant_1d_.size()<<endl;
    for(unsigned int i=0;i<brems.dndx_interpolant_1d_.size();i++)
    {
        os<<"\t\tadress:\t\t"<<brems.dndx_interpolant_1d_.at(i)<<endl;
    }
    os<<"\tdndx_interpolant_2d:\t\t"<<brems.dndx_interpolant_2d_.size()<<endl;
    for(unsigned int i=0;i<brems.dndx_interpolant_2d_.size();i++)
    {
        os<<"\t\tadress:\t\t"<<brems.dndx_interpolant_2d_.at(i)<<endl;
    }
    os<<endl;
    os<<"\tTemp. variables: " << endl;
    os<< "\t\tcomponent:\t\t" << brems.component_<<endl;
    os<< "\t\teLpm:\t\t" << brems.eLpm_<<endl;
    os<<"\t\tprob_for_component:\t"<<brems.prob_for_component_.size()<<endl;
    for(unsigned int i=0;i<brems.prob_for_component_.size();i++)
    {
        os<<"\t\t\tvalue:\t\t"<<brems.prob_for_component_.at(i)<<endl;
    }
    os<<"-----------------------------------------------------------------------------------------------";
    return os;
}

}  // namespace PROPOSAL

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------Parametrizations--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Bremsstrahlung::KelnerKokoulinPetrukhinParametrization(double energy, double v, int i)
{
    double Z3       =   0;
    double result   =   0;
    double Dn       =   0;
    double s1       =   0;

    Components::Component* component = medium_->GetComponents().at(i);

    Z3  =   pow(component->GetNucCharge(), -1./3);

    int step;
    double d, da, dn, Fa, maxV;

    //TODO(mario): Better way? Sat 2017/09/02
    double square_momentum = energy * energy - particle_def_.mass * particle_def_.mass;
    double particle_momentum = sqrt(max(square_momentum, 0.0));

    d       =   particle_def_.mass*particle_def_.mass
                *v/(2*(energy)*(1-v));
    s1      =   (component->GetLogConstant())*Z3;
    da      =   log(1 + ME/(d*SQRTE*s1));
    Dn      =   1.54*pow((component->GetAtomicNum()), 0.27);
    s1      =   ME*Dn/((particle_def_.mass)*s1);
    dn      =   log(Dn/(1 + d*(Dn*SQRTE - 2)/particle_def_.mass));
    maxV    =   ME*(energy - particle_def_.mass)
                /((energy)
                  *(energy - particle_momentum + ME));

    if(v<maxV)
    {
        Fa  =   log(((particle_def_.mass)/d)/(d*(particle_def_.mass)/(ME*ME) + SQRTE)) -
                log(1 + ME/(d*SQRTE*component->GetBPrime()*(pow(component->GetNucCharge() , -2./3))));
    }
    else
    {
        Fa  =   0;
    }

    if((component->GetNucCharge())==1)
    {
        step    =   0;
    }

    else
    {
        step    =   1;
    }


    result = ((4./3)*(1-v) + v*v)
            *(log((particle_def_.mass)/d)
              - 0.5 -da - dn + (dn*step + Fa)/(component->GetNucCharge()));

    return result;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Bremsstrahlung::AndreevBezrukovBugaevParametrization(double energy, double v, int i)
{

    Components::Component* component = medium_->GetComponents().at(i);

    double aux      =   0;
    double Z3       =   0;
    double result   =   0;

    Z3 = pow((component->GetNucCharge()), -1./3);

    double aux1, aux2, a1, a2,zeta, qc, qmin, x1, x2, d1,d2, psi1, psi2;

    a1      =   111.7*Z3/ME;
    a2      =   724.2*Z3*Z3/ME;
    qc      =   1.9*MMU*Z3;
    aux     =   2*(particle_def_.mass)/qc;
    aux     *=  aux;
    zeta    =   sqrt(1+aux);
    qmin    =   pow((particle_def_.mass),2)
                *v/((energy)*(1-v));

    x1      =   a1*qmin;
    x2      =   a2*qmin;

    if((component->GetNucCharge())==1)
    {
        d1  =   0;
        d2  =   0;
    }
    else
    {
        aux1    =   log((particle_def_.mass)/qc);
        aux2    =   (zeta/2)*log((zeta+1)/(zeta-1));
        d1      =   aux1 + aux2;
        d2      =   aux1 + ((3 - pow(zeta , 2))*aux2 + aux)/2;
    }

    aux     =   (particle_def_.mass)*a1;
    aux1    =   log(pow(aux , 2)/(1 + pow(x1 , 2)));
    aux     =   (particle_def_.mass)*a2;
    aux2    =   log(pow(aux , 2)/(1 + pow(x2 , 2)));
    psi1    =   (1+ aux1)/2 + (1 + aux2)/(2*(component->GetNucCharge()));
    psi2    =   (2./3 + aux1)/2 +
                (2./3 + aux2)/(2*(component->GetNucCharge()));

    aux1    =   x1*atan(1/x1);
    aux2    =   x2*atan(1/x2);
    psi1    -=  aux1 + aux2/(component->GetNucCharge());
    aux     =   pow(x1 , 2);
    psi2    +=  2*aux*(1 - aux1 + 3./4*log(aux/(1 + aux)));
    aux     =   pow(x2 , 2);
    psi2    +=  2*aux*(1 - aux2 + 3./4*log(aux/(1 + aux)))
                /(component->GetNucCharge());

    psi1    -=  d1;
    psi2    -=  d2;
    result  =   (2-2*v + pow(v , 2))*psi1 - (2./3)*(1-v)*psi2;

    if(result<0)
    {
        result  =   0;
    }

    return result;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Bremsstrahlung::PetrukhinShestakovParametrization(double energy, double v, int i)
{
    Components::Component* component = medium_->GetComponents().at(i);

    double Z3       =   0;
    double result   =   0;
    double d, Fd;

    Z3  =   pow((component->GetNucCharge()), -1./3);

    d   =   pow((particle_def_.mass) , 2)
            * v/(2*(energy)*(1-v));

    Fd  =   189*Z3/ME;
    Fd  =   (particle_def_.mass)*Fd/(1 + SQRTE*d*Fd);

    if((component->GetNucCharge())>10)
    {
        Fd  *=  (2./3)*Z3;
    }

    result  =   ((4./3)*(1-v) + pow(v , 2))*log(Fd);

    return result;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Bremsstrahlung::CompleteScreeningCase(double energy, double v, int i)
{
    (void) energy;
    Components::Component* component = medium_->GetComponents().at(i);

    double aux      =   0;
    double Z3       =   0;
    double result   =   0;
    double Lr, fZ, Lp;

    Z3  =   pow((component->GetNucCharge()) , -1./3);

    aux =   ALPHA*(component->GetNucCharge());
    aux *=  aux;
    fZ  =   aux*(1/(1 + aux) + 0.20206 + aux*(-0.0369 + aux*(0.0083 - 0.002*aux)));

    //check rounding
    switch((int)((component->GetNucCharge()) + 0.5))
    {

        case 1:
        {
            Lr  =   5.31;
            Lp  =   6.144;
        }break;

        case 2:
        {
            Lr  =   4.79;
            Lp  =   5.621;
        }break;

        case 3:
        {
            Lr  =   4.74;
            Lp  =   5.805;
        }break;

        case 4:
        {
            Lr  =   4.71;
            Lp  =   5.924;
        }break;

        default:
        {
            Lr  =   log(184.15*Z3);
            Lp  =   log (1194*pow(Z3 , 2));
        }break;

    }

    result = (((4./3)*(1-v) + pow(v , 2))*
              ((component->GetNucCharge())*(Lr - fZ) + Lp)
             + (1./9)*(1-v)*((component->GetNucCharge()) + 1))
            /(component->GetNucCharge());

    return result;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------Cross section /lpm / limit----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Bremsstrahlung::CalculateStochasticLoss(double energy, double rnd)
{

    double rand;
    double rsum;

    rand    =   rnd*sum_of_rates_;
    rsum    =   0;

    for(int i=0; i<(medium_->GetNumComponents()); i++)
    {
        rsum    += prob_for_component_.at(i);

        if(rsum > rand)
        {
            if(do_dndx_Interpolation_)
            {
                IntegralLimits limits = SetIntegralLimits(energy, i);

                if(limits.vUp == limits.vMax)
                {
                    return (energy)*limits.vUp;
                }

                return (energy)*(limits.vUp*exp(dndx_interpolant_2d_.at(i)->FindLimit(energy, (rnd)*prob_for_component_.at(i))*log(limits.vMax/limits.vUp)));
            }

            else
            {
                component_ = i;

                return (energy)*dndx_integral_.at(i)->GetUpperLimit();

            }
        }
    }

    //sometime everything is fine, just the probability for interaction is zero
    bool prob_for_all_comp_is_zero=true;
    for(int i=0; i<(medium_->GetNumComponents()); i++)
    {
        IntegralLimits limits = SetIntegralLimits(energy, i);
        if(limits.vUp!= limits.vMax) prob_for_all_comp_is_zero=false;
    }
    if(prob_for_all_comp_is_zero)return 0;

    log_fatal("sum was not initialized correctly");
    return 0;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Bremsstrahlung::ElasticBremsstrahlungCrossSection(double energy, double v, int i)
{
    Components::Component* component = medium_->GetComponents().at(i);

    double aux      =   0;
    double Z3       =   0;
    double result   =   0;
    double Dn       =   0;
    double s1       =   0;

    Z3  =   pow((component->GetNucCharge()), -1./3);

    switch(parametrization_)
    {
        case ParametrizationType::BremsKelnerKokoulinPetrukhin:
            result  =   KelnerKokoulinPetrukhinParametrization(energy, v, i);
            break;
        case ParametrizationType::BremsAndreevBezrukovBugaev:
            result  =   AndreevBezrukovBugaevParametrization(energy, v, i);
            break;
        case ParametrizationType::BremsPetrukhinShestakov:
            result  =   PetrukhinShestakovParametrization(energy, v, i);
            break;
        case ParametrizationType::BremsCompleteScreeningCase:
            result  =   CompleteScreeningCase(energy, v, i);
            break;
        default:
            log_info("Bremsstrahlung: Parametrization type number '%i' is not vaild. \
                Set to default parametrization of KelnerKokoulinPetrukhin with type number '%i' "
                , parametrization_, ParametrizationType::BremsKelnerKokoulinPetrukhin);
    }

    aux =   2*(component->GetNucCharge())*(ME/particle_def_.mass)*RE;
    aux *=  (ALPHA/v)*aux*result;

    if(lpm_effect_enabled_)
    {
        if(parametrization_!=ParametrizationType::BremsKelnerKokoulinPetrukhin)
        {
            s1  =   (component->GetLogConstant())*Z3;
            Dn  =   1.54*pow((component->GetAtomicNum()) , 0.27);
            s1  =   ME*Dn/((particle_def_.mass)*s1);
        }
        aux *=  lpm(energy, v,s1);
    }

    double c2   =   pow(particle_def_.charge , 2);

    return medium_->GetMolDensity()*component->GetAtomInMolecule()*c2*c2*aux;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Bremsstrahlung::lpm(double energy, double v, double s1)
{
    if(init_lpm_effect_)
    {
        Integral* integral_temp = new Integral(IROMB,IMAXS,IPREC);

        lpm_effect_enabled_ = false;
        double sum      =   0;
        init_lpm_effect_    =   false;

        // PROPOSALParticle tmp_particle = PROPOSALParticle(particle_);
        // tmp_particle.SetEnergy(BIGENERGY);

        for(int i=0; i < medium_->GetNumComponents(); i++)
        {

           IntegralLimits limits = SetIntegralLimits(BIGENERGY, i);

           sum +=  integral_temp->Integrate(limits.vMin, limits.vUp, boost::bind(&Bremsstrahlung::FunctionToDEdxIntegral, this, energy,  _1),2);
           sum +=  integral_temp->Integrate(limits.vUp, limits.vMax, boost::bind(&Bremsstrahlung::FunctionToDEdxIntegral, this, energy,  _1),4);
        }

        eLpm_        =   ALPHA*(particle_def_.mass);
        eLpm_        *=  eLpm_/(4*PI*ME*RE*sum);


        //TODO(mario): Ceck Thu 2017/08/24
        // SetIntegralLimits(0);
        lpm_effect_enabled_=true;
        delete integral_temp;

    }

    double G, fi, xi, sp, h, s, s2, s3, ps, Gamma;

    const double fi1    =   1.54954;
    const double G1     =   0.710390;
    const double G2     =   0.904912;
    s1                  *=  s1*SQRT2;
    sp                  =   sqrt(eLpm_*v/(8*(energy)*(1-v)));
    h                   =   log(sp)/log(s1);

    if(sp < s1)
    {
        xi  =   2;
    }
    else if(sp < 1)
    {
        xi  =   1 + h - 0.08*(1 - h)*(1 - (1-h)*(1 - h))/log(s1);
    }
    else
    {
        xi  =   1;
    }

    s       =   sp/sqrt(xi);
    Gamma   =   RE*ME/(ALPHA*(particle_def_.mass)*v);
    Gamma   =   1 +4*PI*(medium_->GetMolDensity())*(medium_->GetSumCharge())*RE*pow(Gamma,2);
    s       *=  Gamma;
    s2      =   pow(s,2);
    s3      =   pow(s,3);

    if(s < fi1)
    {
        fi  =   1-exp(-6*s*(1 + (3-PI)*s) + s3/(0.623 + 0.796*s + 0.658*s2));
    }
    else
    {
        fi  =   1 - 0.012/pow(s2 , 2);
    }

    if(s < G1)
    {
        ps  =   1 - exp(-4*s - 8*s2/(1 + 3.936*s + 4.97*s2 - 0.05*s3 + 7.50*pow(s2 , 2)));
        G   =   3*ps - 2*fi;
    }
    else if (s < G2)
    {
        G   =   36*s2/(36*s2 + 1);
    }
    else
    {
        G   =   1 - 0.022/pow(s2 , 2);
    }

    return ((xi/3)*((v*v)*G/(Gamma*Gamma) + 2*(1 + (1-v)*(1-v))*fi/Gamma))/((4./3)*(1-v) + v*v);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


CrossSections::IntegralLimits Bremsstrahlung::SetIntegralLimits(double energy, int interaction_component)
{
    component_ = interaction_component;
    Components::Component* component = medium_->GetComponents().at(component_);

    IntegralLimits limits;

    limits.vMin = 0.;
    limits.vMax = 1 - (3./4)*SQRTE*(particle_def_.mass/energy)
                *pow((component->GetNucCharge()) , 1./3);

    if(limits.vMax<0)
    {
        limits.vMax   =   0;
    }

    if(lorenz_)
    {
        limits.vMax   =   min(limits.vMax, lorenz_cut_/(energy));
    }

    limits.vMax = min(limits.vMax, (1 - (particle_def_.mass / energy)));
    limits.vUp  = min(limits.vMax, cut_settings_->GetCut(energy));

    return limits;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------Functions to interpolate---------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Bremsstrahlung::FunctionToBuildDNdxInterpolant(double energy)
{
    return dndx_interpolant_2d_.at(component_)->Interpolate(energy, 1.);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Bremsstrahlung::FunctionToBuildDNdxInterpolant2D( double energy , double v)
{
    // double tmp_energy = particle_.GetEnergy();
    // particle_.SetEnergy(energy);

    IntegralLimits limits = SetIntegralLimits(energy, component_);

    if(limits.vUp== limits.vMax)
    {
        return 0;
    }

    v = limits.vUp * exp(v * log(limits.vMax / limits.vUp));

    double dNdx =  dndx_integral_.at(component_)->Integrate(limits.vUp, v, boost::bind(&Bremsstrahlung::FunctionToDNdxIntegral, this, energy, _1),4);
    // particle_.SetEnergy(tmp_energy);

    return dNdx;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Bremsstrahlung::FunctionToBuildDEdxInterpolant( double energy)
{
    // double tmp_energy = particle_.GetEnergy();
    // particle_.SetEnergy(energy);

    // double dEdx =  CalculatedEdx();
    //
    // particle_.SetEnergy(tmp_energy);

    // return dEdx;
    return CalculatedEdx(energy);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Functions to integrate----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Bremsstrahlung::FunctionToDEdxIntegral(double energy, double variable)
{
    return variable * ElasticBremsstrahlungCrossSection(energy, variable, component_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Bremsstrahlung::FunctionToDNdxIntegral(double energy, double variable)
{
    return multiplier_ * ElasticBremsstrahlungCrossSection(energy, variable, component_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// void Bremsstrahlung::SetParametrization(ParametrizationType::Enum parametrization)
// {
//     parametrization_ = parametrization;
//     switch (parametrization_)
//     {
//         case ParametrizationType::BremsKelnerKokoulinPetrukhin:
//         case ParametrizationType::BremsAndreevBezrukovBugaev:
//         case ParametrizationType::BremsPetrukhinShestakov:
//         case ParametrizationType::BremsCompleteScreeningCase:
//             break;
//         default:
//             log_warn("Bremsstrahlung: Parametrization type number '%i' is not vaild. \
//                 Set to default parametrization of KelnerKokoulinPetrukhin with type number '%i' "
//                 , parametrization, ParametrizationType::BremsKelnerKokoulinPetrukhin);
//             parametrization_ = ParametrizationType::BremsKelnerKokoulinPetrukhin;
//             break;
//     }
//
//     if(do_dedx_Interpolation_)
//     {
//         log_warn("dEdx-interpolation enabled before choosing the parametrization. Building the tables again");
//         DisableDEdxInterpolation();
//         EnableDEdxInterpolation();
//     }
//     if(do_dndx_Interpolation_)
//     {
//         log_warn("dNdx-interpolation enabled before choosing the parametrization. Building the tables again");
//         DisableDNdxInterpolation();
//         EnableDNdxInterpolation();
//     }
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Destructor---------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Bremsstrahlung::~Bremsstrahlung()
{
    DisableDNdxInterpolation();
    DisableDEdxInterpolation();

    delete dedx_integral_;
    for(unsigned int i = 0 ; i < dndx_integral_.size() ; i++ ){
        delete dndx_integral_[i];
    }

    dndx_integral_.clear();
}
