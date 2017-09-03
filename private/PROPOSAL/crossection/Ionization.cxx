
// #include <algorithm>

#include <boost/bind.hpp>

#include "PROPOSAL/crossection/Ionization.h"
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


double Ionization::CalculatedEdx(double energy)
{
    if(multiplier_<=0)
    {
        return 0;
    }

    if(do_dedx_Interpolation_)
    {
        return max(dedx_interpolant_->Interpolate(energy), 0.);
    }

    double result, aux;

    CrossSections::IntegralLimits limits = SetIntegralLimits(energy, 0);;

    //TODO(mario): Better way? Sat 2017/09/02
    double square_momentum = energy * energy - particle_def_.mass * particle_def_.mass;
    double particle_momentum = sqrt(max(square_momentum, 0.0));
    double beta    =   particle_momentum/energy;
    double gamma   =   energy/particle_def_.mass;

    aux     =   beta*gamma/(1.e-6*medium_->GetI());
    result  =   log(limits.vUp*(2*ME*energy))+2*log(aux);
    aux     =   limits.vUp/(2*(1 + 1/gamma));
    result  +=  aux*aux;
    aux     =   beta*beta;
    result  -=  aux*(1 + limits.vUp/limits.vMax) + Delta(beta, gamma);

    if(result>0)
    {
        result*=IONK*particle_def_.charge*particle_def_.charge*medium_->GetZA()/(2*aux);
    }
    else
    {
        result=0;
    }
    return multiplier_*(medium_->GetMassDensity()*result
                        + energy*(integral_->Integrate(limits.vMin, limits.vUp, boost::bind(&Ionization::FunctionToDEdxIntegral, this, energy, _1),4)));
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Ionization::CalculatedNdx(double energy)
{
    if(multiplier_<=0)
    {
        return 0;
    }

    if(do_dndx_Interpolation_)
    {
        sum_of_rates_ = max(dndx_interpolant_1d_->Interpolate(energy), 0.);
    }
    else{
        CrossSections::IntegralLimits limits = SetIntegralLimits(energy, 0);;
        sum_of_rates_ = integral_->Integrate(limits.vUp,limits.vMax,boost::bind(&Ionization::FunctionToDNdxIntegral, this, energy, _1),3,1);
    }
    return sum_of_rates_;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Ionization::CalculatedNdx(double energy, double rnd)
{
    if(multiplier_<=0)
    {
        return 0;
    }

    // The random number will be stored to be able
    // to check if dNdx is already calculated for this random number.
    // This avoids a second calculation in CalculateStochasticLoss

    rnd_    =   rnd;

    if(do_dndx_Interpolation_)
    {
        sum_of_rates_ = max(dndx_interpolant_1d_->Interpolate(energy), 0.);
    }
    else
    {
        CrossSections::IntegralLimits limits = SetIntegralLimits(energy, 0);;
        sum_of_rates_ = integral_->IntegrateWithRandomRatio(limits.vUp,limits.vMax,boost::bind(&Ionization::FunctionToDNdxIntegral, this, energy, _1),3,rnd,1);
    }

    return sum_of_rates_;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Ionization::CalculateStochasticLoss(double energy, double rnd1, double rnd2)
{
    if(rnd1 != rnd_ )
    {
        CalculatedNdx( rnd1);
    }

    return CalculateStochasticLoss(energy, rnd2);
}



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-----------------------Enable and disable interpolation---------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Ionization::EnableDNdxInterpolation( std::string path, bool raw)
{

    if(do_dndx_Interpolation_)return;

    bool reading_worked =   true;
    bool storing_failed =   false;

    // charged anti leptons have the same cross sections like charged leptons
    // so they use the same interpolation tables
    string particle_name = particle_def_.name;

    if(!path.empty())
    {
        stringstream filename;
        filename<<path<<"/Ioniz_dNdx"
                <<"_particle_"<<particle_name
                <<"_mass_"<<particle_def_.mass
                <<"_charge_"<<particle_def_.charge
                <<"_lifetime_"<<particle_def_.lifetime
                <<"_med_"<<medium_->GetName()
                <<"_"<<medium_->GetMassDensity()
                <<"_ecut_"<<cut_settings_->GetEcut()
                <<"_vcut_"<<cut_settings_->GetVcut()
                <<"_multiplier_"<<multiplier_;

        if(!raw)
            filename<<".txt";

        if( FileExist(filename.str()) )
        {
            log_debug("Ionization parametrisation tables (dNdx) will be read from file:\t%s",filename.str().c_str());
            ifstream input;

            if(raw)
            {
                input.open(filename.str().c_str(), ios::binary);
            }
            else
            {
                input.open(filename.str().c_str());
            }

            dndx_interpolant_2d_ = new Interpolant();
            dndx_interpolant_1d_ = new Interpolant();
            reading_worked = dndx_interpolant_2d_->Load(input, raw);
            reading_worked = dndx_interpolant_1d_->Load(input, raw);

            input.close();
        }
        if(!FileExist(filename.str()) || !reading_worked )
        {
            if(!reading_worked)
            {
                log_info("File %s is corrupted! Write it again!",filename.str().c_str());
            }

            log_info("Ionization parametrisation tables (dNdx) will be saved to file:\t%s",filename.str().c_str());

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

                dndx_interpolant_2d_ = new Interpolant(NUM1
                                                    , particle_def_.low
                                                    , BIGENERGY
                                                    , NUM1
                                                    , 0
                                                    , 1
                                                    , boost::bind(&Ionization::FunctionToBuildDNdxInterpolant2D, this,  _1, _2)
                                                    , order_of_interpolation_
                                                    , false
                                                    , false
                                                    , true
                                                    , order_of_interpolation_
                                                    , false
                                                    , false
                                                    , false
                                                    , order_of_interpolation_
                                                    , true
                                                    , false
                                                    , false
                                                    );
                dndx_interpolant_1d_ = new Interpolant(NUM1
                                                    , particle_def_.low
                                                    , BIGENERGY
                                                    , boost::bind(&Ionization::FunctionToBuildDNdxInterpolant, this, _1)
                                                    , order_of_interpolation_
                                                    , false
                                                    , false
                                                    , true
                                                    , order_of_interpolation_
                                                    , true
                                                    , false
                                                    , false
                                                    );

                dndx_interpolant_2d_->Save(output, raw);
                dndx_interpolant_1d_->Save(output, raw);

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

        dndx_interpolant_2d_ = new Interpolant(NUM1
                                            , particle_def_.low
                                            , BIGENERGY
                                            , NUM1
                                            , 0
                                            , 1
                                            , boost::bind(&Ionization::FunctionToBuildDNdxInterpolant2D, this,  _1, _2)
                                            , order_of_interpolation_
                                            , false
                                            , false
                                            , true
                                            , order_of_interpolation_
                                            , false
                                            , false
                                            , false
                                            , order_of_interpolation_
                                            , true
                                            , false
                                            , false
                                            );
        dndx_interpolant_1d_ = new Interpolant(NUM1
                                            , particle_def_.low
                                            , BIGENERGY
                                            , boost::bind(&Ionization::FunctionToBuildDNdxInterpolant, this, _1)
                                            , order_of_interpolation_
                                            , false
                                            , false
                                            , true
                                            , order_of_interpolation_
                                            , true
                                            , false
                                            , false
                                            );
    }

    do_dndx_Interpolation_=true;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Ionization::EnableDEdxInterpolation( std::string path, bool raw)
{

    if(do_dedx_Interpolation_)return;

    bool reading_worked =   true;
    bool storing_failed =   false;

    // charged anti leptons have the same cross sections like charged leptons
    // so they use the same interpolation tables
    string particle_name = particle_def_.name;

    if(!path.empty())
    {
        stringstream filename;
        filename<<path<<"/Ioniz_dEdx"
                <<"_particle_"<<particle_name
                <<"_mass_"<<particle_def_.mass
                <<"_charge_"<<particle_def_.charge
                <<"_lifetime_"<<particle_def_.lifetime
                <<"_med_"<<medium_->GetName()
                <<"_"<<medium_->GetMassDensity()
                <<"_ecut_"<<cut_settings_->GetEcut()
                <<"_vcut_"<<cut_settings_->GetVcut()
                <<"_multiplier_"<<multiplier_;

        if(!raw)
            filename<<".txt";

        if( FileExist(filename.str()) )
        {
            log_debug("Ionization parametrisation tables (dEdx) will be read from file:\t%s",filename.str().c_str());
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
            reading_worked = dedx_interpolant_->Load(input, raw);

            input.close();
        }
        if(!FileExist(filename.str()) || !reading_worked )
        {
            if(!reading_worked)
            {
                log_info("File %s is corrupted! Write it again!",filename.str().c_str());
            }

            log_info("Ionization parametrisation tables (dEdx) will be saved to file:\t%s",filename.str().c_str());

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

                dedx_interpolant_ = new Interpolant(NUM1
                                                , particle_def_.low
                                                , BIGENERGY
                                                , boost::bind(&Ionization::FunctionToBuildDEdxInterpolant, this,  _1)
                                                , order_of_interpolation_
                                                , true
                                                , false
                                                , true
                                                , order_of_interpolation_
                                                , false
                                                , false
                                                , true
                                                );
                dedx_interpolant_->Save(output, raw);
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
        dedx_interpolant_ = new Interpolant(NUM1
                                        , particle_def_.low
                                        , BIGENERGY
                                        , boost::bind(&Ionization::FunctionToBuildDEdxInterpolant, this,  _1)
                                        , order_of_interpolation_
                                        , true
                                        , false
                                        , true
                                        , order_of_interpolation_
                                        , false
                                        , false
                                        , true
                                        );
    }

    do_dedx_Interpolation_=true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Ionization::DisableDNdxInterpolation()
{
    delete dndx_interpolant_1d_;
    delete dndx_interpolant_2d_;

    dndx_interpolant_1d_    =   NULL;
    dndx_interpolant_2d_    =   NULL;
    do_dndx_Interpolation_  =   false;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Ionization::DisableDEdxInterpolation()
{
    delete dedx_interpolant_;

    dedx_interpolant_       =   NULL;
    do_dedx_Interpolation_  =   false;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Set and validate options--------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Ionization::ValidateOptions()
{
    if(order_of_interpolation_ < 2)
    {
        order_of_interpolation_ = 5;
        cerr<<"Ionization: Order of Interpolation is not a vaild number\t"<<"Set to 5"<<endl;
    }
    if(order_of_interpolation_ > 6)
    {
        cerr<<"Ionization: Order of Interpolation is set to "<<order_of_interpolation_
            <<".\t Note a order of interpolation > 6 will slow down the program"<<endl;
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Ionization::Ionization()
    :integral_              (  )
{
    name_                  = "Ionization";
    type_                  = ParticleType::DeltaE;

    dedx_interpolant_      = NULL;
    dndx_interpolant_1d_   = NULL;
    dndx_interpolant_2d_   = NULL;
    integral_              = new Integral();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Ionization::Ionization(const Ionization &ioniz)
    :CrossSections          ( ioniz  )
    ,integral_              ( new Integral(*ioniz.integral_) )
{
    if(ioniz.dedx_interpolant_ != NULL)
    {
        dedx_interpolant_ = new Interpolant(*ioniz.dedx_interpolant_) ;
    }
    else
    {
        dedx_interpolant_ = NULL;
    }

    if(ioniz.dndx_interpolant_1d_ != NULL)
    {
        dndx_interpolant_1d_ = new Interpolant(*ioniz.dndx_interpolant_1d_) ;
    }
    else
    {
        dndx_interpolant_1d_ = NULL;
    }

    if(ioniz.dndx_interpolant_2d_ != NULL)
    {
        dndx_interpolant_2d_ = new Interpolant(*ioniz.dndx_interpolant_2d_) ;
    }
    else
    {
        dndx_interpolant_2d_ = NULL;
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Ionization::Ionization(PROPOSALParticle& particle, Medium* medium, EnergyCutSettings* cut_settings)
    : CrossSections(particle, medium, cut_settings)
{
    name_                       = "Ionization";
    type_                       = ParticleType::DeltaE;
    ebig_                       = BIGENERGY;
    do_dedx_Interpolation_      = false;
    do_dndx_Interpolation_      = false;
    multiplier_                 = 1.;


    integral_              = new Integral(IROMB, IMAXS, IPREC);
    dedx_interpolant_      = NULL;
    dndx_interpolant_1d_   = NULL;
    dndx_interpolant_2d_   = NULL;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Ionization& Ionization::operator=(const Ionization &ioniz)
{
    if (this != &ioniz)
    {
      Ionization tmp(ioniz);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Ionization::operator==(const Ionization &ioniz) const
{
    if( this->CrossSections::operator !=(ioniz) )               return false;

    if( *integral_              != *ioniz.integral_)            return false;

    if( dedx_interpolant_ != NULL && ioniz.dedx_interpolant_ != NULL)
    {
        if( *dedx_interpolant_          != *ioniz.dedx_interpolant_)        return false;
    }
    else if( dedx_interpolant_ != ioniz.dedx_interpolant_)                  return false;

    if( dndx_interpolant_1d_ != NULL && ioniz.dndx_interpolant_1d_ != NULL)
    {
        if( *dndx_interpolant_1d_   != *ioniz.dndx_interpolant_1d_)         return false;
    }
    else if( dndx_interpolant_1d_ != ioniz.dndx_interpolant_1d_)            return false;

    if( dndx_interpolant_2d_ != NULL && ioniz.dndx_interpolant_2d_ != NULL)
    {
        if( *dndx_interpolant_2d_   != *ioniz.dndx_interpolant_2d_)         return false;
    }
    else if( dndx_interpolant_2d_ != ioniz.dndx_interpolant_2d_)            return false;


    //else
    return true;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Ionization::operator!=(const Ionization &ioniz) const
{
    return !(*this == ioniz);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Ionization::swap(Ionization &ioniz)
{
    using std::swap;

    this->CrossSections::swap(ioniz);

    integral_->swap( *ioniz.integral_);

    if( dedx_interpolant_ != NULL && ioniz.dedx_interpolant_ != NULL)
    {
        dedx_interpolant_->swap(*ioniz.dedx_interpolant_);
    }
    else if( dedx_interpolant_ == NULL && ioniz.dedx_interpolant_ != NULL)
    {
        dedx_interpolant_ = new Interpolant(*ioniz.dedx_interpolant_);
        ioniz.dedx_interpolant_ = NULL;
    }
    else if( dedx_interpolant_ != NULL && ioniz.dedx_interpolant_ == NULL)
    {
        ioniz.dedx_interpolant_ = new Interpolant(*dedx_interpolant_);
        dedx_interpolant_ = NULL;
    }


    if( dndx_interpolant_1d_ != NULL && ioniz.dndx_interpolant_1d_ != NULL)
    {
        dndx_interpolant_1d_->swap(*ioniz.dndx_interpolant_1d_);
    }
    else if( dndx_interpolant_1d_ == NULL && ioniz.dndx_interpolant_1d_ != NULL)
    {
        dndx_interpolant_1d_ = new Interpolant(*ioniz.dndx_interpolant_1d_);
        ioniz.dndx_interpolant_1d_ = NULL;
    }
    else if( dndx_interpolant_1d_ != NULL && ioniz.dndx_interpolant_1d_ == NULL)
    {
        ioniz.dndx_interpolant_1d_ = new Interpolant(*dndx_interpolant_1d_);
        dndx_interpolant_1d_ = NULL;
    }


    if( dndx_interpolant_2d_ != NULL && ioniz.dndx_interpolant_2d_ != NULL)
    {
        dndx_interpolant_2d_->swap(*ioniz.dndx_interpolant_2d_);
    }
    else if( dndx_interpolant_2d_ == NULL && ioniz.dndx_interpolant_2d_ != NULL)
    {
        dndx_interpolant_2d_ = new Interpolant(*ioniz.dndx_interpolant_2d_);
        ioniz.dndx_interpolant_2d_ = NULL;
    }
    else if( dndx_interpolant_2d_ != NULL && ioniz.dndx_interpolant_2d_ == NULL)
    {
        ioniz.dndx_interpolant_2d_ = new Interpolant(*dndx_interpolant_2d_);
        dndx_interpolant_2d_ = NULL;
    }


}

namespace PROPOSAL
{

ostream& operator<<(std::ostream& os, Ionization const &ioniz)
{
    os<<"---------------------------Ionization( "<<&ioniz<<" )---------------------------"<<std::endl;
    os<< static_cast <const CrossSections &>( ioniz ) << endl;
    os<< "------- Class Specific: " << endl;
    os<<"\tintegral:\t\t"<<ioniz.integral_ << endl;
    os<<endl;
    os<<"\tdo_dedx_Interpolation:\t\t"<<ioniz.do_dedx_Interpolation_<<endl;
    os<<"\tdedx_interpolant:\t\t"<<ioniz.dedx_interpolant_<<endl;
    os<<endl;
    os<<"\tdo_dndx_Interpolation:\t\t"<<ioniz.do_dndx_Interpolation_<<endl;
    os<<"\tdndx_interpolant_1d:\t\t"<<ioniz.dndx_interpolant_1d_<<endl;
    os<<"\tdndx_interpolant_2d:\t\t"<<ioniz.dndx_interpolant_2d_<<endl;
    os<<endl;
    os<<"\tTemp. variables: " << endl;
    os<<"-----------------------------------------------------------------------------------------------";
    return os;
}

}  // namespace PROPOSAL

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------Cross section / limit / private-----------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Ionization::CalculateStochasticLoss(double energy, double rnd)
{
    double rand, rsum;

    rand=medium_->GetSumCharge()*rnd;
    rsum=0;

    for(int i=0; i<medium_->GetNumComponents(); i++){

        Components::Component* component = medium_->GetComponents().at(i);
        rsum+=component->GetAtomInMolecule()* component->GetNucCharge();

        if(rsum>rand)
        {
            if(do_dndx_Interpolation_)
            {
                CrossSections::IntegralLimits limits = SetIntegralLimits(energy, 0);

                if(limits.vUp==limits.vMax)
                {
                    return energy*limits.vUp;
                }
                return energy*(limits.vUp*exp(dndx_interpolant_2d_->FindLimit(energy, rnd*sum_of_rates_)*log(limits.vMax/limits.vUp)));
            }
            else
            {
                return energy*integral_->GetUpperLimit();
            }
        }
    }

    log_fatal("m.totZ was not initialized correctly");

    return 0;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


CrossSections::IntegralLimits Ionization::SetIntegralLimits(double energy, int component)
{
    double aux;

    (void) component;
    CrossSections::IntegralLimits limits;

    double gamma   =   energy/particle_def_.mass;

    limits.vMin    =   (1.e-6*medium_->GetI())/energy;
    aux      =   ME/particle_def_.mass;
    limits.vMax    =   2*ME*(gamma*gamma-1)/((1 + 2*gamma*aux + aux*aux)*energy);
    limits.vMax    =   min(limits.vMax, 1. - particle_def_.mass/energy);

    if(limits.vMax<limits.vMin)
    {
        limits.vMax    =   limits.vMin;
    }

    limits.vUp =   min(limits.vMax, cut_settings_->GetCut(energy));

    if(limits.vUp<limits.vMin)
    {
        limits.vUp    =   limits.vMin;
    }

    return limits;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Ionization::Delta(double beta, double gamma)
{
    double X;

    X   =   log(beta * gamma)/log(10);

    if( X < medium_->GetX0())
    {
        return medium_->GetD0()*pow(10 , 2*(X - medium_->GetX0()));
    }
    else if(X < medium_->GetX1())
    {
        return 2*LOG10 * X + medium_->GetC()
                + medium_->GetA() * pow(medium_->GetX1() - X , medium_->GetM());
    }
    else
    {
        return 2*LOG10 * X + medium_->GetC();
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Ionization::D2Ndvdx(double energy, double v)
{
    double result, aux, aux2;

    CrossSections::IntegralLimits limits = SetIntegralLimits(energy, 0);

    //TODO(mario): Better way? Sat 2017/09/02
    double square_momentum = energy * energy - particle_def_.mass * particle_def_.mass;
    double particle_momentum = sqrt(max(square_momentum, 0.0));
    double beta    =   particle_momentum/energy;
    double gamma   =   energy/particle_def_.mass;

    aux  = beta * beta;
    aux2 = v / (1 + 1 / gamma);
    aux2 *= 0.5 * aux2;
    result = 1 - aux * (v / limits.vMax) + aux2;
    result *= IONK * particle_def_.charge * particle_def_.charge * medium_->GetZA() /
              (2 * aux * energy * v * v);

    return medium_->GetMassDensity() * result;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Ionization::InelCorrection(double energy, double v)
{
    double result, a, b, c;
    CrossSections::IntegralLimits limits = SetIntegralLimits(energy, 0);

    double gamma = energy / particle_def_.mass;

    a      = log(1 + 2 * v * energy / ME);
    b      = log((1 - v / limits.vMax) / (1 - v));
    c      = log((2 * gamma * (1 - v) * ME) / (particle_def_.mass * v));
    result = a * (2 * b + c) - b * b;

    return (ALPHA / (2 * PI)) * result;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------Functions to interpolate---------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Ionization::FunctionToBuildDEdxInterpolant( double energy)
{
    return CalculatedEdx(energy);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Ionization::FunctionToBuildDNdxInterpolant(double energy)
{
    return dndx_interpolant_2d_->Interpolate(energy, 1.0);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Ionization::FunctionToBuildDNdxInterpolant2D( double energy , double v)
{
    CrossSections::IntegralLimits limits = SetIntegralLimits(energy, 0);;

    if(limits.vUp==limits.vMax){
        return 0;
    }

    v=limits.vUp*exp(v*log(limits.vMax/limits.vUp));

    double dNdx = integral_->Integrate(limits.vUp, v,boost::bind(&Ionization::FunctionToDNdxIntegral, this,  energy, _1),3,1);

    return dNdx;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------Functions to integrate----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Ionization::FunctionToDEdxIntegral(double energy, double variable)
{
    return variable*D2Ndvdx(energy, variable)*InelCorrection(energy, variable);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Ionization::FunctionToDNdxIntegral(double energy, double variable)
{
    return multiplier_ * D2Ndvdx(energy, variable) * (1+InelCorrection(energy, variable));
}



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

// void Ionization::SetParametrization(ParametrizationType::Enum parametrization){
//     parametrization_ = parametrization;
//     log_warn("This has no effect. Till now only one parametrization for Ionization implemented");
//     if (parametrization_ != ParametrizationType::IonizBetheBloch)
//         log_warn("The parametrization type number '%i' is different to the one that is implemented with type number '%i' "
//             , parametrization_, ParametrizationType::IonizBetheBloch);
// }
