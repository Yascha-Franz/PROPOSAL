#pragma once

#ifndef Epairproduction_H
#define Epairproduction_H

#include "PROPOSAL/crossection/CrossSections.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"

namespace PROPOSAL
{
    class Epairproduction;
}

std::ostream& operator<<(std::ostream& os, PROPOSAL::Epairproduction const &epair);

namespace PROPOSAL{

class Epairproduction: public CrossSections
{
protected:

    int component_;
    double v_;
    bool reverse_;
    double eLpm_;
    bool do_epair_interpolation_;


    Integral*   integral_;
    Integral*   integral_for_dEdx_;

    Interpolant* dedx_interpolant_;

    std::vector<Integral*>    dndx_integral_;
    std::vector<Interpolant*> dndx_interpolant_1d_; //Stochastic dNdx()
    std::vector<Interpolant*> dndx_interpolant_2d_; //Stochastic dNdx()
    std::vector<Interpolant*> epair_interpolant_;  //!< Interpolates function used by dNdx and dEdx


    std::vector<double> prob_for_component_; //!< probability for each medium component to interact with the particle (formerly H_)


//----------------------------------------------------------------------------//

    /*!
    this is the calculation of the d2Sigma/dvdRo - interface to Integral,
    the function which is defined here is:
    \f[ f(r) =return= \alpha^2r_e^2 \frac{2Z}{1,5\pi}(Z+k)
    \frac{1-v}{v}lpm(r^2,b,s)(F_e+\frac{m_e^2}{m_p^2}F_m)\f]
    */

    double FunctionToIntegral(double energy, double variable);

//----------------------------------------------------------------------------//

    /*!
    Landau Pomeranchuk Migdal effect evaluation,
    if the Landau Pomeranchuk Migdal effect is considered in
    the calculation, function is modified by a factor
    \f[lpm=return=\frac{(1+b)(A+(1+r^2)B)+b(C+(1+r^2)D)+(1-r^2)E}{[(2+r^2)(1+b)
    +x(3+r^2)]\ln\Big(1+\frac{1}{x}\Big)+\frac{1-r^2-b}{1+x}-(3+r^2)}\f]
    */

    double lpm(double energy, double r2, double b, double x);

//----------------------------------------------------------------------------//

    /*!
    this is the calculation of the dSigma/dv:
    \f[e_{Pair}=return = \rho N_Z z^2 \Big[ \int_{1-r_{max}}^{aux}
    f(r)dr + \int^1_{aux} f(r) dr \Big]\f]
    with \f$ aux=max(1-r_{Max}, ComputerPrecision)\f$ and \f$r_{max} =
    \sqrt{1-\frac{4m_e}{E_p v}}\Big(1-\frac{6m_p^2}{E_p^2(1-v)}\Big)\f$
    */

    double EPair(double energy, double v, int component);

//----------------------------------------------------------------------------//

    /*!
    pair production energy losses - interface to Integral;
    \f[ return= v\cdot e_{Pair}(v, components)\f]
    */

    double FunctionToDEdxIntegral(double energy, double variable);

//----------------------------------------------------------------------------//

    /*!
    Set the component of the primary;
    sets the parameters \f$v_{Min}=4 \frac{m_e}{E_p}\f$,
    \f$v_{Max}= max(-\frac{3}{4}\sqrt{e}\frac{m_p}{E_p}Z^{\frac{1}{3}},
    1-6\Big(\frac{m_p}{E_P}\Big)^2, 1-\frac{m_p}{E_p} )\f$
    and \f$v_{Up}=max( v_{max} , v_{Cut}(E_p))\f$
    */

    CrossSections::IntegralLimits SetIntegralLimits(double energy, int component);
//----------------------------------------------------------------------------//

    double FunctionToBuildEpairInterpolant( double energy , double v);

//----------------------------------------------------------------------------//
    double FunctionToBuildDEdxInterpolant( double energy);

//----------------------------------------------------------------------------//
    double FunctionToBuildDNdxInterpolant1D(double energy);

//----------------------------------------------------------------------------//
    double FunctionToBuildDNdxInterpolant2D( double energy, double v);

//----------------------------------------------------------------------------//
    double CalculateStochasticLoss(double energy, double rnd1);

//----------------------------------------------------------------------------//

public:

//----------------------------------------------------------------------------//

    Epairproduction();
    Epairproduction(const Epairproduction&);
    Epairproduction(PROPOSALParticle&, Medium* medium, EnergyCutSettings* cut_settings);
    bool operator==(const Epairproduction &epair) const;
    bool operator!=(const Epairproduction &epair) const;
    Epairproduction& operator=(const Epairproduction&);
    friend std::ostream& operator<<(std::ostream& os, Epairproduction const &epair);
//----------------------------------------------------------------------------//

    /*!
    contribution of pair production to -dE/dx,
    \f[\frac{dE}{dx}=return= multiplier \cdot E_p \sum_{numC}
    \int_{v_{Min}}^{v_{Up}} v e_{Pair}(v, components) dv\f]
    \return [E/cm]
    */

    double CalculatedEdx(double energy);

//----------------------------------------------------------------------------//

    double CalculatedNdx(double energy);

//----------------------------------------------------------------------------//

    double CalculatedNdx(double energy, double rnd);

//----------------------------------------------------------------------------//

    double CalculateStochasticLoss(double energy, double rnd1, double rnd2);

//----------------------------------------------------------------------------//

    void EnableDNdxInterpolation( std::string path ="", bool raw=false);

//----------------------------------------------------------------------------//

    void EnableDEdxInterpolation( std::string path ="", bool raw=false);

//----------------------------------------------------------------------------//

    void EnableEpairInterpolation( std::string path ="", bool raw=false);

//----------------------------------------------------------------------------//

    void DisableDNdxInterpolation();

//----------------------------------------------------------------------------//

    void DisableDEdxInterpolation();

//----------------------------------------------------------------------------//

    void DisableEpairInterpolation();

//----------------------------------------------------------------------------//

    double FunctionToDNdxIntegral(double energy, double variable);

//----------------------------------------------------------------------------//

    void swap(Epairproduction &epair);

//----------------------------------------------------------------------------//

    void ValidateOptions();

//----------------------------------------------------------------------------//
//-----------------------------Getter and Setter------------------------------//
//----------------------------------------------------------------------------//
    //Getter

    int GetComponent() const { return component_; }

    double GetLpm() const { return eLpm_; }

    bool GetReverse() const { return reverse_; }

    std::vector<double> GetProbForComponent() const {
        return prob_for_component_;
    }

    double GetV() const {
        return v_;
    }

	Interpolant* GetDedxInterpolant() const {
		return dedx_interpolant_;
	}

	std::vector<Integral*> GetDndxIntegral() const {
		return dndx_integral_;
	}

	std::vector<Interpolant*> GetDndxInterpolant1d() const {
		return dndx_interpolant_1d_;
	}

	std::vector<Interpolant*> GetDndxInterpolant2d() const {
		return dndx_interpolant_2d_;
	}

	Integral* GetIntegral() const {
		return integral_;
	}

	Integral* GetIntegralForDEdx() const {
		return integral_for_dEdx_;
	}

//----------------------------------------------------------------------------//
    //Setter

    // void SetParametrization(ParametrizationType::Enum parametrization = ParametrizationType::EPairKelnerKokoulinPetrukhin);

	void SetComponent(int component){ component_ = component; }

	void SetLpm(double lpm){ eLpm_ = lpm; }

	void SetReverse(bool reverse){ reverse_ = reverse; }

//----------------------------------------------------------------------------//
    //Destructor
    ~Epairproduction();

};

}

#endif //Epairproduction_H