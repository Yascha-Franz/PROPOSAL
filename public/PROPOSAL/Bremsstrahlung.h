#include <vector>
#include "CrossSections.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Interpolant.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/Particle.h"

class Bremsstrahlung: public CrossSections
{
private:

    bool        lorenz_;        /// enable lorenz cut
    double      lorenz_cut_;  	/// in [MeV] // - set to 1.e6 in Constructor
    int         component_;     /// nucleon in the medium on which the bremsstahlung occur

    Integral*   dedx_integral_;
    Interpolant* dedx_interpolant_;

    std::vector<Integral*>    dndx_integral_;
    std::vector<Interpolant*> dndx_interpolant_1d_; //Stochastic dNdx()
    std::vector<Interpolant*> dndx_interpolant_2d_; //Stochastic dNdx()

    double      eLpm_;



//----------------------------------------------------------------------------//
    //Memberfunctions

    // The following functions define the different
    // parametrization for the Bremsstrahlung Cross Section
    /*!
    \param v relativ energy loss
    \param i the nucleon in the medium on which the bremsstahlung occur
    \return  Calculates \f$ a_{1} \f$ see function Sel in this class
    */
    double KelnerKakoulinPetrukhinParametrization(double v, int i);

//----------------------------------------------------------------------------//
    /*!
    \param v relativ energy loss
    \param i the nucleon in the medium on which the bremsstahlung occur
    \return Calculates \f$ a_{1} \f$ see function Sel in this class
    */

    double AndreevBezrukovBugaevParametrization(double v, int i);

//----------------------------------------------------------------------------//
    /*!
    \param v relativ energy loss
    \param i the nucleon in the medium on which the bremsstahlung occur
    \return Calculates \f$ a_{1} \f$ see function Sel in this class
    */
    double PetrukhinShestakovParametrization(double v, int i);

//----------------------------------------------------------------------------//
    /*!
    \param v relativ energy loss
    \param i the nucleon in the medium on which the bremsstahlung occur
    \return Calculates \f$ a_{1} \f$ see function Sel in this class
    */
    double CompleteScreeningCase(double v, int i);

//----------------------------------------------------------------------------//

    /*!
    Landau Pomeranchuk Migdal effect and dielectric suppression evaluation
    lpm is evaluated:
    \f[ lpm= return=\frac{x_i}{3}\Big[v^2\frac{G(s)}{\gamma^2}+
    2(1+(1-v)^2)\frac{\Phi(s)}{\gamma}\Big] \f]
    \param v relativ energy loss
    \param s1
    \return the lpm correction factor
    */

    double lpm(double v, double s1);

//----------------------------------------------------------------------------//

    double FunctionToDEdxIntegral(double variable);

//----------------------------------------------------------------------------//

    double FunctionToDNdxIntegral(double variable);

//----------------------------------------------------------------------------//

    /*!
    this is what the Elastic Bremsstrahlung Cross Section (EBCS) is equal to
    units are [1/cm] since the multiplication by No*n is done here.
    Corrections for excitations of the nucleus and deep inelastic
    excitations of separate nucleons are
    included (positive term dn/Z), as well as the contribution of
    the mu-diagrams to the inelastic
    bremsstrahlung on the electrons (non-zero only for allowed
    energies of photon after electron recoil).
    four different parametrizations are enumerated, the final result is:
    \f[ \sigma_{el}= \rho_{mol}N(z^2)^2\Big[a_1\Big] \f] and \f$ a_{1} \f$
    depends on the chosen parametrization
    \param i the nucleon in the medium on which the bremsstahlung occur
    \param v relativ energy loss
    \return Elastic Bremsstrahlung Cross Section [1/cm]
    */

    double ElasticBremsstrahlungCrossSection(double v, int i);

//----------------------------------------------------------------------------//

    void SetIntegralLimits(int component);

//----------------------------------------------------------------------------//

public:

//----------------------------------------------------------------------------//

    Bremsstrahlung();
    Bremsstrahlung(const Bremsstrahlung&);
    Bremsstrahlung& operator=(const Bremsstrahlung&);
    Bremsstrahlung(Particle* particle, Medium* medium, EnergyCutSettings* cut_settings);


//----------------------------------------------------------------------------//

    double CalculatedEdx();

//----------------------------------------------------------------------------//

    double CalculatedNdx();


//----------------------------------------------------------------------------//

    double CalculatedNdx(double rnd);

//----------------------------------------------------------------------------//

    double CalculateStochasticLoss();

//----------------------------------------------------------------------------//

    void EnableDNdxInterpolation();

//----------------------------------------------------------------------------//

    void EnableDEdxInterpolation();

//----------------------------------------------------------------------------//

    void DisableDNdxInterpolation();

//----------------------------------------------------------------------------//

    void DisableDEdxInterpolation();

//----------------------------------------------------------------------------//

    double FunctionToBuildDEdxInterpolant(double energy);

//----------------------------------------------------------------------------//

    double FunctionToBuildDNdxInterpolant(double energy);

//----------------------------------------------------------------------------//

    double FunctionToBuildDNdxInterpolant2D(double energy , double v);


    //Getter

    bool GetLorenz() const
    {
        return lorenz_;
    }

//----------------------------------------------------------------------------//
    double GetLorenzCut() const
    {
        return lorenz_cut_;
    }
//----------------------------------------------------------------------------//
    //Setter

    void SetLorenz(bool lorenz);
    void SetLorenzCut(double lorenz_cut);
//----------------------------------------------------------------------------//

    ~Bremsstrahlung(){}

};
