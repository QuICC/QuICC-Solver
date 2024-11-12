/**
 * @file MarginalCurve.cpp
 * @brief Source of the high level simulation
 */

// System includes
//
#include <algorithm>
#include <boost/math/tools/roots.hpp>
#include <iostream>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/QuICCTimer.hpp"
#include "QuICC/Timers/StageTimer.hpp"
#include "Stability/LinearStability.hpp"
#include "Stability/MarginalCurve.hpp"
#include "Types/Math.hpp"

namespace QuICC {

/**
* @brief Functor to wrap calls for stability calculations
*/
template <unsigned int TLevel> class func
{
public:
   /**
    * @brief Constructor
    *
    * @param spStab Linear stability solver
    * @param spLog  Logging stream
    */
   func(std::shared_ptr<LinearStability> spStab, std::shared_ptr<std::ofstream> spLog) : mspStab(spStab), mspLog(spLog) {};

   /**
    * @brief Destructor
    */
   ~func() = default;

   /**
    * @brief Interface to call solver
    */
   MHDFloat operator()(const MHDFloat Ra)
   {
      if constexpr(TLevel == 0)
      {
         return (*this->mspStab)(Ra);
      }
      else if constexpr(TLevel == 1)
      {
         const unsigned int nev = 5;
         std::vector<MHDComplex> evs(nev);
         auto ret = (*this->mspStab)(Ra, evs);

         this->print(Ra, evs);

         return ret;
      }
      else
      {
         throw std::logic_error("Unknown output level");
      }
   };

   /**
    * @brief Logging computed eigenvalues
    */
   void print(const MHDFloat Ra, const std::vector<MHDComplex>& evs) const
   {
      std::ofstream& log = *this->mspLog;
      log << "Ra = " << Ra << std::endl;
      for (unsigned int i = 0; i < evs.size(); i++)
      {
         log << " \t growth: " << evs.at(i) << std::endl;
      }
   }

protected:
   /**
    * @brief Linear stability calculator
    */
   std::shared_ptr<LinearStability> mspStab;

   /**
    * @brief Logging stream
    */
   std::shared_ptr<std::ofstream> mspLog;
};

MarginalCurve::MarginalCurve() : StabilityBase()
{
   // Initialize SLEPc/PETSc
   PetscCallVoid(SlepcInitializeNoArguments());
}

MarginalCurve::~MarginalCurve()
{
   // Finalize SLEPc/PETSc
   PetscCallVoid(SlepcFinalize());
}

void MarginalCurve::mainRun()
{
   std::vector<MHDFloat> eigs = {static_cast<MHDFloat>(
      this->mspRes->sim().dim(Dimensions::Simulation::SIM3D,
         Dimensions::Space::SPECTRAL) -
      1)};

   auto spLinStab = std::make_shared<LinearStability>(eigs, this->mspRes,
      this->mspEqParams->map(), this->createBoundary()->map(),
      this->mspBackend);

   MHDFloat guess = this->mspEqParams->nd(NonDimensional::Rayleigh::id());
   MHDFloat factor = 1.1; // To multiply
   int digits =
      std::numeric_limits<MHDFloat>::digits; // Maximum possible binary digits
                                             // accuracy for type T.
   // digits used to control how accurate to try to make the result.
   int get_digits = (digits * 3) / 4; // Near maximum (3/4) possible accuracy.

   const boost::uintmax_t maxit = 20;
   boost::uintmax_t it =
      maxit; // Initally our chosen max iterations, but updated with actual.
   // We could also have used a maximum iterations provided by any policy:
   // boost::uintmax_t max_it = policies::get_max_root_iterations<Policy>();
   bool is_rising =
      true; // So if result if guess^3 is too low, try increasing guess.
   boost::math::tools::eps_tolerance<double> tol(get_digits);

   auto spLog = std::make_shared<std::ofstream>("marginal.log");

   std::pair<MHDFloat, MHDFloat> r = boost::math::tools::bracket_and_solve_root(
      func<1>(spLinStab, spLog), guess, factor, is_rising, tol, it);

   std::ofstream& logger = *spLog;

   unsigned int prec = (std::numeric_limits<MHDFloat>::digits10*3)/4 + 1;
   logger << std::string(50, '-') << std::endl 
      << "Critical Rayleigh number converged to the bracket: " << std::endl
      << std::setprecision(prec)  << r.first << " < Ra < " << r.second << std::endl;
}

void MarginalCurve::preRun() {}

void MarginalCurve::postRun() {}

} // namespace QuICC
