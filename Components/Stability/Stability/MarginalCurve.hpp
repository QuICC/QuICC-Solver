/**
 * @file MarginalCurve.hpp
 * @brief High level implementation of a stability solver base
 */

#ifndef QUICC_MARGINALCURVE_HPP
#define QUICC_MARGINALCURVE_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Equations/EquationParameters.hpp"
#include "Stability/StabilityBase.hpp"

namespace QuICC {

/**
 * @brief High level implementation of a base for the simulations
 */
class MarginalCurve : public StabilityBase
{
public:
   /**
    * @brief Constructor
    */
   MarginalCurve();

   /**
    * @brief Simple empty destructor
    */
   virtual ~MarginalCurve();

protected:
private:
   /**
    * @brief Do operations required just before starting the main loop
    */
   virtual void preRun() override;

   /**
    * @brief Do operations required during the main loop
    */
   virtual void mainRun() override;

   /**
    * @brief Do operations required just after finishing the main loop
    */
   virtual void postRun() override;
};

/// Typedef for a shared pointer of a MarginalCurve
typedef std::shared_ptr<MarginalCurve> SharedMarginalCurve;
} // namespace QuICC

#endif // QUICC_MARGINALCURVE_HPP
