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
#include "QuICC/Io/Variable/StateFileWriter.hpp"
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

   /**
    * @brief Save eigenfunctions to statefile
    */
   void saveEigenfunction(const int m, const MHDComplex ev,
      const std::vector<MHDComplex>& ef);

   /**
    * @brief state file
    */
   std::shared_ptr<Io::Variable::StateFileWriter> mpH5File;

   /**
    * @brief Scalar fields in eigenfunction
    */
   std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>
      mScalars;

   /**
    * @brief Vector fields in eigenfunction
    */
   std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>
      mVectors;
};

} // namespace QuICC

#endif // QUICC_MARGINALCURVE_HPP
