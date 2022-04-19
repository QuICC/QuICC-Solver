/**
 * @file ICflWrapper.hpp
 * @brief Interface for the CFL constraint
 */

#ifndef QUICC_DIAGNOSTICS_ICFLWRAPPER_HPP
#define QUICC_DIAGNOSTICS_ICFLWRAPPER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Diagnostics/ICflWrapper.hpp"
#include "QuICC/Diagnostics/IVectorWrapper.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief Interface for the CFL constraint
    */
   class ICflWrapper
   {
      public:
         /**
          * @brief Constructor
          *
          * @param Velocity wrapper
          */
         ICflWrapper(const SharedIVectorWrapper spVelocity);

         /**
          * @brief Constructor
          *
          * @param Velocity wrapper
          * @param Magnetic wrapper
          */
         ICflWrapper(const SharedIVectorWrapper spVelocity, const SharedIVectorWrapper spMagnetic);

         /**
          * @brief Destructor
          */
         virtual ~ICflWrapper();

         /**
          * @brief Initialize wrapper
          */
         virtual void init(const std::vector<Array>& mesh) = 0;

         /**
          * @brief Get initial CFL constraint
          */
         virtual Matrix initialCfl() const = 0;

         /**
          * @brief Get CFL constraint
          */
         virtual Matrix cfl() const = 0;

      protected:
         /**
          * @brief Shared velocity wrapper
          */
         SharedIVectorWrapper mspVelocity;

         /**
          * @brief Shared magnetic wrapper
          */
         SharedIVectorWrapper mspMagnetic;

         /**
          * @brief Update minimum CFL in matrix
          */
         void updateCflMatrix(Matrix& cfl) const;

      private:
   };

   /// Typedef for a shared ICflWrapper
   typedef std::shared_ptr<ICflWrapper> SharedICflWrapper;
}
}

#endif // QUICC_DIAGNOSTICS_ICFLWRAPPER_HPP
