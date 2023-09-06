/**
 * @file IVectorWrapper.hpp
 * @brief Interface to the vector field wrappers used by the diagnostics
 */

#ifndef QUICC_DIAGNOSTICS_IVECTORWRAPPER_HPP
#define QUICC_DIAGNOSTICS_IVECTORWRAPPER_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/ScalarFields/ScalarField.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief Interface to the vector field wrappers used by the diagnostics
    */
   class IVectorWrapper
   {
      public:
         /**
          * @brief Constructor
          */
         IVectorWrapper();

         /**
          * @brief Destructor
          */
         virtual ~IVectorWrapper();

         /**
          * @brief Get first vector field component
          */
         virtual const Framework::Selector::PhysicalScalarField& one() const = 0;

         /**
          * @brief Get second vector field component
          */
         virtual const Framework::Selector::PhysicalScalarField& two() const = 0;

         /**
          * @brief Get third vector field component
          */
         virtual const Framework::Selector::PhysicalScalarField& three() const = 0;

         /**
          * @brief Get Resolution
          */
         virtual const Resolution& res() const = 0;

      protected:

      private:
   };

   /// Typedef for a shared IVectorWrapper
   typedef std::shared_ptr<IVectorWrapper> SharedIVectorWrapper;
}
}

#endif // QUICC_DIAGNOSTICS_IVECTORWRAPPER_HPP
