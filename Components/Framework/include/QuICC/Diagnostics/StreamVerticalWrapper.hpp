/**
 * @file StreamVerticalWrapper.hpp
 * @brief Streamfunction and vertical velocity wrapper into velocity field
 */

#ifndef QUICC_DIAGNOSTICS_STREAMVERTICALWRAPPER_HPP
#define QUICC_DIAGNOSTICS_STREAMVERTICALWRAPPER_HPP

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
#include "QuICC/Diagnostics/IVectorWrapper.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief Streamfunction and vertical velocity wrapper into velocity field
    */
   class StreamVerticalWrapper: public IVectorWrapper
   {
      public:
         /**
          * @brief Constructor
          */
         StreamVerticalWrapper(const Framework::Selector::VariantSharedScalarVariable spStream, const Framework::Selector::VariantSharedScalarVariable spVertical);

         /**
          * @brief Destructor
          */
         virtual ~StreamVerticalWrapper();

         /**
          * @brief Get first velocity field component
          */
         virtual const Framework::Selector::PhysicalScalarField& one() const;

         /**
          * @brief Get second velocity field component
          */
         virtual const Framework::Selector::PhysicalScalarField& two() const;

         /**
          * @brief Get third velocity field component
          */
         virtual const Framework::Selector::PhysicalScalarField& three() const;

         /**
          * @brief Get Resolution
          */
         virtual const Resolution& res() const;

      protected:

      private:
         /**
          * @brief Shared streamfunction variable
          */
         Framework::Selector::VariantSharedScalarVariable mspStream;

         /**
          * @brief Shared vertical velocity variable
          */
         Framework::Selector::VariantSharedScalarVariable mspVertical;
   };

   /// Typedef for a shared StreamVerticalWrapper
   typedef std::shared_ptr<StreamVerticalWrapper> SharedStreamVerticalWrapper;
}
}

#endif // QUICC_DIAGNOSTICS_STREAMVERTICALWRAPPER_HPP
