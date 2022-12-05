/**
 * @file IImExPCScheme.hpp
 * @brief Interface for a generic implicit/explicit predictor-corrector scheme
 */

#ifndef QUICC_TIMESTEP_IIMEXPCSCHEME_HPP
#define QUICC_TIMESTEP_IIMEXPCSCHEME_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Timestep/IScheme.hpp"

namespace QuICC {

namespace Timestep {

   /**
    * @brief Interface of generic of an implicit/explicit Runge-Kutta scheme from CB
    */
   class IImExPCScheme: public IScheme
   {
      public:
         /**
          * @brief Constructor
          */
         IImExPCScheme() = default;

         /**
          * @brief Destructor
          */
         virtual ~IImExPCScheme() = default;

         /**
          * @brief Implicit coefficient
          */
         virtual MHDFloat aIm(const int i) const = 0;

         /**
          * @brief Step fraction
          */
         virtual MHDFloat cEx(const int i) const = 0;

      protected:

      private:

   };

   /// Typedef for a shared pointer IImExPCScheme
   typedef std::shared_ptr<IImExPCScheme> SharedIImExPCScheme;

}
}

#endif // QUICC_TIMESTEP_IIMEXPCSCHEME_HPP
