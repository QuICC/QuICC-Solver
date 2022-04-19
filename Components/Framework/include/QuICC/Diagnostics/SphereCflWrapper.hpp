/**
 * @file SphereCflWrapper.hpp
 * @brief CFL constraint in a full sphere
 */

#ifndef QUICC_DIAGNOSTICS_SPHERECFLWRAPPER_HPP
#define QUICC_DIAGNOSTICS_SPHERECFLWRAPPER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Diagnostics/ISphericalCflWrapper.hpp"
#include "QuICC/Diagnostics/IVectorWrapper.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief CFL constraint in a full sphere
    */
   class SphereCflWrapper: public ISphericalCflWrapper
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spVelocity Velocity wrapper
          */
         SphereCflWrapper(const SharedIVectorWrapper spVelocity, const std::map<std::size_t,NonDimensional::SharedINumber>& params);

         /**
          * @brief Constructor
          *
          * @param spVelocity Velocity wrapper
          * @param spMagnetic Magnetic wrapper
          */
         SphereCflWrapper(const SharedIVectorWrapper spVelocity, const SharedIVectorWrapper spMagnetic, const std::map<std::size_t,NonDimensional::SharedINumber>& params);

         /**
          * @brief Destructor
          */
         virtual ~SphereCflWrapper();

      protected:

      private:
         /**
          * @brief Get effective max harmonic degree L
          */
         virtual MHDFloat effectiveMaxL(const MHDFloat r) const;

         /**
          * @brief Compute approximation to first Jacobi root for CFL
          */
         MHDFloat jacobiRoot(const MHDFloat l) const;
   };

   /// Typedef for a shared SphereCflWrapper
   typedef std::shared_ptr<SphereCflWrapper> SharedSphereCflWrapper;
}
}

#endif // QUICC_DIAGNOSTICS_SPHERECFLWRAPPER_HPP
