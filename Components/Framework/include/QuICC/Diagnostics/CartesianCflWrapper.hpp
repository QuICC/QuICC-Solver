/**
 * @file CartesianCflWrapper.hpp
 * @brief CFL constraint in a Cartesian geometry
 */

#ifndef QUICC_DIAGNOSTICS_CARTESIANCFLWRAPPER_HPP
#define QUICC_DIAGNOSTICS_CARTESIANCFLWRAPPER_HPP

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
    * @brief CFL constraint in a Cartesian geometry
    */
   class CartesianCflWrapper: public ICflWrapper
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spVelocity Velocity wrapper
          */
         CartesianCflWrapper(const SharedIVectorWrapper spVelocity);

         /**
          * @brief Destructor
          */
         virtual ~CartesianCflWrapper();

         /**
          * @brief Initialize wrapper
          */
         virtual void init(const std::vector<Array>& mesh);

         /**
          * @brief Get initial CFL constraint
          */
         virtual Matrix initialCfl() const;

         /**
          * @brief Get CFL constraint
          */
         virtual Matrix cfl() const;

      protected:

      private:
         /**
          * @brief Initialise the mesh spacings
          */
         void initMesh(const std::vector<Array>& mesh);

         /**
          * @brief Courant constant used for the CFL computation
          */
         const MHDFloat mcCourant;

         /**
          * @brief Spacing between grid points
          */
         std::vector<Array> mMeshSpacings;
   };

   /// Typedef for a shared CartesianCflWrapper
   typedef std::shared_ptr<CartesianCflWrapper> SharedCartesianCflWrapper;
}
}

#endif // QUICC_DIAGNOSTICS_CARTESIANCFLWRAPPER_HPP
