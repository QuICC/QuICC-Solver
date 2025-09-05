/**
 * @file CartesianCfl.hpp
 * @brief CFL constraint in a Cartesian geometry
 */

#ifndef QUICC_DIAGNOSTICS_CARTESIANCFL_HPP
#define QUICC_DIAGNOSTICS_CARTESIANCFL_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Diagnostics/ICflWrapper.hpp"
#include "QuICC/NonDimensional/INumber.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief CFL constraint in a Cartesian geometry
    */
   class CartesianCfl: public ICflWrapper
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spVelocity Velocity wrapper
          */
         CartesianCfl(const std::map<std::size_t,NonDimensional::SharedINumber>& params, const MHDFloat courant);

         /**
          * @brief Destructor
          */
         virtual ~CartesianCfl() = default;

         /**
          * @brief Define velocity field
          */
         void defineVelocity(const std::size_t velId);

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
          * @brief Velocity ID
          */
         std::size_t mVelId;

         /**
          * @brief Spacing between grid points
          */
         std::vector<Array> mMeshSpacings;
   };

   /// Typedef for a shared CartesianCfl
   typedef std::shared_ptr<CartesianCfl> SharedCartesianCfl;
}
}

#endif // QUICC_DIAGNOSTICS_CARTESIANCFL_HPP
