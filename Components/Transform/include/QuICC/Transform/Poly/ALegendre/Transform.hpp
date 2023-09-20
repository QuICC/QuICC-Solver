/**
 * @file Transform.hpp
 * @brief Implementation of the associated Legendre transform
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_TRANSFORM_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_TRANSFORM_HPP

// Debug includes
//

// Configuration includes
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"

// System includes
//
#include <set>
#include <map>

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/NonDimensional/INumber.hpp"
#include "QuICC/Transform/Poly/ALegendre/IALegendreOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

   /**
    * @brief Implementation of the associated Legendre transform
    */
   class Transform
   {
      public:
         /// Typedef for the configuration class
         typedef IALegendreOperator::SetupType SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef IALegendreOperator::SharedSetupType SharedSetupType;

         /**
          * @brief Constructor
          */
         Transform();

         /**
          * @brief Destructor
          */
         virtual ~Transform();

         /**
          * @brief Initialise the polynomial transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedSetupType spSetup);

         /**
          * @brief Map transform operator to ID
          */
         template <typename TOp> void addOperator(const std::size_t id = 0);

         /**
          * @brief set list of required options
          */
         void requiredOptions(std::set<std::size_t>& list, const Dimensions::Transform::Id dimId) const;

         /**
          * @brief Set the required options
          */
         void setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options, const Dimensions::Transform::Id dimId);

         /**
          * @brief Get the physical grid
          */
         Array meshGrid() const;

         /**
          * @brief Compute transform (Complex -> Complex) using given operator
          *
          * @param rOut Output values
          * @param in   Input values
          * @param op   Transform operator
          */
         void transform(MatrixZ& rOut, const MatrixZ& in, const ITransformOperator& op);

         /**
          * @brief Compute transform (Complex -> Complex) mapped to ID
          *
          * @param rOut Output values
          * @param in   Input values
          * @param op   Transform operator
          */
         void transform(MatrixZ& rOut, const MatrixZ& in, const std::size_t id);

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const;

      protected:

      private:
         /**
          * @brief Generate a physical grid
          */
         static Array generateGrid(const int size);

         /**
          * @brief Initialise the quadrature
          */
         void initQuadrature();

         /**
          * @brief Storage for the quadrature points x = [-1, 1] with internal precision
          */
         internal::Array mIGrid;

         /**
          * @brief Storage for the quadrature weights with internal precision
          */
         internal::Array mIWeights;

         /**
          * @brief Storage for the quadrature points x = [-1, 1] with internal precision
          */
         Array mThGrid;

         /**
          * @brief Polynomial setup object providing the sizes
          */
         SharedSetup    mspSetup;

         /**
          * @brief Store transform operator to ID mapping
          */
         std::map<std::size_t,std::shared_ptr<ITransformOperator> > mOps;
   };

   template <typename TOp> void Transform::addOperator(const std::size_t id)
   {
      std::size_t lid;

      if(id == 0)
      {
         lid = this->mOps.size();
      } else
      {
         lid = id;
      }

      std::shared_ptr<TOp> spOp = std::make_shared<TOp>();
      auto ret = this->mOps.try_emplace(lid, spOp);

      if(!ret.second)
      {
         throw std::logic_error("ALegendre operator could not be added to transform");
      }
   }

}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_TRANSFORM_HPP
