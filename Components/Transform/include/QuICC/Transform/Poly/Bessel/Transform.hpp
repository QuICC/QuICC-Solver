/**
 * @file Transform.hpp
 * @brief Implementation of the spherical Bessel transform
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_TRANSFORM_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_TRANSFORM_HPP


// System includes
//
#include <set>
#include <map>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/NonDimensional/INumber.hpp"
#include "QuICC/Transform/Poly/Bessel/IBesselOperator.hpp"
#include "QuICC/Transform/ITransformMap.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

   /**
    * @brief Implementation of the spherical Bessel transform
    */
   class Transform
   {
      public:
         /// Typedef for the configuration class
         typedef IBesselOperator::SetupType SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef IBesselOperator::SharedSetupType SharedSetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef ITransformMap<IBesselOperator> MapFunctor;

         /// Typedef for the configuration class as a shared pointer
         typedef ITransformMap<IBesselOperator>::MapType MapType;

         /**
          * @brief Constructor
          */
         Transform() = default;

         /**
          * @brief Destructor
          */
         virtual ~Transform() = default;

         /**
          * @brief Initialise the polynomial transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedSetupType spSetup);

         /**
          * @brief Map transform operator to ID
          */
         void addOperator(const MapFunctor& f);

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
         void transform(MatrixZ& rOut, const MatrixZ& in, const IBesselOperator& op);

         /**
          * @brief Compute transform (Complex -> Real) using given operator
          *
          * @param rOut Output values
          * @param in   Input values
          * @param op   Transform operator
          */
         void transform(Matrix& rOut, const MatrixZ& in, const IBesselOperator& op);

         /**
          * @brief Compute transform (Complex -> Complex) mapped to ID
          *
          * @param rOut Output values
          * @param in   Input values
          * @param id   ID of transform operator
          */
         void transform(MatrixZ& rOut, const MatrixZ& in, const std::size_t id);

         /**
          * @brief Compute transform (Complex -> Real) mapped to ID
          *
          * @param rOut Output values
          * @param in   Input values
          * @param id   ID of transform operator
          */
         void transform(Matrix& rOut, const MatrixZ& in, const std::size_t id);

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
         Internal::Array mIGrid;

         /**
          * @brief Storage for the quadrature weights with internal precision
          */
         Internal::Array mIWeights;

         /**
          * @brief Polynomial setup object providing the sizes
          */
         SharedSetup    mspSetup;

         /**
          * @brief Store transform operator to ID mapping
          */
         MapType mOps;
   };

}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_TRANSFORM_HPP
