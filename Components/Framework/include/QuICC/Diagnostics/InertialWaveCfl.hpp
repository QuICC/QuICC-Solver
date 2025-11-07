/**
 * @file InertialWaveCfl.hpp
 * @brief CFL constraint in a spherical geometry
 */

#ifndef QUICC_DIAGNOSTICS_INERTIALWAVECFL_HPP
#define QUICC_DIAGNOSTICS_INERTIALWAVECFL_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/NonDimensional/INumber.hpp"
#include "QuICC/Diagnostics/ICflWrapper.hpp"
#include "QuICC/Diagnostics/IVectorWrapper.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief CFL constraint in a spherical geometry
    */
   class InertialWaveCfl: public ICflWrapper
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spVelocity Velocity wrapper
          */
         InertialWaveCfl(const std::map<std::size_t,NonDimensional::SharedINumber>& params, const MHDFloat courant);

         /**
          * @brief Destructor
          */
         virtual ~InertialWaveCfl() = default;

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
          * @brief CFL conditions
          */
         Matrix mCfl;
   };

   /// Typedef for a shared InertialWaveCfl
   typedef std::shared_ptr<InertialWaveCfl> SharedInertialWaveCfl;
}
}

#endif // QUICC_DIAGNOSTICS_INERTIALWAVECFL_HPP
