/**
 * @file ISphericalHydroCfl.hpp
 * @brief CFL constraint in a spherical geometry
 */

#ifndef QUICC_DIAGNOSTICS_ISPHERICALHYDROCFL_HPP
#define QUICC_DIAGNOSTICS_ISPHERICALHYDROCFL_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/NonDimensional/INumber.hpp"
#include "QuICC/Diagnostics/ISphericalCflWrapper.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief CFL constraint in a spherical geometry
    */
   class ISphericalHydroCfl: public ISphericalCflWrapper
   {
      public:
         /**
          * @brief Constructor
          */
         ISphericalHydroCfl(const std::map<std::size_t,NonDimensional::SharedINumber>& params, const MHDFloat courant);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalHydroCfl() = default;

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
          * @brief Velocity ID
          */
         std::size_t mVelId;
   };

   /// Typedef for a shared ISphericalHydroCfl
   typedef std::shared_ptr<ISphericalHydroCfl> SharedISphericalHydroCfl;
}
}

#endif // QUICC_DIAGNOSTICS_ISPHERICALHYDROCFL_HPP
