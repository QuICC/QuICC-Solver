/**
 * @file ISphericalMagneticInviscidCfl.hpp
 * @brief CFL constraint in a spherical geometry
 */

#ifndef QUICC_DIAGNOSTICS_ISPHERICALMAGNETICINVISCIDCFL_HPP
#define QUICC_DIAGNOSTICS_ISPHERICALMAGNETICINVISCIDCFL_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/NonDimensional/INumber.hpp"
#include "QuICC/Diagnostics/ISphericalInviscidCflWrapper.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief CFL constraint in a spherical geometry
    */
   class ISphericalMagneticInviscidCfl: public ISphericalInviscidCflWrapper
   {
      public:
         /**
          * @brief Constructor
          */
         ISphericalMagneticInviscidCfl(const std::map<std::size_t,NonDimensional::SharedINumber>& params, const MHDFloat courant);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalMagneticInviscidCfl() = default;

         /**
          * @brief Define velocity field
          */
         void defineVelocity(const std::size_t velId);

         /**
          * @brief Define magnetic field
          */
         void defineMagnetic(const std::size_t magId);

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

         /**
          * @brief Magnetic ID
          */
         std::size_t mMagId;
   };

   /// Typedef for a shared ISphericalMagneticInviscidCfl
   typedef std::shared_ptr<ISphericalMagneticInviscidCfl> SharedISphericalMagneticInviscidCfl;
}
}

#endif // QUICC_DIAGNOSTICS_ISPHERICALMAGNETICINVISCIDCFL_HPP
