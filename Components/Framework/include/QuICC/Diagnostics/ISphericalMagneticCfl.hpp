/**
 * @file ISphericalMagneticCfl.hpp
 * @brief CFL constraint in a spherical geometry
 */

#ifndef QUICC_DIAGNOSTICS_ISPHERICALMAGNETICCFL_HPP
#define QUICC_DIAGNOSTICS_ISPHERICALMAGNETICCFL_HPP

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
   class ISphericalMagneticCfl: public ISphericalCflWrapper
   {
      public:
         /**
          * @brief Constructor
          */
         ISphericalMagneticCfl(const std::map<std::size_t,NonDimensional::SharedINumber>& params, const MHDFloat courant);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalMagneticCfl() = default;

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
          * @brief Alfven wave scale
          */
         const MHDFloat mcAlfvenScale;

         /**
          * @brief Alfven wave damping
          */
         const MHDFloat mcAlfvenDamping;

         /**
          * @brief Velocity ID
          */
         std::size_t mVelId;

         /**
          * @brief Magnetic ID
          */
         std::size_t mMagId;
   };

   /// Typedef for a shared ISphericalMagneticCfl
   typedef std::shared_ptr<ISphericalMagneticCfl> SharedISphericalMagneticCfl;
}
}

#endif // QUICC_DIAGNOSTICS_ISPHERICALMAGNETICCFL_HPP
