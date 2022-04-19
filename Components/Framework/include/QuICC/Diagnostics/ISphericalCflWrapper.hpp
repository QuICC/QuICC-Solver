/**
 * @file ISphericalCflWrapper.hpp
 * @brief CFL constraint in a spherical geometry
 */

#ifndef QUICC_DIAGNOSTICS_ISPHERICALCFLWRAPPER_HPP
#define QUICC_DIAGNOSTICS_ISPHERICALCFLWRAPPER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

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
   class ISphericalCflWrapper: public ICflWrapper
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spVelocity Velocity wrapper
          */
         ISphericalCflWrapper(const SharedIVectorWrapper spVelocity, const std::map<std::size_t,NonDimensional::SharedINumber>& params);

         /**
          * @brief Constructor
          *
          * @param spVelocity Velocity wrapper
          * @param spMagnetic Magnetic wrapper
          */
         ISphericalCflWrapper(const SharedIVectorWrapper spVelocity, const SharedIVectorWrapper spMagnetic, const std::map<std::size_t,NonDimensional::SharedINumber>& params);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalCflWrapper();

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
          * @brief Get effective max harmonic degree L
          */
         virtual MHDFloat effectiveMaxL(const MHDFloat r) const = 0;

         /**
          * @brief Initialise the mesh spacings
          */
         void initMesh(const std::vector<Array>& mesh);

         /**
          * @brief Courant constant used for the CFL computation
          */
         const MHDFloat mcCourant;

         /**
          * @brief Alfven wave scale
          */
         const MHDFloat mcAlfvenScale;

         /**
          * @brief Alfven wave damping
          */
         const MHDFloat mcAlfvenDamping;

         /**
          * @brief CFL conditions
          */
         Array mGlobalCfl;

         /**
          * @brief Spacing between grid points
          */
         std::vector<Array> mMeshSpacings;
   };

   /// Typedef for a shared ISphericalCflWrapper
   typedef std::shared_ptr<ISphericalCflWrapper> SharedISphericalCflWrapper;
}
}

#endif // QUICC_DIAGNOSTICS_ISPHERICALCFLWRAPPER_HPP
