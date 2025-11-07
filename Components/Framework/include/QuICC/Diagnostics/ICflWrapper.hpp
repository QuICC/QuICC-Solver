/**
 * @file ICflWrapper.hpp
 * @brief Interface for the CFL constraint
 */

#ifndef QUICC_DIAGNOSTICS_ICFLWRAPPER_HPP
#define QUICC_DIAGNOSTICS_ICFLWRAPPER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Diagnostics/ICflWrapper.hpp"
#include "QuICC/Diagnostics/IVectorWrapper.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief Interface for the CFL constraint
    */
   class ICflWrapper
   {
      public:
         /**
          * @brief Constructor
          */
         ICflWrapper(const MHDFloat courant);

         /**
          * @brief Destructor
          */
         virtual ~ICflWrapper() = default;

         /**
          * @brief Required fields
          */
         std::vector<std::size_t> fieldIds() const;

         /**
          * @brief Set field
          */
         void setField(const std::size_t id, const SharedIVectorWrapper spField);

         /**
          * @brief Initialize wrapper
          */
         virtual void init(const std::vector<Array>& mesh) = 0;

         /**
          * @brief Cfl constraint is active
          */
         bool isActive() const;

         /**
          * @brief Get initial CFL constraint
          */
         virtual Matrix initialCfl() const = 0;

         /**
          * @brief Get CFL constraint
          */
         virtual Matrix cfl() const = 0;

      protected:
         /**
          * @brief Courant constant used for the CFL computation
          */
         const MHDFloat mcCourant;

         /**
          * @brief CFL contraint is active?
          */
         bool mIsActive;

         /**
          * @brief Shared field wrappers
          */
         std::map<std::size_t,SharedIVectorWrapper> mFields;

      private:
   };

   /// Typedef for a shared ICflWrapper
   typedef std::shared_ptr<ICflWrapper> SharedICflWrapper;

} // namespace Diagnostics
} // namespace QuICC

#endif // QUICC_DIAGNOSTICS_ICFLWRAPPER_HPP
