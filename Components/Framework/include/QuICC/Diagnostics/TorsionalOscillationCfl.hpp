/**
 * @file TorsionalOscillationCfl.hpp
 * @brief CFL constraint in a spherical geometry
 */

#ifndef QUICC_DIAGNOSTICS_TORSIONALOSCIALLATIONCFL_HPP
#define QUICC_DIAGNOSTICS_TORSIONALOSCIALLATIONCFL_HPP

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
   class TorsionalOscillationCfl: public ICflWrapper
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spVelocity Velocity wrapper
          */
         TorsionalOscillationCfl(const std::map<std::size_t,NonDimensional::SharedINumber>& params, const MHDFloat courant);

         /**
          * @brief Destructor
          */
         virtual ~TorsionalOscillationCfl() = default;

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

   /// Typedef for a shared TorsionalOscillationCfl
   typedef std::shared_ptr<TorsionalOscillationCfl> SharedTorsionalOscillationCfl;
}
}

#endif // QUICC_DIAGNOSTICS_TORSIONALOSCIALLATIONCFL_HPP
