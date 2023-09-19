/**
 * @file ImExPC2.hpp
 * @brief Implementation of an implicit/explicit predictor-corrector scheme of order 2
 */

#ifndef QUICC_TIMESTEP_IMEXPC2_HPP
#define QUICC_TIMESTEP_IMEXPC2_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Timestep/IImExPCScheme.hpp"
#include "QuICC/Timestep/SparseImExPCTimestepper.hpp"

namespace QuICC {

namespace Timestep {

   /**
    * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order 2
    */
   class ImExPC2: public IImExPCScheme
   {
      public:
         /// Typedef for underlying scheme implementation
         template <typename T1,typename T2,template <typename> class T3>
            using ImplementationType = SparseImExPCTimestepper<T1,T2,T3>;

         /**
          * @brief Constructor
          */
         ImExPC2();

         /**
          * @brief Destructor
          */
         virtual ~ImExPC2() = default;

         /**
          * @brief Butcher's tableau a_ij factor for implicit scheme
          */
         MHDFloat aIm(const int i) const final;

         /**
          * @brief Step fraction
          */
         virtual MHDFloat cEx(const int i) const final;

         /**
          * @brief Number of substeps for final step (this is +1 compared to theoretical value due to implementation)
          */
         int steps() const final;

         /**
          * @brief Order of the scheme
          */
         int order() const final;

         /**
          * @brief Scheme has embedded lower order scheme?
          */
         bool hasEmbedded() const final;

         /**
          * @brief Name of the scheme
          */
         std::string name() const final;

      protected:
         /**
          * @brief Initialize coefficients
          */
         void init();

         /**
          * @brief Storage for the implicit a factors
          */
         std::array<MHDFloat,3> mAIm;

         /**
          * @brief Storage for the step fractions
          */
         std::array<MHDFloat,3> mCEx;

      private:

   };

   /// Typedef for a shared pointer ImExPC2
   typedef std::shared_ptr<ImExPC2> SharedImExPC2;

}
}

#endif // QUICC_TIMESTEP_IMEXPC2_HPP
