/**
 * @file ImExRKCB4.hpp
 * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order 4 (Cavaglieri & Bewley, 2015)
 */

#ifndef QUICC_TIMESTEP_IMEXRKCB4_HPP
#define QUICC_TIMESTEP_IMEXRKCB4_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Timestep/IImExRKCBScheme.hpp"
#include "QuICC/Timestep/SparseImExRK3RTimestepper.hpp"

namespace QuICC {

namespace Timestep {

   /**
    * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order 4
    */
   class ImExRKCB4: public IImExRKCBScheme
   {
      public:
         /// Typedef for underlying scheme implementation
         template <typename T1,typename T2,template <typename> class T3>
            using ImplementationType = SparseImExRK3RTimestepper<T1,T2,T3>;

         /**
          * @brief Constructor
          */
         ImExRKCB4();

         /**
          * @brief Destructor
          */
         virtual ~ImExRKCB4() = default;

         /**
          * @brief Butcher's tableau a_ij factor for implicit scheme
          */
         MHDFloat aIm(const int i, const int j) const;

         /**
          * @brief Butcher's tableau b_i factor for implicit scheme
          */
         MHDFloat bIm(const int i) const;

         /**
          * @brief Butcher's tableau a_ij factor for explicit scheme
          */
         MHDFloat aEx(const int i, const int j) const;

         /**
          * @brief Butcher's tableau b_i factor for explicit scheme
          */
         MHDFloat bEx(const int i) const;

         /**
          * @brief Butcher's tableau c_i factor for explicit scheme
          */
         MHDFloat cEx(const int i) const;

         /**
          * @brief Butcher's tableau b_i factor for implicit embedded lower order scheme
          */
         MHDFloat bImErr(const int i) const;

         /**
          * @brief Butcher's tableau b_i factor for explicit embedded lower order scheme
          */
         MHDFloat bExErr(const int i) const;

         /**
          * @brief Number of substeps for final step (this is +1 compared to theoretical value due to implementation)
          */
         int steps() const;

         /**
          * @brief Order of the scheme
          */
         int order() const;

         /**
          * @brief Scheme has embedded lower order scheme?
          */
         bool hasEmbedded() const;

         /**
          * @brief Name of the scheme
          */
         std::string name() const;

      protected:
         /**
          * @brief Initialize Butcher's tableau
          */
         void init();

         /**
          * @brief Storage for the implicit a factors
          */
         Eigen::Array<MHDFloat,6,6> mAIm;

         /**
          * @brief Storage for the implicit b factors
          */
         Eigen::Array<MHDFloat,6,1> mBIm;

         /**
          * @brief Storage for the explicit a factors
          */
         Eigen::Array<MHDFloat,6,6> mAEx;

         /**
          * @brief Storage for the explicit b factors
          */
         Eigen::Array<MHDFloat,6,1> mBEx;

         /**
          * @brief Storage for the explicit c factors
          */
         Eigen::Array<MHDFloat,6,1> mCEx;

         /**
          * @brief Storage for the implicit embedded scheme b factors
          */
         Eigen::Array<MHDFloat,6,1> mBImErr;

         /**
          * @brief Storage for the explicit embedded scheme b factors
          */
         Eigen::Array<MHDFloat,6,1> mBExErr;

      private:

   };

   /// Typedef for a shared pointer ImExRKCB4
   typedef std::shared_ptr<ImExRKCB4> SharedImExRKCB4;
}
}

#endif // QUICC_TIMESTEP_IMEXRKCB4_HPP
