/**
 * @file ThreeTermRecurrence.hpp
 * @brief General implementation for three term recurrence relations
 */

#ifndef QUICC_POLYNOMIAL_THREETERMRECURRENCE_HPP
#define QUICC_POLYNOMIAL_THREETERMRECURRENCE_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//
#include <functional>

// External includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace Polynomial {

   /**
    * @brief Implementation of the three term recurrence relations
    */
   class ThreeTermRecurrence
   {
      public:
         /// Typedef for the function signature of an n independent constant normalizer
         typedef std::function<Internal::Array()> NormalizerC;

         /// Typedef for the function signature of an n dependent constant normalizer
         typedef std::function<Internal::Array(const Internal::MHDFloat)> NormalizerNC;

         /// Typedef for the function signature of an n independent one parameter normalizer
         typedef std::function<Internal::Array(const Internal::MHDFloat)> NormalizerA;

         /// Typedef for the function signature of an n dependent one parameter normalizer
         typedef std::function<Internal::Array(const Internal::MHDFloat, const Internal::MHDFloat)> NormalizerNA;

         /// Typedef for the function signature of an n independent two parameter normalizer
         typedef std::function<Internal::Array(const Internal::MHDFloat, const Internal::MHDFloat)> NormalizerAB;

         /// Typedef for the function signature of an n dependent two parameter normalizer
         typedef std::function<Internal::Array(const Internal::MHDFloat, const Internal::MHDFloat, const Internal::MHDFloat)> NormalizerNAB;

         /**
          * @brief Compute three term recurrence for $P_n(x)$ normalizer
          *
          * Internal computation can be done in multiple precision
          */
         static void Pn(Eigen::Ref<Internal::Matrix> ipn, const int n, const Eigen::Ref<const Internal::Matrix>& ipn_1, const Eigen::Ref<const Internal::Matrix>& ipn_2, const Internal::Array& igrid, NormalizerNC norm);

         /**
          * @brief Compute three term recurrence for $P_1(x)$ normalizer
          *
          * Internal computation can be done in multiple precision
          */
         static void P1(Eigen::Ref<Internal::Matrix> ip1, const Eigen::Ref<const Internal::Matrix>& ip0, const Internal::Array& igrid, NormalizerC norm);

         /**
          * @brief Compute three term recurrence for $P_0(x)$ normalizer
          *
          * Internal computation can be done in multiple precision
          */
         static void P0(Eigen::Ref<Internal::Matrix> ip0, const Internal::Array& igrid, NormalizerC norm);

         /**
          * @brief Compute three term recurrence for $P_n(x)^{(\alpha)}$ with one parameter normalizer
          *
          * Internal computation can be done in multiple precision
          */
         static void Pn(Eigen::Ref<Internal::Matrix> ipn, const int n, const Internal::MHDFloat alpha, const Eigen::Ref<const Internal::Matrix>& ipn_1, const Eigen::Ref<const Internal::Matrix>& ipn_2, const Internal::Array& igrid, NormalizerNA norm);

         /**
          * @brief Compute three term recurrence for $P_1(x)^{(\alpha)}$ with one parameter normalizer
          *
          * Internal computation can be done in multiple precision
          */
         static void P1(Eigen::Ref<Internal::Matrix> ip1, const Internal::MHDFloat alpha, const Eigen::Ref<const Internal::Matrix>& ip0, const Internal::Array& igrid, NormalizerA norm);

         /**
          * @brief Compute three term recurrence for $P_0(x)^{(\alpha)}$ with one parameter normalizer
          *
          * Internal computation can be done in multiple precision
          */
         static void P0(Eigen::Ref<Internal::Matrix> ip0, const Internal::MHDFloat alpha, const Internal::Array& igrid, NormalizerA norm);

         /**
          * @brief Compute three term recurrence for $P_n(x)^{(\alpha,\beta)}$ with two parameter normalizer
          *
          * Internal computation can be done in multiple precision
          */
         static void Pn(Eigen::Ref<Internal::Matrix> ipn, const int n, const Internal::MHDFloat alpha, const Internal::MHDFloat beta, const Eigen::Ref<const Internal::Matrix>& ipn_1, const Eigen::Ref<const Internal::Matrix>& ipn_2, const Internal::Array& igrid, NormalizerNAB norm);

         /**
          * @brief Compute three term recurrence for $dP_n(x)^{(\alpha,\beta)}$ with two parameter normalizer
          *
          * Internal computation can be done in multiple precision
          */
         static void dPn(Eigen::Ref<Internal::Matrix> idpn, const int n, const Internal::MHDFloat alpha, const Internal::MHDFloat beta, const Eigen::Ref<const Internal::Matrix>& ipn, const Eigen::Ref<const Internal::Matrix>& ipn_1, const Internal::Array& igrid, NormalizerNAB norm);

         /**
          * @brief Compute \f$P_1^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void P1(Eigen::Ref<Internal::Matrix> ip1, const Internal::MHDFloat alpha, const Internal::MHDFloat beta, const Eigen::Ref<const Internal::Matrix>& ip0, const Internal::Array& igrid, NormalizerAB norm);

          /**
          * @brief Compute \f$P_1^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void P1(Eigen::Ref<Internal::Matrix> ip1, const Internal::MHDFloat alpha, const Internal::MHDFloat beta, const Internal::Array& igrid, NormalizerAB norm);


         /**
          * @brief Compute \f$P_0^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void P0(Eigen::Ref<Internal::Matrix> ip0, const Internal::MHDFloat alpha, const Internal::MHDFloat beta, const Internal::Array& igrid, NormalizerAB norm);

      private:
         /**
          * @brief Constructor
          */
         ThreeTermRecurrence();

         /**
          * @brief Destructor
          */
         ~ThreeTermRecurrence();

   };
}
}

#endif // QUICC_POLYNOMIAL_THREETERMRECURRENCE_HPP
