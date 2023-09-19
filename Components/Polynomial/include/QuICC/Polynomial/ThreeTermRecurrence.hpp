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
#include "Types/Precision.hpp"

namespace QuICC {

namespace Polynomial {

   /**
    * @brief Implementation of the three term recurrence relations
    */
   class ThreeTermRecurrence
   {
      public:
         /// Typedef for the function signature of an n independent constant normalizer
         typedef std::function<internal::Array()> NormalizerC;

         /// Typedef for the function signature of an n dependent constant normalizer
         typedef std::function<internal::Array(const internal::MHDFloat)> NormalizerNC;

         /// Typedef for the function signature of an n independent one parameter normalizer
         typedef std::function<internal::Array(const internal::MHDFloat)> NormalizerA;

         /// Typedef for the function signature of an n dependent one parameter normalizer
         typedef std::function<internal::Array(const internal::MHDFloat, const internal::MHDFloat)> NormalizerNA;

         /// Typedef for the function signature of an n independent two parameter normalizer
         typedef std::function<internal::Array(const internal::MHDFloat, const internal::MHDFloat)> NormalizerAB;

         /// Typedef for the function signature of an n dependent two parameter normalizer
         typedef std::function<internal::Array(const internal::MHDFloat, const internal::MHDFloat, const internal::MHDFloat)> NormalizerNAB;

         /**
          * @brief Compute three term recurrence for $P_n(x)$ normalizer
          *
          * Internal computation can be done in multiple precision
          */
         static void Pn(Eigen::Ref<internal::Matrix> ipn, const int n, const Eigen::Ref<const internal::Matrix>& ipn_1, const Eigen::Ref<const internal::Matrix>& ipn_2, const internal::Array& igrid, NormalizerNC norm);

         /**
          * @brief Compute three term recurrence for $P_1(x)$ normalizer
          *
          * Internal computation can be done in multiple precision
          */
         static void P1(Eigen::Ref<internal::Matrix> ip1, const Eigen::Ref<const internal::Matrix>& ip0, const internal::Array& igrid, NormalizerC norm);

         /**
          * @brief Compute three term recurrence for $P_0(x)$ normalizer
          *
          * Internal computation can be done in multiple precision
          */
         static void P0(Eigen::Ref<internal::Matrix> ip0, const internal::Array& igrid, NormalizerC norm);

         /**
          * @brief Compute three term recurrence for $P_n(x)^{(\alpha)}$ with one parameter normalizer
          *
          * Internal computation can be done in multiple precision
          */
         static void Pn(Eigen::Ref<internal::Matrix> ipn, const int n, const internal::MHDFloat alpha, const Eigen::Ref<const internal::Matrix>& ipn_1, const Eigen::Ref<const internal::Matrix>& ipn_2, const internal::Array& igrid, NormalizerNA norm);

         /**
          * @brief Compute three term recurrence for $P_1(x)^{(\alpha)}$ with one parameter normalizer
          *
          * Internal computation can be done in multiple precision
          */
         static void P1(Eigen::Ref<internal::Matrix> ip1, const internal::MHDFloat alpha, const Eigen::Ref<const internal::Matrix>& ip0, const internal::Array& igrid, NormalizerA norm);

         /**
          * @brief Compute three term recurrence for $P_0(x)^{(\alpha)}$ with one parameter normalizer
          *
          * Internal computation can be done in multiple precision
          */
         static void P0(Eigen::Ref<internal::Matrix> ip0, const internal::MHDFloat alpha, const internal::Array& igrid, NormalizerA norm);

         /**
          * @brief Compute three term recurrence for $P_n(x)^{(\alpha,\beta)}$ with two parameter normalizer
          *
          * Internal computation can be done in multiple precision
          */
         static void Pn(Eigen::Ref<internal::Matrix> ipn, const int n, const internal::MHDFloat alpha, const internal::MHDFloat beta, const Eigen::Ref<const internal::Matrix>& ipn_1, const Eigen::Ref<const internal::Matrix>& ipn_2, const internal::Array& igrid, NormalizerNAB norm);

         /**
          * @brief Compute three term recurrence for $dP_n(x)^{(\alpha,\beta)}$ with two parameter normalizer
          *
          * Internal computation can be done in multiple precision
          */
         static void dPn(Eigen::Ref<internal::Matrix> idpn, const int n, const internal::MHDFloat alpha, const internal::MHDFloat beta, const Eigen::Ref<const internal::Matrix>& ipn, const Eigen::Ref<const internal::Matrix>& ipn_1, const internal::Array& igrid, NormalizerNAB norm);

         /**
          * @brief Compute \f$P_1^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void P1(Eigen::Ref<internal::Matrix> ip1, const internal::MHDFloat alpha, const internal::MHDFloat beta, const Eigen::Ref<const internal::Matrix>& ip0, const internal::Array& igrid, NormalizerAB norm);

          /**
          * @brief Compute \f$P_1^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void P1(Eigen::Ref<internal::Matrix> ip1, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid, NormalizerAB norm);


         /**
          * @brief Compute \f$P_0^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void P0(Eigen::Ref<internal::Matrix> ip0, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid, NormalizerAB norm);

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
