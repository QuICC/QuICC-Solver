/**
 * @file Precision.hpp
 * @brief Small wrapper class for generic normal of multiple precision internal computations
 */

#ifndef QUICC_PRECISION_HPP
#define QUICC_PRECISION_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"

#ifdef QUICC_MULTPRECISION
   #include "Types/MpTypedefs.hpp"

   /// Define a small macro to replace float constants to strings in the case of MP computations
   #define MHD_MP(c) MHDMpFloat(#c)

   /// Define a small macro to replace float constants to strings in the case of MP computations
   #define MHD_MP_LONG(c) MHDMpFloat(#c)
#else
   #include <cmath>
   /// For normal computations the macro does nothing
   #define MHD_MP(c) c
   /// For normal computations the macro does nothing
   #define MHD_MP_LONG(c) c ## L
#endif // QUICC_MULTPRECISION

namespace QuICC {

   namespace internal {

      #ifdef QUICC_MULTPRECISION
      //////////////////////////////////////////////////////////////////
         /// Typedef for the internal float type
         typedef QuICC::MHDMpFloat MHDFloat;

         /// Typedef for the internal long double type
         typedef QuICC::MHDMpFloat MHDLong;

         /// Typedef for the internal Array type
         typedef QuICC::MpArray Array;

         /// Typedef for the internal coefficient array type
         typedef QuICC::MpACoeff ACoeff;

         /// Typedef for the internal Matrix type
         typedef QuICC::MpMatrix Matrix;

         /// Typedef for the internal Matrix type
         typedef QuICC::MpSparseMatrix SparseMatrix;

         /// Typedef for the internal Array type
         typedef Array   ArrayL;

         /// Typedef for the internal Matrix type
         typedef Matrix   MatrixL;

         /// Typedef for the smart internal Array type
         typedef QuICC::SharedMpArray SharedArray;

         /// Typedef for the smart internal Matrix type
         typedef QuICC::SharedMpMatrix SharedMatrix;
      //////////////////////////////////////////////////////////////////
      #else
      //////////////////////////////////////////////////////////////////
         /// Typedef for the internal float type
         typedef QuICC::MHDFloat MHDFloat;

         /// Typedef for the internal long double type
         typedef QuICC::MHDLong MHDLong;

         /// Typedef for the internal Array type
         typedef QuICC::Array Array;

         /// Typedef for the internal Matrix type
         typedef QuICC::ACoeff ACoeff;

         /// Typedef for the internal Matrix type
         typedef QuICC::Matrix Matrix;

         /// Typedef for the internal Matrix type
         typedef QuICC::SparseMatrix SparseMatrix;

         /// Typedef for the internal Array type
         typedef Eigen::Matrix<MHDLong, Eigen::Dynamic, 1>   ArrayL;

         /// Typedef for the internal Matrix type
         typedef Eigen::Matrix<MHDLong, Eigen::Dynamic, Eigen::Dynamic>   MatrixL;

         /// Typedef for the smart internal Array type
         typedef QuICC::SharedArray SharedArray;

         /// Typedef for the smart internal Matrix type
         typedef QuICC::SharedMatrix SharedMatrix;
      //////////////////////////////////////////////////////////////////
      #endif // QUICC_MULTPRECISION
   }

   /**
    * @brief Simple class holding some typedefs to allow for internal MP computations
    */
   class Precision
   {
      public:

         /**
          * @brief Precision dependent mathematical constant \f$\pi\f$
          */
         static const internal::MHDFloat PI;

         /**
          * @brief Long double or MP precision mathematical constant \f$\pi\f$
          */
         static const internal::MHDLong PI_long;

         /**
          * @brief Flag to test whether Multiple Precision is enabled
          */
         static const bool hasMP;

         /**
          * @brief Initialise the precision setup
          */
         static void init();

         /**
          * @brief Cast the internal value to an external one
          *
          * @param val Internal value to cast
          */
         static MHDFloat cast(internal::MHDFloat val);

         /**
          * @brief Cast the internal smart Array to an external one
          *
          * @param spIArr Internal smart Array to cast
          */
         static SharedArray cast(internal::SharedArray spIArr);

         /**
          * @brief Cast the internal smart Matrix to an external one
          *
          * @param spIMat Internal smart Matrix to cast
          */
         static SharedMatrix cast(internal::SharedMatrix spIMat);

         /**
          * @brief Cast the internal smart Array to an external one
          *
          * @param rIArr Internal Array to cast
          */
         static Array cast(const internal::Array& rIArr);

         /**
          * @brief Cast the internal smart Matrix to an external one
          *
          * @param rIMat Internal Matrix to cast
          */
         static Matrix cast(const internal::Matrix& rIMat);

         /**
          * @brief Cast the internal Eigen::Ref Matrix to an external one
          *
          * @param rIMat Internal Matrix to cast
          */
         static Matrix cast(const Eigen::Ref<const internal::Matrix>& rIMat);

      private:
         /**
         * @brief Empty constructor
         */
         Precision();

         /**
         * @brief Simple empty destructor
         */
         ~Precision();
   };

   inline void Precision::init()
   {
   }

   inline MHDFloat Precision::cast(internal::MHDFloat val)
   {
      #ifdef QUICC_MULTPRECISION
         return val.convert_to<MHDFloat>();
      #else
         // flush to zero subnormal values
         if(std::abs(val) < std::numeric_limits<MHDFloat>::min())
         {
            val = 0.0;
         }
         return val;
      #endif // QUICC_MULTPRECISION
   }

   namespace implementationDetail
   {
   template <class T1, class T2, std::enable_if_t<std::is_same_v<typename T1::Scalar, MHDFloat> &&
      std::is_same_v<typename T2::Scalar, internal::MHDFloat>, bool> = true>
   inline void cast(Eigen::DenseBase<T1>& OMat, const  Eigen::DenseBase<T2>& IMat)
   {
      // Loop over whole matrix
      for(int j=0; j < IMat.cols(); ++j)
      {
         for(int i=0; i < IMat.rows(); ++i)
         {
            OMat(i,j) = Precision::cast(IMat(i,j));
         }
      }
   }
   }

   inline SharedArray Precision::cast(internal::SharedArray spIArr)
   {
      auto spArr = std::make_shared<Array>(spIArr->size());

      implementationDetail::cast(*spArr, *spIArr);

      return spArr;
   }

   inline SharedMatrix Precision::cast(internal::SharedMatrix spIMat)
   {
      auto spMat = std::make_shared<Matrix>(spIMat->rows(),spIMat->cols());

      implementationDetail::cast(*spMat, *spIMat);

      return spMat;
   }

   inline Array Precision::cast(const internal::Array& rIArr)
   {
      Array arr(rIArr.size());

      implementationDetail::cast(arr, rIArr);

      return arr;
   }

   inline Matrix Precision::cast(const internal::Matrix& rIMat)
   {
      Matrix mat(rIMat.rows(),rIMat.cols());

      implementationDetail::cast(mat, rIMat);

      return mat;
   }

   inline Matrix Precision::cast(const Eigen::Ref<const internal::Matrix>& rIMat)
   {
      Matrix mat(rIMat.rows(),rIMat.cols());

      implementationDetail::cast(mat, rIMat);

      return mat;
   }

#ifdef QUICC_MULTPRECISION
   /// Create a namespace alias for the internal precision stuff pointing to mpfr namespace
   namespace  precision = boost::multiprecision;
#else
   /// Create a namespace alias for the internal precision stuff pointing to std namespace
   namespace  precision = std;
#endif // QUICC_MULTPRECISION
}

#endif // QUICC_PRECISION_HPP
