/**
 * @file MpTypedefs.hpp
 * @brief Definition of typedefs for muliple precision computations
 */

#ifndef QUICC_MPTYPEDEFS_HPP
#define QUICC_MPTYPEDEFS_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//
#if defined QUICC_MPBACKEND_BOOST
#include <boost/multiprecision/cpp_dec_float.hpp>
#elif defined QUICC_MPBACKEND_GMP
#include <boost/multiprecision/gmp.hpp>
#elif defined QUICC_MPBACKEND_MPFR
#include <boost/multiprecision/mpfr.hpp>
#elif defined QUICC_MPBACKEND_QUAD
#include <boost/multiprecision/float128.hpp>
#endif // defined QUICC_MPBACKEND_BOOST
#include <boost/multiprecision/eigen.hpp>
#include <Eigen/Core>

// Project includes
//

namespace QuICC {

   /**
    * @name Basic scalar types typedefs
    */
   //@{
   /// Typedef for multiple precision floating point type value
#if defined QUICC_MPBACKEND_BOOST
   typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<QUICC_MULTPRECISION_DIGITS> > MHDMpFloat;
#elif defined QUICC_MPBACKEND_GMP
   typedef boost::multiprecision::number<boost::multiprecision::gmp_float<QUICC_MULTPRECISION_DIGITS> > MHDMpFloat;
#elif defined QUICC_MPBACKEND_MPFR
   typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<QUICC_MULTPRECISION_DIGITS> > MHDMpFloat;
#elif defined QUICC_MPBACKEND_QUAD
   typedef boost::multiprecision::float128 MHDMpFloat;
#endif // defined QUICC_MPBACKEND_BOOST
   //@}

   /**
    * @name Array types typedefs
    */
   //@{
   /// Typedef for an array of multiple precision float values
   typedef Eigen::Matrix<MHDMpFloat, Eigen::Dynamic, 1>   MpArray;
   //@}

   /**
    * @name Coefficient wise array types typedefs
    */
   //@{
   /// Typedef for an array of float values
   typedef Eigen::Array<MHDMpFloat, Eigen::Dynamic, 1>   MpACoeff;
   //@}

   /**
    * @name Matrix types typedefs
    */
   //@{
   /// Typedef for a matrix of multiple precision float values
   typedef Eigen::Matrix<MHDMpFloat, Eigen::Dynamic, Eigen::Dynamic>  MpMatrix;
   //@}
   //
   /**
    * @name Sparse Matrix types typedefs
    */
   //@{
   /// Typedef for a sparse matrix of float values
   typedef Eigen::SparseMatrix<MHDMpFloat>   MpSparseMatrix;
   //@}

   /**
    * @name Shared pointer typedefs
    */
   //@{
   /// Typedef for an smart reference counting pointer on an array of multiple precision real values
   typedef std::shared_ptr<MpArray>   SharedMpArray;
   /// Typedef for an smart reference counting pointer on an matrix of real values
   typedef std::shared_ptr<MpMatrix>   SharedMpMatrix;
   //@}

}

#endif // QUICC_MPTYPEDEFS_HPP
