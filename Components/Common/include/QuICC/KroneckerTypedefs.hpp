/**
 * @file KroneckerTypedefs.hpp
 * @brief Some general typedefs used for kronecker products
 */

#ifndef QUICC_KRONECKERTYPEDEFS_HPP
#define QUICC_KRONECKERTYPEDEFS_HPP

// Configuration includes
//

// System includes
//
#include <tuple>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"

namespace QuICC {

   /**
    * @name Kronecker product of real operators
    */
   //@{
   // No eigen dimension
   typedef std::tuple<SparseMatrix,SparseMatrix,SparseMatrix> KronNoEigenRProduct;
   // Single eigen dimension
   typedef std::tuple<SparseMatrix,SparseMatrix> KronEigen1DRProduct;
   // Two eigen dimension
   typedef SparseMatrix KronEigen2DRProduct;
   // Three eigen dimension
   typedef MHDFloat KronEigen3DRProduct;
   //@}

   /**
    * @name Sum of Kronecker products of real operators
    */
   //@{
   // No eigen dimension
   typedef std::vector<KronNoEigenRProduct> KronNoEigenRSum;
   // Single eigen dimension
   typedef std::vector<KronEigen1DRProduct> KronEigen1DRSum;
   // Two eigen dimension
   typedef SparseMatrix KronEigen2DRSum;
   // Three eigen dimension
   typedef MHDFloat KronEigen3DRSum;
   //@}

   /**
    * @name Kronecker product of complex operators in decoupled form
    */
   //@{
   // No eigen dimension
   typedef std::tuple<DecoupledZSparse,DecoupledZSparse,DecoupledZSparse> KronNoEigenZProduct;
   // Single eigen dimension
   typedef std::tuple<DecoupledZSparse,DecoupledZSparse> KronEigen1DZProduct;
   // Two eigen dimension
   typedef DecoupledZSparse KronEigen2DZProduct;
   // Three eigen dimension
   typedef MHDComplex KronEigen3DZProduct;
   //@}

   /**
    * @name Sum of Kronecker products of complex operators in decoupled form
    */
   //@{
   // No eigen dimension
   typedef std::vector<KronNoEigenZProduct> KronNoEigenZSum;
   // Single eigen dimension
   typedef std::vector<KronEigen1DZProduct> KronEigen1DZSum;
   // Two eigen dimension
   typedef DecoupledZSparse KronEigen2DZSum;
   // Three eigen dimension
   typedef MHDComplex KronEigen3DZSum;
   //@}
}

#endif // QUICC_KRONECKERTYPEDEFS_HPP
