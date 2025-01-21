/**
 * @file ITripleHarmonicOperator.hpp
 * @brief Implementation of the generic spectral triple harmonic operator
 */

#ifndef QUICC_DENSESM_CHEBYSHEV_LINEARMAP_ITRIPLEHARMONICOPERATOR_HPP
#define QUICC_DENSESM_CHEBYSHEV_LINEARMAP_ITRIPLEHARMONICOPERATOR_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/ILinearMapOperator.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

/**
 * @brief Implementation of the generic spectral triple harmonic operator
 */
class ITripleHarmonicOperator : public ILinearMapOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param rows    Number of rows
    * @param cols    Number of cols
    * @param p       Power of radial prefactor
    * @param lOut    Output harmonic degree
    * @param mOut    Output harmonic degree
    * @param lF      harmonic degree of f
    * @param mF      harmonic degree of f
    * @param lIn     Input harmonic degree
    * @param mIn     Input harmonic degree
    * @param pF      Radial function pointer
    * @param lower   Lower boundary
    * @param upper   Upper boundary
    */
   ITripleHarmonicOperator(const int rows, const int cols, const int p, const int lOut,
      const int mOut, const int lF, const int mF, const int lIn, const int mIn,
      std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower,
      const Scalar_t upper) :
       ILinearMapOperator(rows, cols, lower, upper),
       mP(p),
       mLout(lOut),
       mMout(mOut),
       mLf(lF),
       mMf(mF),
       mLin(lIn),
       mMin(mIn),
       mpF(pF) {};

   /**
    * @brief Destructor
    */
   virtual ~ITripleHarmonicOperator() = default;

protected:
   /**
    * @brief Power of radial prefactor
    */
   const int mP;

   /**
    * @brief Harmonic degree of output
    */
   int mLout;

   /**
    * @brief Harmonic order of output
    */
   int mMout;

   /**
    * @brief Harmonic degree of f
    */
   int mLf;

   /**
    * @brief Harmonic order of f
    */
   int mMf;

   /**
    * @brief Harmonic degree of input
    */
   int mLin;

   /**
    * @brief Harmonic order of input
    */
   int mMin;

   /**
    * @brief Functor for f
    */
   std::shared_ptr<RadialTorPolFunction> mpF;

private:
};

   /**
    * @brief Compute operator for product with spectral expansion
    *
    * @param mat   Matrix to store operator
    * @param rows  Rows of operator
    * @param cols  Cols of operator
    * @param spec  Spectral expansion coefficients
    * @param spec  Size of expansion
    */
template <typename TMat, typename TSpec> void expansionProduct(TMat& mat, const int rows, const int cols, const TSpec& spec, const int nN)
{
   assert(mat.rows() >= rows);
   assert(mat.cols() >= cols);
   assert(spec.rows() >= nN);
   assert(spec.cols() >= 1);

   for(int i = 0; i < rows; i++)
   {
      for(int j = 0; j < nN; j++)
      {
         if(i-j >= 0)
         {
            mat(i,i-j) += spec(j);
         }
         else if(j-i < nN)
         {
            mat(i,j-i) += spec(j);
         }
         if(j > 0 && i+j < cols)
         {
            mat(i,i+j) += spec(j);
         }
      }
   }
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_ITRIPLEHARMONICOPERATOR_HPP
