/**
 * @file PardisoLU.hpp
 * @brief Declarations needed to use the Pardiso routines in the C++ 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PARDISOLU_HPP
#define PARDISOLU_HPP

#include "Eigen/SparseCore"
#include "Eigen/src/Core/util/DisableStupidWarnings.h"

#include "Eigen/src/misc/Solve.h"
#include "Eigen/src/misc/SparseSolve.h"

#include "../External/Interfaces/Pardiso_Interface.hpp"
#include <complex>

namespace Eigen { 

namespace Pardiso { 
   /**
    * @brief Declaration for the real valued pardisoinit routine
    */
   void pardisoinit(void *, int *, int *, double *, int *, double);

   /**
    * @brief Declaration for the complex valued pardisoinit routine
    */
   void pardisoinit(void *, int *, int *, double *, int *, std::complex<double>);

   /**
    * @brief Declaration for the complex valued pardiso routine
    */
   void pardiso(void *, const int *, const int *, const int *, const int *, const std::complex<double> *, const int *, const int *, const int *, const int *, int *, const int *, std::complex<double> *, std::complex<double> *, int *, double *, std::complex<double>);

   /**
    * @brief Declaration for the real valued pardiso routine
    */
   void pardiso(void *, const int *, const int *, const int *, const int *, const double *, const int *, const int *, const int *, const int *, int *, const int *, double *, double *, int *, double *, double);

   /**
    * @brief Declaration for the real valued pardiso_chkmatrix routine
    */
   inline void pardiso_chkmatrix(int *mtype, int *n, double *a, int *ia, int *ja, int *error, double)
   {
      ::pardiso_chkmatrix(mtype,n,a,ia,ja,error);
   }

   /**
    * @brief Declaration for the complex valued pardiso_chkmatrix routine
    */
   inline void pardiso_chkmatrix(int *mtype, int *n, pardiso_complex *a, int *ia, int *ja, int *error, std::complex<double>)
   {
      ::pardiso_chkmatrix_z(mtype,n,a,ia,ja,error);
   }

   /**
    * @brief Declaration for the real valued pardiso_chkvec routine
    */
   inline void pardiso_chkvec(int *n, int *nrhs, double *b, int *error, double)
   {
      ::pardiso_chkvec(n,nrhs,b,error);
   }

   /**
    * @brief Declaration for the complex valued pardiso_chkvec routine
    */
   inline void pardiso_chkvec(int *n, int *nrhs, pardiso_complex *b, int *error, std::complex<double>)
   {
      ::pardiso_chkvec_z(n,nrhs,b,error);
   }

   /**
    * @brief Declaration for the real valued pardiso_printstats routine
    */
   inline void pardiso_printstats(int *mtype, int *n, double *a, int *ia, int *ja, int *nrhs, double *b, int *error, double)
   {
      ::pardiso_printstats(mtype,n,a,ia,ja,nrhs,b,error);
   }

   /**
    * @brief Declaration for the complex valued pardiso_printstats routine
    */
   inline void pardiso_printstats(int *mtype, int *n, pardiso_complex *a, int *ia, int *ja, int *nrhs, pardiso_complex *b, int *error, std::complex<double>)
   {
      ::pardiso_printstats_z(mtype,n,a,ia,ja,nrhs,b,error);
   }

} // end namespace Pardiso

template <typename TMatrixType>
class PardisoLU: internal::noncopyable
{
  public:
    typedef TMatrixType MatrixType;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::RealScalar RealScalar;
    typedef typename MatrixType::Index Index;
    typedef Matrix<Scalar,Dynamic,1> Vector;
    typedef Matrix<int, 1, MatrixType::ColsAtCompileTime> IntRowVectorType;
    typedef Matrix<int, MatrixType::RowsAtCompileTime, 1> IntColVectorType;
    typedef SparseMatrix<Scalar,RowMajor,int> PardisoMatrixType;
    typedef Array<Index,64,1> PardisoIParmType;
    typedef Array<double,64,1> PardisoDParmType;

  public:

    PardisoLU() { init(); }

    PardisoLU(const MatrixType& matrix)
    {
      init();
      compute(matrix);
    }

    ~PardisoLU()
    {
       int phase = -1; // Release internal memory.
       int n = m_copyMatrix.rows();
       int nrhs = 1;
       int idum;
       Scalar ddum;
       int msglvl = 0;

       Pardiso::pardiso(m_solverPtr, &m_maxfct, &m_mnum, &phase,
             &n, &ddum, m_outerIndexPtr, m_innerIndexPtr, &idum, &nrhs,
             m_iparm.data(), &msglvl, &ddum, &ddum, &m_error, m_dparm.data(), Scalar());
    }

    inline Index rows() const { return m_copyMatrix.rows(); }
    inline Index cols() const { return m_copyMatrix.cols(); }

    inline void setIterativeRefinement(const int iterations) { m_iparm(7) = iterations; }
    inline void setAnalyzeMsgLvl(const int lvl) { m_analyzeMsgLvl = lvl; }
    inline void setFactorMsgLvl(const int lvl) { m_factorMsgLvl = lvl; }
    inline void setSolveMsgLvl(const int lvl) { m_solveMsgLvl = lvl; }
    inline int  getPeakMemory() const { return m_peakmem; }

    /** \brief Reports whether previous computation was successful.
      *
      * \returns \c Success if computation was succesful,
      *          \c NumericalIssue if the matrix.appears to be negative.
      */
    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "Decomposition is not initialized.");
      return m_info;
    }

    /** Computes the sparse LU decomposition of \a matrix 
     *  Note that the matrix should be column-major, and in compressed format for best performance.
     *  \sa SparseMatrix::makeCompressed().
     */
    void compute(const MatrixType& matrix)
    {
      analyzePattern(matrix);
      factorize(matrix);
    }

    /** \returns the solution x of \f$ A x = b \f$ using the current decomposition of A.
      *
      * \sa compute()
      */
    template<typename Rhs>
    inline const internal::solve_retval<PardisoLU, Rhs> solve(const MatrixBase<Rhs>& b) const
    {
      eigen_assert(m_isInitialized && "PardisoLU is not initialized.");
      eigen_assert(rows()==b.rows()
                && "PardisoLU::solve(): invalid number of rows of the right hand side matrix b");
      return internal::solve_retval<PardisoLU, Rhs>(*this, b.derived());
    }

    /** \returns the solution x of \f$ A x = b \f$ using the current decomposition of A.
      *
      * \sa compute()
      */
    template<typename Rhs>
    inline const internal::sparse_solve_retval<PardisoLU, Rhs> solve(const SparseMatrixBase<Rhs>& b) const
    {
      eigen_assert(m_isInitialized && "PardisoLU is not initialized.");
      eigen_assert(rows()==b.rows()
                && "PardisoLU::solve(): invalid number of rows of the right hand side matrix b");
      return internal::sparse_solve_retval<PardisoLU, Rhs>(*this, b.derived());
    }

    /** Performs a symbolic decomposition on the sparcity of \a matrix.
      *
      * This function is particularly useful when solving for several problems having the same structure.
      *
      * \sa factorize(), compute()
      */
    void analyzePattern(const MatrixType& matrix)
    {
      copyInput(matrix);

      int phase = 11;
      int n = m_copyMatrix.rows();
      int nrhs = 1;
      int idum;
      Scalar ddum;

      Pardiso::pardiso (m_solverPtr, &m_maxfct, &m_mnum, &phase,
            &n, m_valuePtr, m_outerIndexPtr, m_innerIndexPtr, &idum, &nrhs,
            m_iparm.data(), &m_analyzeMsgLvl, &ddum, &ddum, &m_error,  m_dparm.data(), Scalar());

      // Extract peak memory usage
      m_peakmem = std::max(m_iparm(14), m_iparm(15) + m_iparm(16));

      m_isInitialized = true;
      m_info = m_error ? InvalidInput : Success;
      m_analysisIsOk = m_error ? false : true;
      m_factorizationIsOk = false;
    }

    /** Performs a numeric decomposition of \a matrix
      *
      * The given matrix must have the same sparcity than the matrix on which the pattern analysis has been performed.
      *
      * \sa analyzePattern(), compute()
      */
    void factorize(const MatrixType& matrix)
    {
      eigen_assert(m_analysisIsOk && "PardisoLU: you must first call analyzePattern()");

      copyInput(matrix);

      int phase = 22;
      int n = m_copyMatrix.rows();
      int nrhs = 1;
      int idum;
      Scalar ddum;

      Pardiso::pardiso(m_solverPtr, &m_maxfct, &m_mnum, &phase,
            &n, m_valuePtr, m_outerIndexPtr, m_innerIndexPtr, &idum, &nrhs,
            m_iparm.data(), &m_factorMsgLvl, &ddum, &ddum, &m_error, m_dparm.data(), Scalar());

      m_info = m_error ? NumericalIssue : Success;
      m_factorizationIsOk = m_error ? false : true;
    }

    #ifndef EIGEN_PARSED_BY_DOXYGEN
    /** \internal */
    template<typename BDerived,typename XDerived>
    bool _solve(const MatrixBase<BDerived> &b, MatrixBase<XDerived> &x) const;
    #endif

  protected:

    void init()
    {
      m_info = InvalidInput;
      m_isInitialized = false;
      m_outerIndexPtr = 0;
      m_innerIndexPtr = 0;
      m_valuePtr      = 0;
      m_error         = 0;
      m_maxfct        = 1;
      m_mnum          = 1;
      m_analyzeMsgLvl = 0;
      m_factorMsgLvl  = 0;
      m_solveMsgLvl   = 0;
      m_peakmem       = 0;

      /* -------------------------------------------------------------------- */
      /* ..  Setup Pardiso control parameters and initialize the solvers      */
      /*     internal adress pointers. This is only necessary for the FIRST   */
      /*     call of the PARDISO solver.                                      */
      /* ---------------------------------------------------------------------*/
      int solver = 0; // use sparse direct solver

      Pardiso::pardisoinit(m_solverPtr, &solver, m_iparm.data(), m_dparm.data(), &m_error, Scalar());
      if(m_error != 0)
      {
         m_info = InvalidInput;
         eigen_assert(m_error != -10 && "[Pardiso]: No license file found");
         eigen_assert(m_error != -11 && "[Pardiso]: License is expired");
         eigen_assert(m_error != -12 && "[Pardiso]: Wrong username or hostname");
      }

      m_iparm(2) = 1; // Number of processors
      m_iparm(7) = 1;
    }
    
    void copyInput(const MatrixType& mat)
    {
       // non supported input -> copy
       m_copyMatrix = mat;
       m_outerIndexPtr = m_copyMatrix.outerIndexPtr();
       m_innerIndexPtr = m_copyMatrix.innerIndexPtr();
       m_valuePtr      = m_copyMatrix.valuePtr();

      /* -------------------------------------------------------------------- */    
      /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
      /*     notation.                                                        */
      /* -------------------------------------------------------------------- */ 
      for (int i = 0; i < mat.rows()+1; i++)
      {
         m_outerIndexPtr[i] += 1;
      }
      for (int i = 0; i < mat.nonZeros(); i++)
      {
         m_innerIndexPtr[i] += 1;
      }
    }

    PardisoMatrixType m_copyMatrix;
    const Scalar* m_valuePtr;
    int* m_outerIndexPtr;
    int* m_innerIndexPtr;
    mutable void* m_solverPtr[64];
    mutable PardisoIParmType m_iparm;
    mutable PardisoDParmType m_dparm;

    mutable ComputationInfo m_info;
    bool m_isInitialized;
    int m_factorizationIsOk;
    int m_analysisIsOk;

    mutable int m_error;
    int m_maxfct;
    int m_mnum;
    int m_analyzeMsgLvl;
    int m_factorMsgLvl;
    int m_solveMsgLvl;
    int m_peakmem;
    
  private:
    PardisoLU(PardisoLU& ) { }

};

template<typename MatrixType>
template<typename BDerived,typename XDerived>
bool PardisoLU<MatrixType>::_solve(const MatrixBase<BDerived> &b, MatrixBase<XDerived> &x) const
{
  int rhsCols = b.cols();
  eigen_assert((BDerived::Flags&RowMajorBit)==0 && "PardisoLU backend does not support non col-major rhs yet");
  eigen_assert((XDerived::Flags&RowMajorBit)==0 && "PardisoLU backend does not support non col-major result yet");
 
  int phase = 33;
  int n = m_copyMatrix.rows();
  int idum;

  if(b.derived().data() == x.derived().data())
  {
     Scalar ddum;
     m_iparm(6) = 1;
     Pardiso::pardiso(m_solverPtr, &m_maxfct, &m_mnum, &phase,
                        &n, m_valuePtr, m_outerIndexPtr, m_innerIndexPtr, &idum, &rhsCols,
                        m_iparm.data(), &m_solveMsgLvl, b.const_cast_derived().data(), &ddum, &m_error, m_dparm.data(), Scalar());
  } else
  {
     Pardiso::pardiso(m_solverPtr, &m_maxfct, &m_mnum, &phase,
                      &n, m_valuePtr, m_outerIndexPtr, m_innerIndexPtr, &idum, &rhsCols,
                      m_iparm.data(), &m_solveMsgLvl, b.const_cast_derived().data(), x.derived().data(), &m_error, m_dparm.data(), Scalar());
  }

  if (m_error!=0)
    return false;

  return true;
}

namespace internal {

template<typename TMatrixType, typename Rhs>
struct solve_retval<PardisoLU<TMatrixType>, Rhs>
  : solve_retval_base<PardisoLU<TMatrixType>, Rhs>
{
  typedef PardisoLU<TMatrixType> Dec;
  EIGEN_MAKE_SOLVE_HELPERS(Dec,Rhs)

  template<typename Dest> void evalTo(Dest& dst) const
  {
    dec()._solve(rhs(),dst);
  }
};

template<typename TMatrixType, typename Rhs>
struct sparse_solve_retval<PardisoLU<TMatrixType>, Rhs>
  : sparse_solve_retval_base<PardisoLU<TMatrixType>, Rhs>
{
  typedef PardisoLU<TMatrixType> Dec;
  EIGEN_MAKE_SPARSE_SOLVE_HELPERS(Dec,Rhs)

  template<typename Dest> void evalTo(Dest& dst) const
  {
    this->defaultEvalTo(dst);
  }
};

} // end namespace internal

} // end namespace Eigen

#include "Eigen/src/Core/util/ReenableStupidWarnings.h"

#endif // PARDISOLU_HPP
