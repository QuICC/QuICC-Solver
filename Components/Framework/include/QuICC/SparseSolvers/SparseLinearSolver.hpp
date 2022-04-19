/**
 * @file SparseLinearSolver.hpp
 * @brief Implementation of a templated (coupled) linear solver structure
 */

#ifndef QUICC_SOLVER_SPARSELINEARSOLVER_HPP
#define QUICC_SOLVER_SPARSELINEARSOLVER_HPP

// Configuration includes
//

// System includes
//
#include <memory>
#include <stdexcept>

// External includes
//

// Project includes
//
#include "QuICC/Math/Constants.hpp"
#include "QuICC/Framework/MpiFramework.hpp"
#include "QuICC/SparseSolvers/SparseSolverBase.hpp"
#include "QuICC/SparseSolvers/SparseLinearSolverTools.hpp"

namespace QuICC {

namespace Solver {

   namespace internal {

      void addOperators(SparseMatrix& mat, const MHDFloat c, const DecoupledZSparse& decMat);

      void addOperators(SparseMatrixZ& mat, const MHDFloat c, const DecoupledZSparse& decMat);

      void addCorrection(DecoupledZMatrix& rVal, const SparseMatrixZ& corr);

      template <typename TData, typename TCorr> void addCorrection(TData& rVal, const TCorr& corr);
   }

   /**
    * @brief Implementation of a templated (coupled) linear solver structure
    */
   template <typename TOperator, typename TData, template <typename> class TSolver> class SparseLinearSolver: public SparseSolverBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          * @param timeId  Solver timing with respect to timestepping
          */
         SparseLinearSolver(const int start, const std::size_t timeId);

         /**
          * @brief Destructor
          */
         virtual ~SparseLinearSolver();

         /**
          * @brief Initialise the solver matrices storage
          *
          * @param n Size of matrices
          */
         virtual void initMatrices(const int n);

         /**
          * @brief Initialise solution after data was copied
          */
         virtual void initSolutions();

         /**
          * @brief Initialise solver
          */
         void initSolver();

         /**
          * @brief Update solver
          */
         void updateSolver();

         /**
          * @brief Set solver RHS data to zero
          */
         void zeroSolver();

         /**
          * @brief Prepare fields for implicit solve
          */
         virtual bool preSolve();

         /**
          * @brief Solve linear systems
          */
         void solve();

         /**
          * @brief Work on fields after implicit solve
          */
         virtual bool postSolve();

         /**
          * @brief Get the number of linear systems in solver
          */
         int nSystem() const;

         /**
          * @brief Set LHS matrix
          *
          * @param idx Index of the matrix
          */
         TOperator& rLHSMatrix(const MHDFloat id, const int idx);

         /**
          * @brief Add RHS and solution data storage
          *
          * @param rows Number of rows of matrix
          * @param cols Number of columns required
          */
         virtual void addStorage(const int rows, const int cols);

         /**
          * @brief Build the scheme operators
          */
         void buildOperators(const int idx, const std::map<std::size_t, DecoupledZSparse>& ops, const int size);

         /**
          * @brief Set RHS data
          *
          * @param idx   Index of the data
          */
         TData& rRHSData(const int idx);

         /**
          * @brief Get solution data
          *
          * @param idx   Index of the data
          */
         const TData& solution(const int idx) const;

         /**
          * @brief Set solution data
          *
          * @param idx   Index of the data
          */
         TData& rSolution(const int idx);

         /**
          * @brief Get inhomogeneous boundary condition
          *
          * @param idx   Index of the condition
          */
         const Eigen::SparseMatrix<typename TData::Scalar>& inhomogeneous(const int idx) const;

         /**
          * @brief Set inhomogeneous boundary condition
          *
          * @param idx   Index of the condition
          */
         Eigen::SparseMatrix<typename TData::Scalar>& rInhomogeneous(const int idx);

         /**
          * @brief Add inhomogeneous boundary condition to RHS data
          */
         void addInhomogeneous();

         /**
          * @brief Computation error
          */
         virtual MHDFloat error() const;

         /**
          * @brief Finished computation?
          */
         virtual bool finished();

      protected:
         /// Typedef for shared solver
         typedef std::shared_ptr< TSolver<TOperator> > SharedSolverType;
         /// Typedef for vector of shared solvers
         typedef std::vector<SharedSolverType, Eigen::aligned_allocator<SharedSolverType> > VectorSharedSolverType;

         /**
          * @brief Initialise the solver matrices storage
          *
          * @param n Size of matrices
          */
         void initMatrices(const MHDFloat id, const int n);

         /**
          * @brief Additional updates after solver update
          */
         virtual void postSolverUpdate();

         /**
          * @brief Correct solution obtained from linear solver
          */
         virtual int correctSolution(const int iteration);

         /**
          * @brief Current ID of the solvers
          */
         MHDFloat mId;

         /**
          * @brief LHS operators
          */
         std::map<MHDFloat, std::vector<TOperator> >  mLHSMatrix;

         /**
          * @brief Storage for linear solve's RHS
          */
         std::vector<TData>  mRHSData;

         /**
          * @brief Storage for solution of linear solve
          */
         std::vector<TData>  mSolution;

         /**
          * @brief Storage for inhomogeneous boundary conditions
          */
         std::vector<Eigen::SparseMatrix<typename TData::Scalar> >  mInhomogeneous;

         /**
          * @brief Create sparse solvers
          */
         std::map<MHDFloat,VectorSharedSolverType>  mSolver;

      private:
   };

   template <typename TOperator,typename TData,template <typename> class TSolver> SparseLinearSolver<TOperator,TData,TSolver>::SparseLinearSolver(const int start, const std::size_t timeId)
      : SparseSolverBase(start, timeId), mId(0.0)
   {
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> SparseLinearSolver<TOperator,TData,TSolver>::~SparseLinearSolver()
   {
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> MHDFloat SparseLinearSolver<TOperator,TData,TSolver>::error() const
   {
      return -1.0;
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> bool SparseLinearSolver<TOperator,TData,TSolver>::finished()
   {
      return true;
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> bool SparseLinearSolver<TOperator,TData,TSolver>::preSolve()
   {
      return true;
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::addInhomogeneous()
   {
      for(size_t i = 0; i < this->mRHSData.size(); i++)
      {
         if(this->mInhomogeneous.at(i).nonZeros() > 0)
         {
            Solver::internal::addCorrection(this->mRHSData.at(i), this->mInhomogeneous.at(i));
         }
      }
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::solve()
   {
      int iteration = 0;

      while(iteration >= 0)
      {
         // Set unused modes to zero
         for(int i = 0; i < this->mZeroIdx; ++i)
         {
            this->mSolution.at(i).setZero();
         }

         // Solve other modes
         typename std::map<MHDFloat, VectorSharedSolverType>::iterator sIt = this->mSolver.find(this->mId);
         for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
         {
            #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
               MpiFramework::syncSubComm(MpiFramework::SPECTRAL, i);
            #endif //define QUICC_MPI && defined QUICC_MPISPSOLVE

            //internal::solveWrapper<TOperator,TData>(this->mSolution.at(i), this->mSolver.at(i+start), this->mRHSData.at(i));
            internal::solveWrapper(this->mSolution.at(i), sIt->second.at(i), this->mRHSData.at(i));

            // Stop simulation if solve failed
            if(sIt->second.at(i)->info() != Eigen::Success)
            {
               throw std::logic_error("Sparse direct solve failed!");
            }
         }

         // Callback site for correcting solve solution (for example for influence matrix approach)
         iteration = this->correctSolution(iteration);
      }
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> int SparseLinearSolver<TOperator,TData,TSolver>::correctSolution(const int iteration)
   {
      return -1;
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> bool SparseLinearSolver<TOperator,TData,TSolver>::postSolve()
   {
      return false;
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::zeroSolver()
   {
      // Set solver RHS to zero
      for(unsigned int i = 0; i < this->mRHSData.size(); ++i)
      {
         this->mRHSData.at(i).setZero();
      }
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::initSolver()
   {
      // Loop over matrices
      for(auto it = this->mLHSMatrix.begin(); it != this->mLHSMatrix.end(); ++it)
      {
         this->mSolver.insert(std::make_pair(it->first, VectorSharedSolverType()));

         typename std::map<MHDFloat, VectorSharedSolverType>::iterator sIt = this->mSolver.find(it->first);
         sIt->second.reserve(it->second.size());

         for(size_t i = 0; i < it->second.size(); ++i)
         {
            Eigen::aligned_allocator< TSolver<TOperator> > alloc;
            #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
               MpiFramework::syncSubComm(MpiFramework::SPECTRAL, i);

               std::shared_ptr< TSolver<TOperator> >  solver = std::allocate_shared< TSolver<TOperator> >(alloc, MpiFramework::getSubComm(MpiFramework::SPECTRAL, i));
            #else
               std::shared_ptr< TSolver<TOperator> >  solver = std::allocate_shared< TSolver<TOperator> >(alloc);
            #endif //define QUICC_MPI && defined QUICC_MPISPSOLVE

            sIt->second.push_back(solver);
         }
      }

      // Compute pattern and factorisation
      this->updateSolver();
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::updateSolver()
   {
      for(auto it = this->mLHSMatrix.begin(); it != this->mLHSMatrix.end(); ++it)
      {
         // Compute factorisation
         for(size_t i = 0; i < it->second.size(); i++)
         {
            if(static_cast<int>(i) % this->nSystem() >= this->mZeroIdx)
            {
               #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
                  MpiFramework::syncSubComm(MpiFramework::SPECTRAL, i);
               #endif //define QUICC_MPI && defined QUICC_MPISPSOLVE

               typename std::map<MHDFloat, VectorSharedSolverType>::iterator sIt = this->mSolver.find(it->first);
               // Safety assert to make sur matrix is compressed
               assert(it->second.at(i).isCompressed());

               sIt->second.at(i)->compute(it->second.at(i));

               // Stop simulation if factorization failed
               if(sIt->second.at(i)->info() != Eigen::Success)
               {
                  throw std::logic_error("Matrix factorization failed!");
               }
            }
         }
      }

      // Provide call site for additional updates
      this->postSolverUpdate();
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::postSolverUpdate()
   {
      // Default implementation does nothing
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::initMatrices(const int n)
   {
      this->initMatrices(0, n);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::initMatrices(const MHDFloat id, const int n)
   {
      // Do not reinitialise if work already done by other field
      if(this->mLHSMatrix.count(id) == 0 || this->mLHSMatrix.find(id)->second.size() == 0)
      {
         if(this->mLHSMatrix.count(id) == 0)
         {
            this->mLHSMatrix.insert(std::make_pair(id, std::vector<TOperator>()));
         }

         typename std::map<MHDFloat,std::vector<TOperator> >::iterator it = this->mLHSMatrix.find(id);

         // Reserve space for the LHS matrices
         it->second.reserve(n);

         // Initialise storage
         for(int i = 0; i < n; ++i)
         {
            // Create storage for LHS matrices
            it->second.push_back(TOperator());
         }
      }
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::initSolutions()
   {
      // Nothing to be done in general.
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      // Add storage for RHS data
      this->mRHSData.push_back(TData(rows,cols));
      this->mRHSData.back().setZero();

      // Add storage for solution
      this->mSolution.push_back(TData(rows,cols));
      this->mSolution.back().setZero();

      // Add storage for inhomogeneous boundary value
      this->mInhomogeneous.push_back(Eigen::SparseMatrix<typename TData::Scalar>(rows,cols));
      this->mInhomogeneous.back().setZero();
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::buildOperators(const int idx, const std::map<std::size_t,DecoupledZSparse>& ops, const int size)
   {
      std::map<std::size_t,DecoupledZSparse>::const_iterator iOpA = ops.find(ModelOperator::ImplicitLinear::id());
      this->rLHSMatrix(0, idx).resize(size, size);
      Solver::internal::addOperators(this->rLHSMatrix(0, idx), 1.0, iOpA->second);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> int SparseLinearSolver<TOperator,TData,TSolver>::nSystem() const
   {
      return this->mRHSData.size();
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> TOperator& SparseLinearSolver<TOperator,TData,TSolver>::rLHSMatrix(const MHDFloat id, const int idx)
   {
      return this->mLHSMatrix.find(id)->second.at(idx);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> TData& SparseLinearSolver<TOperator,TData,TSolver>::rRHSData(const int idx)
   {
      return this->mRHSData.at(idx);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> const TData& SparseLinearSolver<TOperator,TData,TSolver>::solution(const int idx) const
   {
      return this->mSolution.at(idx);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> TData& SparseLinearSolver<TOperator,TData,TSolver>::rSolution(const int idx)
   {
      return this->mSolution.at(idx);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> const Eigen::SparseMatrix<typename TData::Scalar>& SparseLinearSolver<TOperator,TData,TSolver>::inhomogeneous(const int idx) const
   {
      return this->mInhomogeneous.at(idx);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> Eigen::SparseMatrix<typename TData::Scalar>& SparseLinearSolver<TOperator,TData,TSolver>::rInhomogeneous(const int idx)
   {
      return this->mInhomogeneous.at(idx);
   }

   namespace internal {

      inline void addOperators(SparseMatrix& mat, const MHDFloat c, const DecoupledZSparse& decMat)
      {
         assert(decMat.real().rows() > 0);
         assert(decMat.real().cols() > 0);
         assert(decMat.imag().size() == 0 || decMat.imag().nonZeros() == 0);

         if(c != 1.0)
         {
            mat += c*decMat.real();
         } else
         {
            mat += decMat.real();
         }
      }

      inline void addOperators(SparseMatrixZ& mat, const MHDFloat c, const DecoupledZSparse& decMat)
      {
         assert(decMat.real().rows() > 0);
         assert(decMat.real().cols() > 0);
         assert(decMat.imag().rows() > 0);
         assert(decMat.imag().cols() > 0);
         assert(decMat.real().rows() == decMat.imag().rows());
         assert(decMat.real().cols() == decMat.imag().cols());

         if(c != 1.0)
         {
            mat += c*decMat.real().cast<MHDComplex>() + c*Math::cI*decMat.imag();
         } else
         {
            mat += decMat.real().cast<MHDComplex>() + Math::cI*decMat.imag();
         }
      }

      inline void addCorrection(DecoupledZMatrix& rVal, const SparseMatrixZ& corr)
      {
         assert(rVal.real().rows() > 0);
         assert(rVal.real().cols() > 0);
         assert(rVal.imag().rows() > 0);
         assert(rVal.imag().cols() > 0);
         assert(rVal.real().rows() == rVal.imag().rows());
         assert(rVal.real().cols() == rVal.imag().cols());

         rVal.real() += corr.real();
         rVal.imag() += corr.imag();
      }

      template <typename TData, typename TCorr> inline void addCorrection(TData& rVal, const TCorr& corr)
      {
         rVal += corr;
      }
   }
}
}

#endif // QUICC_SOLVER_SPARSELINEARSOLVER_HPP
