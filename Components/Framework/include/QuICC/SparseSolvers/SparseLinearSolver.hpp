/**
 * @file SparseLinearSolver.hpp
 * @brief Implementation of a templated (coupled) linear solver structure
 */

#ifndef QUICC_SOLVER_SPARSELINEARSOLVER_HPP
#define QUICC_SOLVER_SPARSELINEARSOLVER_HPP

// System includes
//
#include <memory>
#include <stdexcept>

// Project includes
//
#include "QuICC/Math/Constants.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/Framework/MpiFramework.hpp"
#include "QuICC/SparseSolvers/SparseSolverBase.hpp"
#include "QuICC/SparseSolvers/SparseLinearSolverTools.hpp"
#include "QuICC/Register/Solution.hpp"
#include "QuICC/Register/Rhs.hpp"
#include "QuICC/Tag/Operator/Lhs.hpp"

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
         virtual ~SparseLinearSolver() = default;

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
          * @brief Update solver after updated solution was copied
          */
         virtual void updateSolutions();

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
         std::size_t nSystem() const;

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
         typedef std::vector<SharedSolverType> VectorSharedSolverType;

         /**
          * @brief Has solver matrix?
          *
          * @param opId Operator ID
          */
         bool hasSolverMatrix(const std::size_t opId) const;

         /**
          * @brief Has solver matrix?
          *
          * @param opId Operator ID
          * @param id   ID of matrix
          */
         bool hasSolverMatrix(const std::size_t opId, const MHDFloat id) const;

         /**
          * @brief Set solver matrix
          *
          * @param opId Operator ID
          * @param id   ID of matrix
          * @param idx Index of the matrix
          */
         TOperator& solverMatrix(const std::size_t, const MHDFloat id, const int idx);

         /**
          * @brief Get solver matrix
          *
          * @param opId Operator ID
          * @param id   ID of matrix
          * @param idx Index of the matrix
          */
         const TOperator& solverMatrix(const std::size_t, const MHDFloat id, const int idx) const;

         /**
          * @brief Get storage register
          */
         std::vector<TData>& reg(const std::size_t);

         /**
          * @brief Get storage register
          */
         const std::vector<TData>& reg(const std::size_t) const;

         /**
          * @brief Add RHS and solution data storage
          *
          * @param rows Number of rows of matrix
          * @param cols Number of columns required
          * @param ids  Register IDs
          */
         void addRegister(const int rows, const int cols, const std::vector<std::size_t>& ids);

         /**
          * @brief Initialise the solver matrices storage
          *
          * @param n Size of matrices
          */
         void initMatrices(const std::size_t opId, const MHDFloat id, const int n);

         /**
          * @brief Initialise the LHS matrices storage
          *
          * @param id Matrix id of matrices
          * @param n  Size of matrices
          */
         void initMatrices(const MHDFloat id, const int n);

         /**
          * @brief Initialise solver for given ID
          *
          * @param opId Operator ID
          */
         void initSolver(const std::size_t opId);

         /**
          * @brief Update solver for given ID
          *
          * @param opId Operator ID
          */
         void updateSolver(const std::size_t opId);

         /**
          * @brief Additional updates after solver update
          */
         virtual void postSolverUpdate();

         /**
          * @brief Correct solution obtained from linear solver
          */
         virtual int correctSolution(const int iteration);

         /**
          * @brief Current operator ID of the solvers
          */
         std::size_t mOpId;

         /**
          * @brief Current ID of the solvers
          */
         MHDFloat mId;

         /**
          * @brief Operators used by solver
          */
         std::map<std::size_t, std::map<MHDFloat, std::vector<TOperator> > >  mSolverMatrix;

         /**
          * @brief Storage for field
          */
         std::map<std::size_t, std::vector<TData> > mStorage;

         /**
          * @brief Storage for inhomogeneous boundary conditions
          */
         std::vector<Eigen::SparseMatrix<typename TData::Scalar> >  mInhomogeneous;

         /**
          * @brief Create sparse solvers
          */
         std::map<std::size_t, std::map<MHDFloat,VectorSharedSolverType> >  mSolver;

      private:
   };

   template <typename TOperator,typename TData,template <typename> class TSolver> SparseLinearSolver<TOperator,TData,TSolver>::SparseLinearSolver(const int start, const std::size_t timeId)
      : SparseSolverBase(start, timeId), mOpId(Tag::Operator::Lhs::id()), mId(0.0)
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

   template <typename TOperator,typename TData,template <typename> class TSolver> std::vector<TData>& SparseLinearSolver<TOperator,TData,TSolver>::reg(const std::size_t id)
   {
      assert(this->mStorage.count(id) > 0);

      return this->mStorage.at(id);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> const std::vector<TData>& SparseLinearSolver<TOperator,TData,TSolver>::reg(const std::size_t id) const
   {
      assert(this->mStorage.count(id) > 0);

      return this->mStorage.at(id);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::addRegister(const int rows, const int cols, const std::vector<std::size_t>& ids)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      for(auto id: ids)
      {
         if(this->mStorage.count(id) == 0)
         {
            this->mStorage.emplace(id, std::vector<TData>());
         }

         this->reg(id).push_back(TData(rows,cols));
         this->reg(id).back().setZero();
      }
    }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::addInhomogeneous()
   {
      for(size_t i = 0; i < this->nSystem(); i++)
      {
         if(this->mInhomogeneous.at(i).nonZeros() > 0)
         {
            Solver::internal::addCorrection(this->reg(Register::Rhs::id()).at(i), this->mInhomogeneous.at(i));
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
            this->reg(Register::Solution::id()).at(i).setZero();
         }

         // Solve other modes
         assert(this->mSolver.count(this->mOpId) > 0);
         assert(this->mSolver.at(this->mOpId).count(this->mId) > 0);
         typename std::map<MHDFloat, VectorSharedSolverType>::iterator sIt = this->mSolver.at(this->mOpId).find(this->mId);
         for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
         {
            #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
               MpiFramework::syncSubComm(MpiFramework::SPECTRAL, i);
            #endif //define QUICC_MPI && defined QUICC_MPISPSOLVE

            assert(this->reg(Register::Solution::id()).size() > i);
            assert(sIt->second.size() > i);
            assert(this->nSystem() > i);

            internal::solveWrapper(this->reg(Register::Solution::id()).at(i), sIt->second.at(i), this->reg(Register::Rhs::id()).at(i));

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
      for(std::size_t i = 0; i < this->nSystem(); ++i)
      {
         this->reg(Register::Rhs::id()).at(i).setZero();
      }
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::initSolver()
   {
      for(auto it = this->mSolverMatrix.begin(); it != this->mSolverMatrix.end(); ++it)
      {
         this->initSolver(it->first);
      }
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::initSolver(const std::size_t opId)
   {
      // Loop over matrices
      assert(this->mSolverMatrix.count(opId) > 0);
      auto&& lhsMatrix = this->mSolverMatrix.at(opId);

      // Create operator map
      if(this->mSolver.count(opId) == 0)
      {
         this->mSolver.emplace(opId, std::map<MHDFloat, VectorSharedSolverType>());
      }

      // Initialize solvers with matrices
      for(auto it = lhsMatrix.begin(); it != lhsMatrix.end(); ++it)
      {
         auto&& solver = this->mSolver.at(opId);
         solver.insert(std::make_pair(it->first, VectorSharedSolverType()));

         typename std::map<MHDFloat, VectorSharedSolverType>::iterator sIt = solver.find(it->first);
         sIt->second.reserve(it->second.size());

         for(size_t i = 0; i < it->second.size(); ++i)
         {
            #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
               MpiFramework::syncSubComm(MpiFramework::SPECTRAL, i);

               std::shared_ptr< TSolver<TOperator> >  solver = std::make_shared< TSolver<TOperator> >(MpiFramework::getSubComm(MpiFramework::SPECTRAL, i));
            #else
               std::shared_ptr< TSolver<TOperator> >  solver = std::make_shared< TSolver<TOperator> >();
            #endif //define QUICC_MPI && defined QUICC_MPISPSOLVE

            sIt->second.push_back(solver);
         }
      }

      // Compute pattern and factorisation
      this->updateSolver(opId);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::updateSolver()
   {
      for(auto it = this->mSolverMatrix.begin(); it != this->mSolverMatrix.end(); ++it)
      {
         this->updateSolver(it->first);
      }
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::updateSolver(const std::size_t opId)
   {
      assert(this->mSolverMatrix.count(opId) > 0);

      auto&& lhsMatrix = this->mSolverMatrix.at(opId);
      for(auto it = lhsMatrix.begin(); it != lhsMatrix.end(); ++it)
      {
         // Compute factorisation
         for(std::size_t i = 0; i < it->second.size(); i++)
         {
            if(i % this->nSystem() >= static_cast<std::size_t>(this->mZeroIdx))
            {
               #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
                  MpiFramework::syncSubComm(MpiFramework::SPECTRAL, i);
               #endif //define QUICC_MPI && defined QUICC_MPISPSOLVE

               assert(this->mSolver.count(opId) > 0);
               auto&& solver = this->mSolver.at(opId);
               typename std::map<MHDFloat, VectorSharedSolverType>::iterator sIt = solver.find(it->first);
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
      this->initMatrices(Tag::Operator::Lhs::id(), 0, n);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::initMatrices(const MHDFloat id, const int n)
   {
      this->initMatrices(Tag::Operator::Lhs::id(), id, n);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::initMatrices(const std::size_t opId, const MHDFloat id, const int n)
   {
      // Create operator ID
      if(this->mSolverMatrix.count(opId) == 0)
      {
         this->mSolverMatrix.emplace(opId, std::map<MHDFloat, std::vector<TOperator>>());
      }

      // Do not reinitialise if work already done by other field
      auto&& lhsMatrix = this->mSolverMatrix.at(opId);
      if(lhsMatrix.count(id) == 0 || lhsMatrix.find(id)->second.size() == 0)
      {
         if(lhsMatrix.count(id) == 0)
         {
            lhsMatrix.insert(std::make_pair(id, std::vector<TOperator>()));
         }

         typename std::map<MHDFloat,std::vector<TOperator> >::iterator it = lhsMatrix.find(id);

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

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::updateSolutions()
   {
      // Nothing to be done in general.
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      // Add additional registers for RHS data and solution
      std::vector<std::size_t> ids = {Register::Solution::id(), Register::Rhs::id()};
      this->addRegister(rows, cols, ids);

      // Add storage for inhomogeneous boundary value
      this->mInhomogeneous.push_back(Eigen::SparseMatrix<typename TData::Scalar>(rows,cols));
      this->mInhomogeneous.back().setZero();
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseLinearSolver<TOperator,TData,TSolver>::buildOperators(const int idx, const std::map<std::size_t,DecoupledZSparse>& ops, const int size)
   {
      std::map<std::size_t,DecoupledZSparse>::const_iterator iOpA = ops.find(ModelOperator::ImplicitLinear::id());
      std::map<std::size_t,DecoupledZSparse>::const_iterator iOpC = ops.find(ModelOperator::Boundary::id());

      auto&& lhsMat = this->rLHSMatrix(0, idx);
      lhsMat.resize(size, size);
      Solver::internal::addOperators(lhsMat, 1.0, iOpA->second);
      Solver::internal::addOperators(lhsMat, 1.0, iOpC->second);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> std::size_t SparseLinearSolver<TOperator,TData,TSolver>::nSystem() const
   {
      if(this->mStorage.count(Register::Rhs::id()) > 0)
      {
         return this->reg(Register::Rhs::id()).size();
      }
      else
      {
         return 0;
      }
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> bool SparseLinearSolver<TOperator,TData,TSolver>::hasSolverMatrix(const std::size_t opId) const
   {
      bool res = (this->mSolverMatrix.count(opId) > 0);
      return res;
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> bool SparseLinearSolver<TOperator,TData,TSolver>::hasSolverMatrix(const std::size_t opId, const MHDFloat id) const
   {
      bool res = (this->hasSolverMatrix(opId) && this->mSolverMatrix.at(opId).count(id) > 0);
      return res;
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> TOperator& SparseLinearSolver<TOperator,TData,TSolver>::solverMatrix(const std::size_t opId, const MHDFloat id, const int idx)
   {
      assert(this->mSolverMatrix.count(opId) > 0);
      assert(this->mSolverMatrix.at(opId).count(id) > 0);
      assert(this->mSolverMatrix.at(opId).at(id).size() > static_cast<std::size_t>(idx));

      return this->mSolverMatrix.at(opId).at(id).at(idx);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> const TOperator& SparseLinearSolver<TOperator,TData,TSolver>::solverMatrix(const std::size_t opId, const MHDFloat id, const int idx) const
   {
      assert(this->mSolverMatrix.count(opId) > 0);
      assert(this->mSolverMatrix.at(opId).count(id) > 0);
      assert(this->mSolverMatrix.at(opId).at(id).size() > static_cast<std::size_t>(idx));

      return this->mSolverMatrix.at(opId).at(id).at(idx);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> TOperator& SparseLinearSolver<TOperator,TData,TSolver>::rLHSMatrix(const MHDFloat id, const int idx)
   {
      return this->solverMatrix(Tag::Operator::Lhs::id(), id, idx);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> TData& SparseLinearSolver<TOperator,TData,TSolver>::rRHSData(const int idx)
   {
      return this->reg(Register::Rhs::id()).at(idx);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> const TData& SparseLinearSolver<TOperator,TData,TSolver>::solution(const int idx) const
   {
      return this->reg(Register::Solution::id()).at(idx);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> TData& SparseLinearSolver<TOperator,TData,TSolver>::rSolution(const int idx)
   {
      return this->reg(Register::Solution::id()).at(idx);
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
} // Solver
} // QuICC

#endif // QUICC_SOLVER_SPARSELINEARSOLVER_HPP
