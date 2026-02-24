/**
 * @file ITimestepperBase.hpp
 * @brief Implementation of a templated (coupled) linear solver structure
 */

#ifndef QUICC_TIMESTEP_PREDICTORCORRECTOR_VIEWS_ITIMESTEPPERBASE_HPP
#define QUICC_TIMESTEP_PREDICTORCORRECTOR_VIEWS_ITIMESTEPPERBASE_HPP

// System includes
//
#include <memory>
#include <stdexcept>

// Project includes
//
#include "QuICC/Tag/Operator/Influence.hpp"
#include "Types/Math.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/Framework/MpiFramework.hpp"
#include "QuICC/Solver/SparseSolver.hpp"
#include "Timestep/PredictorCorrector/Views/details/SparseLinearSolverTools.hpp"
#include "Timestep/PredictorCorrector/Views/details/TimesteppperTools.hpp"
#include "QuICC/Register/Solution.hpp"
#include "QuICC/Register/Rhs.hpp"
#include "QuICC/Tag/Operator/Lhs.hpp"

namespace QuICC {

namespace Timestep {

namespace PredictorCorrector {

namespace Views {

   template <typename TOperator, typename TData, typename TImpl> class ITimestepper;

   /**
    * @brief Implementation of a generic timestepper
    */
   template <typename TOperator, typename TData, typename TImpl> class ITimestepperBase
   {
      public:
         /**
          * @brief Constructor
          */
         ITimestepperBase();

         /**
          * @brief Destructor
          */
         virtual ~ITimestepperBase() = default;

         /**
          * @brief Is operator initialized?
          */
         bool isInitialized() const;

         /**
          * @brief Set operator to initialized
          */
         void setInitialized();

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
          * @brief Solve linear systems
          */
         void solve();

         /**
          * @brief Add RHS and solution data storage
          *
          * @param rows Number of rows of matrix
          * @param cols Number of columns required
          */
         virtual void addStorage(const int rows, const int cols) = 0;

         /**
          * @brief Init storage for inhomogeneous BC
          *
          * @param rows Number of rows of matrix
          * @param cols Number of columns required
          */
         virtual void initInhomogeneous(const int rows, const int cols);

         /**
          * @brief Get inhomogeneous boundary condition
          */
         const Eigen::SparseMatrix<typename TData::Scalar>& inhomogeneous() const;

         /**
          * @brief Set inhomogeneous boundary condition
          */
         Eigen::SparseMatrix<typename TData::Scalar>& rInhomogeneous();

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
         typedef Framework::Selector::SparseSolver<TOperator> SolverType;

         /// Typedef for shared solver
         typedef std::shared_ptr<SolverType> SharedSolverType;

         /**
          * @brief Has solver matrix?
          *
          * @param opId Operator ID
          */
         bool hasLinearOperator(const std::size_t opId) const;

         /**
          * @brief Has solver matrix?
          *
          * @param opId Operator ID
          * @param id   ID of matrix
          */
         bool hasLinearOperator(const std::size_t opId, const MHDFloat id) const;

         /**
          * @brief Set linear operator
          *
          * @param opId Operator ID
          * @param id   ID of matrix
          */
         TOperator& linearOperator(const std::size_t opId, const MHDFloat id);

         /**
          * @brief Get linear operator
          *
          * @param opId Operator ID
          * @param id   ID of matrix
          */
         const TOperator& linearOperator(const std::size_t opId, const MHDFloat id) const;

         /**
          * @brief Get storage register
          *
          * @param id  Register ID
          */
         TData& reg(const std::size_t regId);

         /**
          * @brief Get storage register
          *
          * @param id  Register ID
          */
         const TData& reg(const std::size_t regId) const;

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
          * @param opId Operator ID
          * @param id   Id of matrix
          */
         void initMatrices(const std::size_t opId, const MHDFloat id);

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
          * @brief Flag for operator initialization
          */
         bool mIsInitialized;

         /**
          * @brief Current operator ID of the solvers
          */
         std::size_t mOpId;

         /**
          * @brief Current ID of the solvers
          */
         MHDFloat mId;

      private:
         /// ITimestepper class is a friend
         friend class ITimestepper<TOperator, TData, TImpl>;

         /**
          * @brief Operators used by solver
          */
         std::map<std::size_t, std::map<MHDFloat, TOperator> >  mOperators;

         /**
          * @brief Create sparse solvers
          */
         std::map<std::size_t, std::map<MHDFloat,SharedSolverType> >  mSolver;

         /**
          * @brief Storage for field
          */
         std::map<std::size_t, TData> mStorage;

         /**
          * @brief Storage for inhomogeneous boundary conditions
          */
         Eigen::SparseMatrix<typename TData::Scalar>  mInhomogeneous;
   };

   template <typename TOperator,typename TData,typename TImpl> ITimestepperBase<TOperator,TData,TImpl>::ITimestepperBase()
      : mIsInitialized(false), mOpId(Tag::Operator::Lhs::id()), mId(0.0)
   {
   }

   template <typename TOperator,typename TData,typename TImpl> MHDFloat ITimestepperBase<TOperator,TData,TImpl>::error() const
   {
      return -1.0;
   }

   template <typename TOperator,typename TData,typename TImpl> bool ITimestepperBase<TOperator,TData,TImpl>::finished()
   {
      return true;
   }

   template <typename TOperator,typename TData,typename TImpl> TData& ITimestepperBase<TOperator,TData,TImpl>::reg(const std::size_t id)
   {
      assert(this->mStorage.count(id) > 0);

      return this->mStorage.at(id);
   }

   template <typename TOperator,typename TData,typename TImpl> const TData& ITimestepperBase<TOperator,TData,TImpl>::reg(const std::size_t id) const
   {
      assert(this->mStorage.count(id) > 0);

      return this->mStorage.at(id);
   }

   template <typename TOperator,typename TData,typename TImpl> void ITimestepperBase<TOperator,TData,TImpl>::addRegister(const int rows, const int cols, const std::vector<std::size_t>& ids)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      for(auto id: ids)
      {
         this->mStorage.emplace(id, TData(rows, cols));
      }
    }

   template <typename TOperator,typename TData,typename TImpl> void ITimestepperBase<TOperator,TData,TImpl>::addInhomogeneous()
   {
      if(this->mInhomogeneous.nonZeros() > 0)
      {
         details::addCorrection(this->reg(Register::Rhs::id()), this->mInhomogeneous);
      }
   }

   template <typename TOperator,typename TData,typename TImpl> void ITimestepperBase<TOperator,TData,TImpl>::solve()
   {
      int iteration = 0;

      while(iteration >= 0)
      {
         // Solve other modes
         assert(this->mSolver.count(this->mOpId) > 0);
         assert(this->mSolver.at(this->mOpId).count(this->mId) > 0);
         auto&& spSolver = this->mSolver.at(this->mOpId).find(this->mId)->second;
         details::solveWrapper(this->reg(Register::Solution::id()), spSolver, this->reg(Register::Rhs::id()));

         // Stop simulation if solve failed
         if(spSolver->info() != Eigen::Success)
         {
            throw std::logic_error("Sparse direct solve failed!");
         }

         // Callback site for correcting solve solution (for example for influence matrix approach)
         iteration = this->correctSolution(iteration);
      }
   }

   template <typename TOperator,typename TData,typename TImpl> int ITimestepperBase<TOperator,TData,TImpl>::correctSolution(const int iteration)
   {
      return -1;
   }

   template <typename TOperator,typename TData,typename TImpl> void ITimestepperBase<TOperator,TData,TImpl>::zeroSolver()
   {
      this->reg(Register::Rhs::id()).setZero();
   }

   template <typename TOperator,typename TData,typename TImpl> void ITimestepperBase<TOperator,TData,TImpl>::initSolver()
   {
      this->initSolver(Tag::Operator::Lhs::id());
      if(this->hasLinearOperator(Tag::Operator::Influence::id()))
      {
         this->initSolver(Tag::Operator::Influence::id());
      }
   }

   template <typename TOperator,typename TData,typename TImpl> void ITimestepperBase<TOperator,TData,TImpl>::initSolver(const std::size_t opId)
   {
      // Loop over matrices
      assert(this->mOperators.count(opId) > 0);
      auto&& lhsMatrix = this->mOperators.at(opId);

      // Create operator map
      if(this->mSolver.count(opId) == 0)
      {
         this->mSolver.emplace(opId, std::map<MHDFloat, SharedSolverType>());
      }

      // Initialize solvers with matrices
      for(auto it = lhsMatrix.begin(); it != lhsMatrix.end(); ++it)
      {
         auto&& solvers = this->mSolver.at(opId);

         auto  spSolver = std::make_shared<SolverType>();
         solvers.emplace(it->first, spSolver);
      }

      // Compute pattern and factorisation
      this->updateSolver(opId);
   }

   template <typename TOperator,typename TData,typename TImpl> void ITimestepperBase<TOperator,TData,TImpl>::updateSolver()
   {
      this->updateSolver(Tag::Operator::Lhs::id());
   }

   template <typename TOperator,typename TData,typename TImpl> void ITimestepperBase<TOperator,TData,TImpl>::updateSolver(const std::size_t opId)
   {
      assert(this->mOperators.count(opId) > 0);
      assert(this->mSolver.count(opId) > 0);

      auto&& lhsMatrix = this->mOperators.at(opId);
      for(auto it = lhsMatrix.begin(); it != lhsMatrix.end(); ++it)
      {
         // Compute factorisation
         auto&& solver = this->mSolver.at(opId);
         typename std::map<MHDFloat, SharedSolverType>::iterator sIt = solver.find(it->first);
         // Safety assert to make sur matrix is compressed
         assert(it->second.isCompressed());

         sIt->second->compute(it->second);

         // Stop simulation if factorization failed
         if(sIt->second->info() != Eigen::Success)
         {
            throw std::logic_error("Matrix factorization failed!");
         }
      }

      // Provide call site for additional updates
      this->postSolverUpdate();
   }

   template <typename TOperator,typename TData,typename TImpl> void ITimestepperBase<TOperator,TData,TImpl>::postSolverUpdate()
   {
      // Default implementation does nothing
   }

   template <typename TOperator,typename TData,typename TImpl> void ITimestepperBase<TOperator,TData,TImpl>::initMatrices(const std::size_t opId, const MHDFloat id)
   {
      // Create operator ID
      if(this->mOperators.count(opId) == 0)
      {
         this->mOperators.emplace(opId, std::map<MHDFloat, TOperator>());
      }

      // Do not reinitialise if work already done by other field
      auto&& lhsMatrix = this->mOperators.at(opId);
      if(lhsMatrix.count(id) == 0)
      {
         lhsMatrix.emplace(id, TOperator());
      }
   }

   template <typename TOperator,typename TData,typename TImpl> void ITimestepperBase<TOperator,TData,TImpl>::initInhomogeneous(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      // Add storage for inhomogeneous boundary value
      this->mInhomogeneous = Eigen::SparseMatrix<typename TData::Scalar>(rows,cols);
      this->mInhomogeneous.setZero();
   }

   template <typename TOperator,typename TData,typename TImpl> bool ITimestepperBase<TOperator,TData,TImpl>::hasLinearOperator(const std::size_t opId) const
   {
      bool res = (this->mOperators.count(opId) > 0);
      return res;
   }

   template <typename TOperator,typename TData,typename TImpl> bool ITimestepperBase<TOperator,TData,TImpl>::hasLinearOperator(const std::size_t opId, const MHDFloat id) const
   {
      bool res = (this->hasLinearOperator(opId) && this->mOperators.at(opId).count(id) > 0);
      return res;
   }

   template <typename TOperator,typename TData,typename TImpl> TOperator& ITimestepperBase<TOperator,TData,TImpl>::linearOperator(const std::size_t opId, const MHDFloat id)
   {
      assert(this->mOperators.count(opId) > 0);
      assert(this->mOperators.at(opId).count(id) > 0);

      return this->mOperators.at(opId).at(id);
   }

   template <typename TOperator,typename TData,typename TImpl> const TOperator& ITimestepperBase<TOperator,TData,TImpl>::linearOperator(const std::size_t opId, const MHDFloat id) const
   {
      assert(this->mOperators.count(opId) > 0);
      assert(this->mOperators.at(opId).count(id) > 0);

      return this->mOperators.at(opId).at(id);
   }

   template <typename TOperator,typename TData,typename TImpl> const Eigen::SparseMatrix<typename TData::Scalar>& ITimestepperBase<TOperator,TData,TImpl>::inhomogeneous() const
   {
      return this->mInhomogeneous;
   }

   template <typename TOperator,typename TData,typename TImpl> Eigen::SparseMatrix<typename TData::Scalar>& ITimestepperBase<TOperator,TData,TImpl>::rInhomogeneous()
   {
      return this->mInhomogeneous;
   }

   template <typename TOperator,typename TData,typename TImpl> bool ITimestepperBase<TOperator,TData,TImpl>::isInitialized() const
   {
      return this->mIsInitialized;
   }

   template <typename TOperator,typename TData,typename TImpl> void  ITimestepperBase<TOperator,TData,TImpl>::setInitialized()
   {
      this->mIsInitialized = true;
   }

} // Views
} // PredictorCorrector
} // Timestep
} // QuICC

#endif // QUICC_TIMESTEP_PREDICTORCORRECTOR_VIEWS_ITIMESTEPPERBASE_HPP
