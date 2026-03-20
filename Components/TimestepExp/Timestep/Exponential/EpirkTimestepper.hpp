/**
 * @file EpirkTimestepper.hpp
 * @brief Implementation of a templated (coupled) equation timestepper for
 * EPIRK exponential schemes.
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_EPIRKTIMESTEPPER_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_EPIRKTIMESTEPPER_HPP

// System includes
//
#include <atomic>
#include <memory>

// Project includes
//
#include "QuICC/Register/Implicit.hpp"
#include "QuICC/Register/Solution.hpp"
#include "QuICC/Register/Rhs.hpp"
#include "QuICC/Register/Intermediate.hpp"
#include "QuICC/Register/Temporary.hpp"
#include "QuICC/Register/Coordinator.hpp"
#include "QuICC/Tag/Operator/Lhs.hpp"
#include "QuICC/Tag/Operator/Qi.hpp"
#include "QuICC/Solver/SparseSolver.hpp"
#include "Timestep/Exponential/IExpScheme.hpp"
#include "Timestep/Exponential/Kiops.hpp"
#include "Timestep/Exponential/IomKrylov.hpp"
#include "Timestep/Exponential/HighamExponential.hpp"
#include "Timestep/Exponential/AugmentedJacobianFunctor.hpp"
#include "Timestep/Exponential/details/TimesteppperTools.hpp"
#include "View/ViewDense.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

/**
 * @brief Implementation of a templated (coupled) equation timestepper for
 * EPIRK exponential schemes.
 */
template <typename TOperator, typename TData, typename TImpl>
class EpirkTimestepper
{
public:
   /// Typedef for explicit exponential functor
   typedef HighamExponential ExpFunctor;

   /// Typedef for action of augmented Jacobiab functor
   typedef AugmentedJacobianFunctor JacobianFunctor;

   /// Typedef for Krylov functor
   typedef IomKrylov<JacobianFunctor> KrylovFunctor;

   /// Typedef for Phi functor
   typedef Kiops<KrylovFunctor, ExpFunctor> PhiFunctor;

   /**
    * @brief Constructor
    */
   EpirkTimestepper();

   /**
    * @brief Destructor
    */
   virtual ~EpirkTimestepper() = default;

   /**
    * @brief Is operator initialized?
    */
   bool isInitialized() const;

   /**
    * @brief Set operator to initialized
    */
   void setInitialized();

   /**
    * @brief Set timestepper scheme
    */
   void setScheme(std::shared_ptr<IExpScheme> spScheme);

   /**
    * @brief Initialize the Phi solver
    */
   void initPhi(std::shared_ptr<AugmentedJacobianFunctor> pJac);

   /**
    * @brief Add data storage
    *
    * @param rows Number of rows of matrix
    * @param cols Number of columns required
    */
   void addStorage(const int rows, const int cols);

   /**
    * @brief Initialise the solver matrices storage
    */
   void initMatrices(const std::vector<std::size_t>& matIds, const std::size_t startRow);

   /**
    * @brief Initialise solver
    */
   void initSolver();

   /**
    * @brief Reset solver for new calculation
    */
   void resetSolver();

   /**
    * @brief Add data to register
    *
    * @param data New values
    * @param startRow Start row
    */
   void addData(const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& rhs, const std::size_t startRow, const std::size_t regId, const std::size_t col);

   /**
    * @brief Apply quasi-inverse and add to register
    *
    * @param data New values
    * @param startRow Start row
    */
   void addQiData(const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& rhs, const std::size_t startRow, const std::size_t regId, const std::size_t col);

   /**
    * @brief Enforce boundary conditions on register
    *
    * @param rows Rows of data
    * @param cols Cols of data
    * @param startRow Start row
    */
   void enforceBoundaryConditions(const std::size_t rows, const std::size_t cols, const std::size_t startRow, const std::size_t regId, const std::size_t col);

   /**
    * @brief Get data from register
    *
    * @param sol solution values
    * @param startRow Start row
    */
   void getData(View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& sol, const std::size_t startRow, const std::size_t regId, const std::size_t col);

   /**
    * @brief Correct data from register
    *
    * @param corr corrections to solution data
    * @param startRow Start row
    */
   void correctData(const std::vector<std::tuple<MHDComplex,int,int>>& sol, const int rows, const int cols, const std::size_t startRow, const std::size_t regId, const std::size_t col);

   /**
    * @brief Set initial solution
    *
    * @param sol initial solution data
    * @param startRow Start row
    */
   void setSolution(const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& sol, const std::size_t startRow);

   /**
    * @brief Update solver after solution was updated
    */
   void updateSolutions();

   /**
    * @brief Compute timestep
    */
   void stepForward();

   /**
    * @brief Build the scheme operators
    *
    * @param ops  Operators for the timestepper
    */
   void buildOperators(const std::map<std::size_t, std::map<std::size_t, std::pair<int, DecoupledZSparse>>>& ops,
      const MHDFloat dt, const std::size_t startRow);

protected:
   /// Typedef for shared solver
   typedef Framework::Selector::SparseSolver<TOperator> SolverType;

   /// Typedef for shared solver
   typedef std::shared_ptr<SolverType> SharedSolverType;

   /**
    * @brief Get storage register
    *
    * @param id  Register ID
    */
   TData& reg(const std::size_t id);

   /**
    * @brief Get storage register
    *
    * @param id  Register ID
    */
   const TData& reg(const std::size_t id) const;

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
    * @param opId    Operator ID
    * @param matId   Id of matrix
    */
   void initMatrices(const std::size_t opId, const std::size_t matId);

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
    * @brief Timestepping scheme
    */
   std::shared_ptr<IExpScheme> mspScheme;

   /**
    * @brief KIOPS algorithm
    */
   std::unique_ptr<PhiFunctor> mpPhi;

   /**
    * @brief Flag for operator initialization
    */
   bool mIsInitialized;

private:
   /**
    * @brief Krylov convergence tolerance
    */
   const double mcKrylovTol;

   /**
    * @brief Krylov orthogonalization order
    */
   const int mcKrylovOrder;

   /**
    * @brief Phi convergence tolerance
    */
   const double mcPhiTol;

   /**
    * @brief Phi delta
    */
   const double mcPhiDelta;

   /**
    * @brief Phi minimum size of Krylov subspace
    */
   const int mcPhiMmin;

   /**
    * @brief Phi maximum size of Krylov subspace
    */
   const int mcPhiMmax;

   /**
    * @brief Current timestep
    */
   MHDFloat mDt;

   /**
    * @brief RHS storage ID
    */
   std::size_t mRhsId;

   /**
    * @brief Column to use in RHS storage
    */
   std::size_t mRhsCol;

   /**
    * @brief Column to use in RHS storage
    */
   int mKrylovM;

   /**
    * @brief Storage for field
    */
   std::map<std::size_t, TData> mStorage;

   /**
    * @brief Operators: opID -> startRow -> mat
    */
   std::map<std::size_t, std::map<std::size_t, TOperator> >  mOperators;

   /**
    * @brief Create sparse solvers: opID -> startRow -> solver
    */
   std::map<std::size_t, std::map<std::size_t, SharedSolverType> >  mSolver;

   /**
    * @brief Shared Augmented Jacobian
    */
   std::shared_ptr<AugmentedJacobianFunctor> mpJac;
};

template <typename TOperator, typename TData, typename TImpl>
EpirkTimestepper<TOperator, TData, TImpl>::EpirkTimestepper()
   : mIsInitialized(false), mcKrylovTol(1e-12), mcKrylovOrder(2), mcPhiTol(1e-12), mcPhiDelta(1.4), mcPhiMmin(10), mcPhiMmax(128), mDt(-1), mRhsId(Register::Rhs::id()), mRhsCol(1), mKrylovM(10)
{}

template <typename TOperator, typename TData, typename TImpl>
void EpirkTimestepper<TOperator, TData, TImpl>::setScheme(
   std::shared_ptr<IExpScheme> spScheme)
{
   this->mspScheme = spScheme;
}

template <typename TOperator, typename TData, typename TImpl>
void EpirkTimestepper<TOperator, TData, TImpl>::stepForward()
{
   std::vector<double> ts = {1};
   auto& matW = this->reg(Register::Intermediate::id());
   auto& matU = this->reg(Register::Rhs::id());
   this->mpJac->setMatrixHandle(this->reg(Register::Temporary::id()));
   this->mpPhi->setup(ts, matU);
   this->mKrylovM = this->mpPhi->compute(matW, ts, matU, this->mKrylovM, PhiFunctor::Task::I);
   this->mpPhi->printInfo();

   details::computeAXPY(this->reg(Register::Solution::id()), this->mDt, matW);
}

template <typename TOperator,typename TData,typename TImpl> bool EpirkTimestepper<TOperator,TData,TImpl>::isInitialized() const
{
   return this->mIsInitialized;
}

template <typename TOperator,typename TData,typename TImpl> void  EpirkTimestepper<TOperator,TData,TImpl>::setInitialized()
{
   this->mIsInitialized = true;
}

template <typename TOperator,typename TData,typename TImpl> void  EpirkTimestepper<TOperator,TData,TImpl>::initPhi(std::shared_ptr<AugmentedJacobianFunctor> pJac)
{
   assert(this->mDt > 0);

   this->mpJac = pJac;

   auto eFunc = std::make_unique<ExpFunctor>();
   auto kFunc = std::make_unique<KrylovFunctor>(this->mpJac, this->mcKrylovOrder, this->mcKrylovTol);

   this->mpPhi = std::make_unique<PhiFunctor>(std::move(kFunc), std::move(eFunc), this->mcPhiTol, this->mcPhiDelta, this->mcPhiMmin, this->mcPhiMmax);
}

template <typename TOperator, typename TData, typename TImpl>
void EpirkTimestepper<TOperator, TData, TImpl>::addStorage(
   const int rows, const int cols)
{
   // Assert for non zero rows and columns
   assert(rows > 0);
   assert(cols > 0);

   // Add additional registers
   std::vector<std::size_t> ids = {
      Register::Solution::id()
   };
   this->addRegister(rows, 1, ids);

   ids = {
      Register::Rhs::id()
   };
   this->addRegister(rows, 2, ids);

   ids = {
      Register::Intermediate::id()
   };
   this->addRegister(rows, 1, ids);

   ids = {
      Register::Temporary::id()
   };
   this->addRegister(rows, 1, ids);
}

template <typename TOperator,typename TData,typename TImpl> void EpirkTimestepper<TOperator,TData,TImpl>::resetSolver()
{
   this->reg(this->mRhsId).setZero();
}

template <typename TOperator,typename TData,typename TImpl> void EpirkTimestepper<TOperator,TData,TImpl>::initSolver()
{
   this->initSolver(Tag::Operator::Lhs::id());
}

template <typename TOperator,typename TData,typename TImpl> void EpirkTimestepper<TOperator,TData,TImpl>::updateSolutions()
{
}

template <typename TOperator, typename TData, typename TImpl>
void EpirkTimestepper<TOperator, TData, TImpl>::buildOperators(
   const std::map<std::size_t, std::map<std::size_t, std::pair<int, DecoupledZSparse>>>& ops,
   const MHDFloat dt, const std::size_t startRow)
{
   // Update timestep
   this->mDt = dt;

   auto setOps = [](const std::size_t tId, const std::size_t startRow, const auto& ops, auto& tOps)
   {
      // Set RHS Quasi-inverse
      assert(ops.count(tId) > 0);
      auto&& opQs = ops.at(tId);
      assert(tOps.count(tId) > 0);
      auto&& _Qs = tOps.at(tId);
      for(auto&& opData: opQs)
      {
         auto&& matId = opData.first;
         auto&& size = opData.second.first;
         auto&& mat = opData.second.second;

         auto tMatId = matId + startRow;
         assert(_Qs.count(tMatId));
         auto&& _op = _Qs.at(tMatId);
         _op.resize(size, size);
         details::addOperators(_op, 1.0, mat);
      }
   };

   // Set Quasi-Inverse
   setOps(Tag::Operator::Qi::id(),
         startRow,
         ops,
         this->mOperators);

   // Set LHS
   setOps(Tag::Operator::Lhs::id(),
         startRow,
         ops,
         this->mOperators);
}

template <typename TOperator, typename TData, typename TImpl>
void EpirkTimestepper<TOperator, TData, TImpl>::initMatrices(const std::vector<std::size_t>& matIds, const std::size_t startRow)
{
   // Initialise matrices
   for(auto&& id: matIds)
   {
      this->initMatrices(Tag::Operator::Lhs::id(), id + startRow);
      this->initMatrices(Tag::Operator::Qi::id(), id + startRow);
   }
}

template <typename TOperator,typename TData,typename TImpl> TData& EpirkTimestepper<TOperator,TData,TImpl>::reg(const std::size_t id)
{
   assert(this->mStorage.count(id) > 0);

   return this->mStorage.at(id);
}

template <typename TOperator,typename TData,typename TImpl> const TData& EpirkTimestepper<TOperator,TData,TImpl>::reg(const std::size_t id) const
{
   assert(this->mStorage.count(id) > 0);

   return this->mStorage.at(id);
}

template <typename TOperator,typename TData,typename TImpl> void EpirkTimestepper<TOperator,TData,TImpl>::addRegister(const int rows, const int cols, const std::vector<std::size_t>& ids)
{
   // Assert for non zero rows and columns
   assert(rows > 0);
   assert(cols > 0);

   for(auto id: ids)
   {
      this->mStorage.emplace(id, TData(rows, cols));
   }

   for(auto id: ids)
   {
      this->mStorage.at(id).setZero();
   }
}

template <typename TOperator, typename TData, typename TImpl>
void EpirkTimestepper<TOperator, TData, TImpl>::setSolution(const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& sol, const std::size_t startRow)
{
   details::flatten2Real(this->reg(Register::Solution::id()),
         sol, startRow, 0);
}

template <typename TOperator, typename TData, typename TImpl>
void EpirkTimestepper<TOperator, TData, TImpl>::addData(const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& rhs, const std::size_t startRow, const std::size_t regId, const std::size_t col)
{
   details::flattenAXPY(this->reg(regId), 1.0,
         rhs, startRow, col);
}

template <typename TOperator, typename TData, typename TImpl>
void EpirkTimestepper<TOperator, TData, TImpl>::addQiData(const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& rhs, const std::size_t startRow, const std::size_t regId, const std::size_t col)
{
   const auto& rows = rhs.dims()[0];
   const auto& cols = rhs.dims()[1];
   Matrix tmp(2*rows*cols, 1);
   details::flatten2Real(tmp, rhs, 0, 0);

   Eigen::Map<Matrix> tmpRe(tmp.data(), rows, 2*cols);
   auto& regMat = this->reg(regId);
   Eigen::Map<Matrix> out(regMat.data() + col*regMat.rows() + startRow, rows, 2*cols);
   out.noalias() += this->mOperators.at(Tag::Operator::Qi::id()).at(startRow)*tmpRe;
}

template <typename TOperator, typename TData, typename TImpl>
void EpirkTimestepper<TOperator, TData, TImpl>::enforceBoundaryConditions(const std::size_t rows, const std::size_t cols, const std::size_t startRow, const std::size_t regId, const std::size_t col)
{
   auto& regMat = this->reg(regId);
   Eigen::Map<Matrix> regMap(regMat.data() + col*regMat.rows() + startRow, rows, 2*cols);
   Matrix tmp = regMap;

   regMap = this->mSolver.at(Tag::Operator::Lhs::id()).at(startRow)->solve(tmp);
}

template <typename TOperator, typename TData, typename TImpl>
void EpirkTimestepper<TOperator, TData, TImpl>::getData(View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& sol, const std::size_t startRow, const std::size_t regId, const std::size_t col)
{
   details::unflatten2Complex(sol, this->reg(regId), startRow, col);
}

template <typename TOperator, typename TData, typename TImpl>
void EpirkTimestepper<TOperator, TData, TImpl>::correctData(const std::vector<std::tuple<MHDComplex,int,int>>& sol, const int rows, const int cols, const std::size_t startRow, const std::size_t regId, const std::size_t col)
{
   details::addCorrection(this->reg(regId),
         sol, rows, cols, startRow, col);
}

template <typename TOperator,typename TData,typename TImpl> void EpirkTimestepper<TOperator,TData,TImpl>::initMatrices(const std::size_t opId, const std::size_t matId)
{
   // Create operator ID
   if(this->mOperators.count(opId) == 0)
   {
      this->mOperators.emplace(opId, std::map<std::size_t, TOperator>());
   }

   // Do not reinitialise if work already done by other field
   auto&& lhsMatrix = this->mOperators.at(opId);
   if(lhsMatrix.count(matId) == 0)
   {
      lhsMatrix.emplace(matId, TOperator());
   }
}

template <typename TOperator,typename TData,typename TImpl> void EpirkTimestepper<TOperator,TData,TImpl>::initSolver(const std::size_t opId)
{
   // Loop over matrices
   assert(this->mOperators.count(opId) > 0);
   auto&& lhsMatrix = this->mOperators.at(opId);

   // Create solver vector
   if(this->mSolver.count(opId) == 0)
   {
      this->mSolver.emplace(opId, std::map<std::size_t, SharedSolverType>());
   }

   // Initialize solvers with matrices
   auto&& solvers = this->mSolver.at(opId);
   for(auto&& op: lhsMatrix)
   {
      auto  spSolver = std::make_shared<SolverType>();
      solvers.emplace(op.first, spSolver);
   }

   // Compute pattern and factorisation
   this->updateSolver(opId);
}

template <typename TOperator,typename TData,typename TImpl> void EpirkTimestepper<TOperator,TData,TImpl>::updateSolver(const std::size_t opId)
{
   assert(this->mOperators.count(opId) > 0);
   assert(this->mSolver.count(opId) > 0);

   auto&& lhsMatrix = this->mOperators.at(opId);
   auto&& solver = this->mSolver.at(opId);
   for(auto&& op: lhsMatrix)
   {
      auto&& matId = op.first;
      auto&& mat = op.second;
      auto&& sol = *solver.at(matId);

      // Safety assert to make sur matrix is compressed
      assert(mat.isCompressed());

      // Compute factorisation
      sol.compute(mat);

      // Stop simulation if factorization failed
      if(sol.info() != Eigen::Success)
      {
         throw std::logic_error("Matrix factorization failed!");
      }
   }
}

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_EPIRKTIMESTEPPER_HPP
