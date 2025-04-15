/**
 * @file ITimestepper.hpp
 * @brief Implementation of base for the templated (coupled) equation
 * timestepper
 */

#ifndef QUICC_TIMESTEP_PREDICTORCORRECTOR_VIEWS_ITIMESTEPPER_HPP
#define QUICC_TIMESTEP_PREDICTORCORRECTOR_VIEWS_ITIMESTEPPER_HPP

// System includes
//
#include <Eigen/Dense>
#include <map>
#include <memory>
#include <set>

// Project includes
//
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/SplitBoundary.hpp"
#include "QuICC/ModelOperator/SplitBoundaryValue.hpp"
#include "QuICC/ModelOperator/SplitImplicitLinear.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/Register/Implicit.hpp"
#include "QuICC/Register/Influence.hpp"
#include "QuICC/Tag/Operator/Rhs.hpp"
#include "Timestep/PredictorCorrector/Views/ITimestepperBase.hpp"
#include "Timestep/PredictorCorrector/Views/details/TimesteppperTools.hpp"
#include "QuICC/Tag/Operator/Influence.hpp"

namespace QuICC {

namespace Timestep {

namespace PredictorCorrector {

namespace Views {

/**
 * @brief Implementation of a templated (coupled) equation timestepper
 */
template <typename TOperator, typename TData, typename TImpl>
class ITimestepper
    : public ITimestepperBase<TOperator, TData, TImpl>
{
public:
   /**
    * @brief Constructor
    *
    * @param timeId  Solver timing with respect to timestepping
    */
   ITimestepper(const std::size_t timeId);

   /**
    * @brief Destructor
    */
   virtual ~ITimestepper() = default;

   /**
    * @brief Initialise the solver matrices storage
    */
   virtual void initMatrices();

   /**
    * @brief Update the LHS matrix with new timedependence
    */
   void updateTimeMatrix(const MHDFloat dt);

   /**
    * @brief Build the scheme operators
    *
    * @param ops  Operators for the timestepper
    */
   virtual void buildOperators(const std::map<std::size_t, DecoupledZSparse>& ops, const MHDFloat dt,
      const int size);

   /**
    * @brief Finished timestep?
    */
   MHDFloat error() const;

   /**
    * @brief Finished timestep?
    */
   bool finished();

   /**
    * @brief Get current timestep fraction
    */
   virtual MHDFloat stepFraction() const = 0;

protected:
   using ITimestepperBase<TOperator, TData, TImpl>::initMatrices;

   /**
    * @brief Number of substeps
    */
   virtual int steps() const = 0;

   /**
    * @brief
    */
   virtual void postSolverUpdate();

   /**
    * @brief Implicit coefficient a for linear operator
    *
    * A = (T + a L)
    */
   virtual MHDFloat aIm(const int step) const = 0;

   /**
    * @brief Explicit calculation took place?
    */
   bool mHasExplicit;

   /**
    * @brief Current substep
    */
   int mStep;

   /**
    * @brief Current timestep
    */
   MHDFloat mDt;

   /**
    * @brief Timestep error
    */
   MHDFloat mError;

   /**
    * @brief ID of the register to use
    */
   std::size_t mRegisterId;

   /**
    * @brief Mass matrix operator
    */
   SparseMatrix mMassMatrix;

private:
};

template <typename TOperator, typename TData, typename TImpl>
ITimestepper<TOperator, TData, TImpl>::ITimestepper(
   const std::size_t timeId) :
    ITimestepperBase<TOperator, TData, TImpl>(timeId),
    mHasExplicit(true),
    mStep(0),
    mDt(-1.0),
    mError(-1.0),
    mRegisterId(Register::Implicit::id())
{}

template <typename TOperator, typename TData, typename TImpl>
MHDFloat ITimestepper<TOperator, TData, TImpl>::error() const
{
   return this->mError;
}

template <typename TOperator, typename TData, typename TImpl>
bool ITimestepper<TOperator, TData, TImpl>::finished()
{
   return (this->mStep == 0);
}

template <typename TOperator, typename TData, typename TImpl>
void ITimestepper<TOperator, TData, TImpl>::updateTimeMatrix(
   const MHDFloat dt)
{
   // Update stored timestep
   MHDFloat oldDt = this->mDt;
   this->mDt = dt;

   // Get list of different IDs
   std::set<MHDFloat> filter;
   for (int step = 0; step < this->steps(); ++step)
   {
      const auto a = this->aIm(step);
      filter.insert(a);
   }

   // Loop over all operator IDs
   for (auto opIt = this->mOperators.begin();
        opIt != this->mOperators.end(); ++opIt)
   {
      // Loop of step IDs
      for (auto a: filter)
      {
         // Update is only required if aIm is not zero
         if (opIt->second.count(a) > 0 && a != 0.0)
         {
            // Get the number of nonzero elements in time dependence
            size_t nnz = this->linearOperator(Tag::Operator::Rhs::id(), 0).nonZeros();

            // Update LHS and RHS matrices
            for (size_t k = 0;
                 k < static_cast<size_t>(this->linearOperator(Tag::Operator::Rhs::id(), 0).outerSize());
                 ++k)
            {
               typename TOperator::InnerIterator lhsIt(
                  this->linearOperator(opIt->first, a), k);
               for (typename TOperator::InnerIterator timeIt(
                       this->linearOperator(Tag::Operator::Rhs::id(), 0), k);
                    timeIt; ++timeIt)
               {
                  // Only keep going if nonzero elements are left
                  if (nnz > 0)
                  {
                     assert(lhsIt.col() == timeIt.col());
                     assert(lhsIt.row() <= timeIt.row());

                     // LHS matrix might have additional nonzero entries
                     while (lhsIt.row() < timeIt.row() && lhsIt)
                     {
                        ++lhsIt;
                     }

                     // Update LHS matrix
                     if (timeIt.row() == lhsIt.row())
                     {
                        // Update values
                        lhsIt.valueRef() +=
                           a * (oldDt - this->mDt) * timeIt.value();

                        // Update nonzero counter
                        nnz--;
                     }

                     // Update LHS iterators and counters
                     ++lhsIt;
                  }
                  else
                  {
                     break;
                  }
               }
            }

            // Abort if some nonzero entries where not updated
            if (nnz != 0)
            {
               throw std::logic_error(
                  "Update of timestepping matrices failed");
            }
         }
      }
   }
}

template <typename TOperator, typename TData, typename TImpl>
void ITimestepper<TOperator, TData, TImpl>::initMatrices()
{
   // Initialise base matrices
   for (int i = 0; i < this->steps(); i++)
   {
      this->initMatrices(Tag::Operator::Lhs::id(), this->aIm(i));
   }

   this->initMatrices(Tag::Operator::Rhs::id(), 0);
}

template <typename TOperator, typename TData, typename TImpl>
void ITimestepper<TOperator, TData, TImpl>::buildOperators(
   const std::map<std::size_t, DecoupledZSparse>& ops,
   const MHDFloat dt, const int size)
{
   std::map<std::size_t, DecoupledZSparse>::const_iterator iOpA =
      ops.find(ModelOperator::ImplicitLinear::id());
   std::map<std::size_t, DecoupledZSparse>::const_iterator iOpB =
      ops.find(ModelOperator::Time::id());
   std::map<std::size_t, DecoupledZSparse>::const_iterator iOpC =
      ops.find(ModelOperator::Boundary::id());

   // Check if equation is solved in split form
   bool isSplit = (ops.count(ModelOperator::SplitImplicitLinear::id()) > 0);
   std::map<std::size_t, DecoupledZSparse>::const_iterator iOpSC =
      ops.find(ModelOperator::SplitBoundary::id());
   std::map<std::size_t, DecoupledZSparse>::const_iterator iOpSA =
      ops.find(ModelOperator::SplitImplicitLinear::id());
   std::map<std::size_t, DecoupledZSparse>::const_iterator iOpSCV =
      ops.find(ModelOperator::SplitBoundaryValue::id());

   // Update timestep
   this->mDt = dt;

   // Set explicit matrix
   this->linearOperator(Tag::Operator::Rhs::id(), 0).resize(size, size);
   details::addOperators(this->linearOperator(Tag::Operator::Rhs::id(), 0), 1.0, iOpA->second);

   // Set mass matrix
   this->mMassMatrix.resize(size, size);
   details::addOperators(this->mMassMatrix, 1.0, iOpB->second);

   // Set implicit matrix
   for (int i = 0; i < this->steps(); ++i)
   {
      MHDFloat a = this->aIm(i);

      // Set LHS matrix
      auto&& lhsMat = this->linearOperator(Tag::Operator::Lhs::id(), a);
      lhsMat.resize(size, size);
      std::cerr << "a = " << a << std::endl;
      details::addOperators(lhsMat, 1.0, iOpB->second);
      std::cerr << "B:" << std::endl;
      std::cerr << lhsMat << std::endl;
      details::addOperators(lhsMat, -a * this->mDt, iOpA->second);
      std::cerr << "B - a A:" << std::endl;
      std::cerr << lhsMat << std::endl;
      details::addOperators(lhsMat, 1.0, iOpC->second);
      std::cerr << "B - aA + BC:" << std::endl;
      std::cerr << lhsMat << std::endl;

      if (isSplit)
      {
         // Initialise influence matrices
         MHDFloat aInf = 0.0;
         this->initMatrices(Tag::Operator::Influence::id(), aInf);

         // Set other LHS matrix
         auto&& infMatrix =
            this->linearOperator(Tag::Operator::Influence::id(), aInf);
         infMatrix.resize(size, size);
         details::addOperators(infMatrix, 1.0, iOpSA->second);
         details::addOperators(infMatrix, 1.0, iOpSC->second);
      }
   }

   if (isSplit)
   {
      // Store information for particular solution
      auto&& infRhs = this->reg(Register::Influence::id());
      details::initInfluence(infRhs, iOpSCV->second, iOpSC->second);
   }
}

template <typename TOperator, typename TData, typename TImpl>
void ITimestepper<TOperator, TData, TImpl>::postSolverUpdate()
{
   const auto lhsId = Tag::Operator::Lhs::id();
   const auto opId = Tag::Operator::Influence::id();

   if (this->mSolver.count(lhsId) > 0 && this->mSolver.count(opId) > 0)
   {
      // Solver for both stages
      auto&& sIt1 = this->mSolver.at(opId).find(0.0)->second;
      auto&& sIt2 = this->mSolver.at(lhsId).find(0.0)->second;

      // Compute green's functions
      auto&& infKernel = this->reg(Register::Influence::id());
      assert(infKernel.real().rows() == infKernel.imag().rows());
      assert(infKernel.real().cols() == infKernel.imag().cols());
      auto rows = infKernel.real().rows();
      auto cols = infKernel.real().cols() / 3;
      TData rhs(rows, cols);
      for (int i = 0; i < cols; i++)
      {
         rhs.real().col(i) = infKernel.real().col(3 * i + 2);
         rhs.imag().col(i) = infKernel.imag().col(3 * i + 2);
      }
      TData sol(rows, cols);
      details::solveWrapper(sol, sIt1, rhs);
      details::computeMV(rhs, this->mMassMatrix, sol);
      details::solveWrapper(sol, sIt2, rhs);
      details::computeSetInfluence(infKernel, sol);
   }
}

} // namespace Views
} // namespace PredictorCorrector
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_PREDICTORCORRECTOR_VIEWS_ITIMESTEPPER_HPP
