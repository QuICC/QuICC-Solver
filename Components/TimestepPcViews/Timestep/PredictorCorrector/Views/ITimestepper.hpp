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
#include "QuICC/Tag/Operator/Influence.hpp"

namespace QuICC {

namespace Timestep {

namespace PredictorCorrector {

namespace Views {

namespace details {
/**
 * @brief Compute z = y
 */
template <typename TData> void computeSet(TData& z, const TData& y);

/**
 * @brief Compute z = y
 */
void computeSet(DecoupledZMatrix& z, const DecoupledZMatrix& y);

/**
 * @brief Compute z = a*y
 */
template <typename TData>
void computeSet(TData& z, const MHDFloat a, const TData& y);

/**
 * @brief Compute z = a*y
 */
void computeSet(DecoupledZMatrix& z, const MHDFloat a,
   const DecoupledZMatrix& y);

/**
 * @brief Compute z = A*y
 */
template <typename TOperator, typename TData>
void computeMV(TData& z, const TOperator& mat, const TData& y);

/**
 * @brief Compute z = A*y
 */
void computeMV(DecoupledZMatrix& z, const SparseMatrix& mat,
   const DecoupledZMatrix& y);

/**
 * @brief Compute z = a*x + b*y + z
 */
template <typename TData>
void computeAXPBYPZ(TData& z, const MHDFloat a, const TData& x,
   const MHDFloat b, const TData& y);

/**
 * @brief Compute z = a*x + b*y + z
 */
void computeAXPBYPZ(DecoupledZMatrix& z, const MHDFloat a,
   const DecoupledZMatrix& x, const MHDFloat b, const DecoupledZMatrix& y);

/**
 * @brief Compute y = a*M*x + y
 */
template <typename TOperator, typename TData>
void computeAMXPY(TData& y, const TOperator& mat, const MHDFloat a,
   const TData& x);

/**
 * @brief Compute y = a*M*x + y
 */
void computeAMXPY(DecoupledZMatrix& y, const SparseMatrix& mat,
   const MHDFloat a, const DecoupledZMatrix& x);

/**
 * @brief Compute z = a*M*x + y + b*z
 */
template <typename TData>
void computeAMXPYPBZ(TData& z, const SparseMatrix& mat, const MHDFloat a,
   const TData& x, const TData& y, const MHDFloat b);

/**
 * @brief Compute z = a*M*x + y + b*z
 */
void computeAMXPYPBZ(DecoupledZMatrix& z, const SparseMatrix& mat,
   const MHDFloat a, const DecoupledZMatrix& x, const DecoupledZMatrix& y,
   const MHDFloat b);

/**
 * @brief Compute z = a*M*x + b*y + z
 */
template <typename TData>
void computeAMXPBYPZ(TData& z, const SparseMatrix& mat, const MHDFloat a,
   const TData& x, const MHDFloat b, const TData& y);

/**
 * @brief Compute z = a*M*x + b*y + z
 */
void computeAMXPBYPZ(DecoupledZMatrix& z, const SparseMatrix& mat,
   const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b,
   const DecoupledZMatrix& y);

/**
 * @brief Compute z = a*M*x + b*y + M*z
 */
template <typename TData>
void computeAMXPBYPMZ(TData& z, const SparseMatrix& mat, const MHDFloat a,
   const TData& x, const MHDFloat b, const TData& y);

/**
 * @brief Compute z = a*M*x + b*y + M*z
 */
void computeAMXPBYPMZ(DecoupledZMatrix& z, const SparseMatrix& mat,
   const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b,
   const DecoupledZMatrix& y);

/**
 * @brief Compute y = a*x + y
 */
template <typename TData>
void computeAXPY(TData& y, const MHDFloat a, const TData& x);

/**
 * @brief Compute y = a*x + y
 */
void computeAXPY(DecoupledZMatrix& y, const MHDFloat a,
   const DecoupledZMatrix& x);

/**
 * @brief Compute y = x + a*y
 */
template <typename TData>
void computeXPAY(TData& y, const TData& x, const MHDFloat a);

/**
 * @brief Compute y = x + a*y
 */
void computeXPAY(DecoupledZMatrix& y, const DecoupledZMatrix& x,
   const MHDFloat a);

/**
 * @brief Compute error
 */
template <typename TData>
void computeErrorFromDiff(MHDFloat& err, const TData& diff, const TData& ref);

/**
 * @brief Compute error
 */
void computeErrorFromDiff(MHDFloat& err, const DecoupledZMatrix& diff,
   const DecoupledZMatrix& ref);

/**
 * @brief Compute error
 */
template <typename TData>
void computeError(MHDFloat& err, const TData& x, const TData& y);

/**
 * @brief Compute error
 */
void computeError(MHDFloat& err, const DecoupledZMatrix& x,
   const DecoupledZMatrix& y);

/**
 * @brief Initialize influence matrix kernel and boundary condition
 *
 * @param kernel     3 parts, initialize boundary condition and boundary value
 * @param val        Boundary value of inhomogeneous problem
 * @param bc         Boundary conditions
 */
template <typename TData>
void initInfluence(TData& kernel, const DecoupledZSparse& val,
   const DecoupledZSparse& bc);

/**
 * @brief Initialize influence matrix kernel and boundary condition
 *
 * @param kernel     3 parts, initialize boundary condition and boundary value
 * @param val        Boundary value of inhomogeneous problem
 * @param bc         Boundary conditions
 */
void initInfluence(DecoupledZMatrix& kernel, const DecoupledZSparse& val,
   const DecoupledZSparse& bc);

/**
 * @brief Store influence matrix kernel and boundary condition
 *
 * @param reg  Kernel storage with 3 parts: kernel, boundary value and boundary
 * condition
 * @param x    Kernel solution to store
 */
template <typename TData> void computeSetInfluence(TData& reg, const TData& x);

/**
 * @brief Store influence matrix kernel and boundary condition
 *
 * @param reg  Kernel storage with 3 parts: kernel, boundary value and boundary
 * condition
 * @param x    Kernel solution to store
 */
void computeSetInfluence(DecoupledZMatrix& reg, const DecoupledZMatrix& x);

/**
 * @brief Compute correction from influence matrix (Green's function)
 *
 * @param y Split solution to correct
 * @param x 3 parts of influence matrix solution
 */
template <typename TData>
void computeInfluenceCorrection(TData& y, const TData& x);

/**
 * @brief Compute correction from influence matrix (Green's function)
 *
 * @param y Split solution to correct
 * @param x 3 parts of influence matrix solution
 */
void computeInfluenceCorrection(DecoupledZMatrix& y, const DecoupledZMatrix& x);

} // namespace details

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
    * @brief Set RHS matrix at t_n
    */
   TOperator& rRHSMatrix();

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
            size_t nnz = this->rRHSMatrix().nonZeros();

            // Update LHS and RHS matrices
            for (size_t k = 0;
                 k < static_cast<size_t>(this->rRHSMatrix().outerSize());
                 ++k)
            {
               typename TOperator::InnerIterator lhsIt(
                  this->linearOperator(opIt->first, a), k);
               for (typename TOperator::InnerIterator timeIt(
                       this->rRHSMatrix(), k);
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
      this->initMatrices(this->aIm(i));
   }

   if(this->hasLinearOperator(Tag::Operator::Rhs::id()))
   {
      this->initMatrices(Tag::Operator::Rhs::id(), 0);
   }
}

template <typename TOperator, typename TData, typename TImpl>
TOperator& ITimestepper<TOperator, TData, TImpl>::rRHSMatrix()
{
   return this->linearOperator(Tag::Operator::Rhs::id());
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
   this->rRHSMatrix().resize(size, size);
   details::addOperators(this->rRHSMatrix(), 1.0, iOpA->second);

   // Set mass matrix
   this->mMassMatrix.resize(size, size);
   details::addOperators(this->mMassMatrix, 1.0, iOpB->second);

   // Set implicit matrix
   for (int i = 0; i < this->steps(); ++i)
   {
      MHDFloat a = this->aIm(i);

      // Set LHS matrix
      auto&& lhsMat = this->rLHSMatrix(a);
      lhsMat.resize(size, size);
      details::addOperators(lhsMat, 1.0, iOpB->second);
      details::addOperators(lhsMat, -a * this->mDt, iOpA->second);
      details::addOperators(lhsMat, 1.0, iOpC->second);

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

namespace details {
template <typename TData> inline void computeSet(TData& y, const TData& x)
{
   y = x;
}

inline void computeSet(DecoupledZMatrix& y, const DecoupledZMatrix& x)
{
   y.real() = x.real();

   y.imag() = x.imag();
}

template <typename TData>
inline void computeSet(TData& y, const MHDFloat a, const TData& x)
{
   y = a * x;
}

inline void computeSet(DecoupledZMatrix& y, const MHDFloat a,
   const DecoupledZMatrix& x)
{
   y.real() = a * x.real();

   y.imag() = a * x.imag();
}

template <typename TData>
inline void computeAXPY(TData& y, const MHDFloat a, const TData& x)
{
   if (a != 0.0)
   {
      y += a * x;
   }
}

inline void computeAXPY(DecoupledZMatrix& y, const MHDFloat a,
   const DecoupledZMatrix& x)
{
   if (a != 0.0)
   {
      y.real() += a * x.real();

      y.imag() += a * x.imag();
   }
}

template <typename TData>
inline void computeXPAY(TData& y, const TData& x, const MHDFloat a)
{
   if (a != 0.0)
   {
      y = x + a * y;
   }
   else
   {
      computeSet<TData>(y, x);
   }
}

inline void computeXPAY(DecoupledZMatrix& y, const DecoupledZMatrix& x,
   const MHDFloat a)
{
   if (a != 0.0)
   {
      y.real() = x.real() + a * y.real();

      y.imag() = x.imag() + a * y.imag();
   }
   else
   {
      computeSet(y, x);
   }
}

template <typename TData>
inline void computeXPAYPBZ(TData& z, const TData& x, const MHDFloat a,
   const TData& y, const MHDFloat b)
{
   if (a == 0.0)
   {
      z = x + b * z;
   }
   else if (b == 0.0)
   {
      z = x + a * y;
   }
   else
   {
      z = x + a * y + b * z;
   }
}

inline void computeXPAYPBZ(DecoupledZMatrix& z, const DecoupledZMatrix& x,
   const MHDFloat a, const DecoupledZMatrix& y, const MHDFloat b)
{
   if (a == 0.0)
   {
      z.real() = x.real() + b * z.real();

      z.imag() = x.imag() + b * z.imag();
   }
   else if (b == 0.0)
   {
      z.real() = x.real() + a * y.real();

      z.imag() = x.imag() + a * y.imag();
   }
   else
   {
      z.real() = x.real() + a * y.real() + b * z.real();

      z.imag() = x.imag() + a * y.imag() + b * z.imag();
   }
}

template <typename TData>
inline void computeAXPBYPZ(TData& z, const MHDFloat a, const TData& x,
   const MHDFloat b, const TData& y)
{
   if (a == 0.0)
   {
      z += b * y;
   }
   else if (b == 0.0)
   {
      z += a * x;
   }
   else
   {
      z += a * x + b * y;
   }
}

inline void computeAXPBYPZ(DecoupledZMatrix& z, const MHDFloat a,
   const DecoupledZMatrix& x, const MHDFloat b, const DecoupledZMatrix& y)
{
   if (a == 0.0)
   {
      z.real() += b * y.real();

      z.imag() += b * y.imag();
   }
   else if (b == 0.0)
   {
      z.real() += a * x.real();

      z.imag() += a * x.imag();
   }
   else
   {
      z.real() += a * x.real() + b * y.real();

      z.imag() += a * x.imag() + b * y.imag();
   }
}

template <typename TOperator, typename TData>
void computeAMXPY(TData& y, const TOperator& mat, const MHDFloat a,
   const TData& x)
{
   if (a != 0.0)
   {
      y += mat * (a * x);
   }
}

inline void computeAMXPY(DecoupledZMatrix& y, const SparseMatrix& mat,
   const MHDFloat a, const DecoupledZMatrix& x)
{
   if (a != 0.0)
   {
      y.real() += mat * (a * x.real());

      y.imag() += mat * (a * x.imag());
   }
}

template <typename TData>
void computeAMXPYPBZ(TData& z, const SparseMatrix& mat, const MHDFloat a,
   const TData& x, const TData& y, const MHDFloat b)
{
   if (a == 0.0)
   {
      z = y + b * z;
   }
   else
   {
      z = mat * (a * x) + y + b * z;
   }
}

inline void computeAMXPYPBZ(DecoupledZMatrix& z, const SparseMatrix& mat,
   const MHDFloat a, const DecoupledZMatrix& x, const DecoupledZMatrix& y,
   const MHDFloat b)
{
   if (a == 0.0)
   {
      z.real() = y.real() + b * z.real();

      z.imag() = y.imag() + b * z.imag();
   }
   else
   {
      z.real() = mat * (a * x.real()) + y.real() + b * z.real();

      z.imag() = mat * (a * x.imag()) + y.imag() + b * z.imag();
   }
}

template <typename TData>
void computeAMXPBYPZ(TData& z, const SparseMatrix& mat, const MHDFloat a,
   const TData& x, const MHDFloat b, const TData& y)
{
   if (a == 0.0)
   {
      z += b * y;
   }
   else if (b == 0.0)
   {
      z += mat * (a * x);
   }
   else
   {
      z += mat * (a * x) + b * y;
   }
}

inline void computeAMXPBYPZ(DecoupledZMatrix& z, const SparseMatrix& mat,
   const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b,
   const DecoupledZMatrix& y)
{
   if (a == 0.0)
   {
      z.real() += b * y.real();

      z.imag() += b * y.imag();
   }
   else if (b == 0.0)
   {
      z.real() += mat * (a * x.real());

      z.imag() += mat * (a * x.imag());
   }
   else
   {
      z.real() += mat * (a * x.real()) + b * y.real();

      z.imag() += mat * (a * x.imag()) + b * y.imag();
   }
}

template <typename TData>
void computeAMXPBYPMZ(TData& z, const SparseMatrix& mat, const MHDFloat a,
   const TData& x, const MHDFloat b, const TData& y)
{
   if (a == 0.0)
   {
      z = b * y + mat * z;
   }
   else if (b == 0.0)
   {
      z = mat * (a * x + z);
   }
   else
   {
      z = mat * (a * x + z) + b * y;
   }
}

inline void computeAMXPBYPMZ(DecoupledZMatrix& z, const SparseMatrix& mat,
   const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b,
   const DecoupledZMatrix& y)
{
   if (a == 0.0)
   {
      z.real() = b * y.real() + mat * z.real();

      z.imag() = b * y.imag() + mat * z.imag();
   }
   else if (b == 0.0)
   {
      z.real() = mat * (a * x.real() + z.real());

      z.imag() = mat * (a * x.imag() + z.imag());
   }
   else
   {
      z.real() = mat * (a * x.real() + z.real()) + b * y.real();

      z.imag() = mat * (a * x.imag() + z.imag()) + b * y.imag();
   }
}

template <typename TOperator, typename TData>
inline void computeMV(TData& y, const TOperator& A, const TData& x)
{
   y = A * x;
}

inline void computeMV(DecoupledZMatrix& y, const SparseMatrix& A,
   const DecoupledZMatrix& x)
{
   y.real() = A * x.real();

   y.imag() = A * x.imag();
}

template <typename TData>
inline void computeErrorFromDiff(MHDFloat& err, const TData& diff,
   const TData& ref)
{
   err = std::max(err,
      (diff.array() / (1.0 + ref.array().abs())).abs().maxCoeff());
}

inline void computeErrorFromDiff(MHDFloat& err, const DecoupledZMatrix& diff,
   const DecoupledZMatrix& ref)
{
   err = std::max(err, (diff.real().array() / (1.0 + ref.real().array().abs()))
                          .abs()
                          .maxCoeff());

   err = std::max(err, (diff.imag().array() / (1.0 + ref.imag().array().abs()))
                          .abs()
                          .maxCoeff());
}

template <typename TData>
inline void computeError(MHDFloat& err, const TData& x, const TData& y)
{
   err = std::max(err,
      ((x.array() - y.array()) / (1.0 + x.array().abs())).abs().maxCoeff());
}

inline void computeError(MHDFloat& err, const DecoupledZMatrix& x,
   const DecoupledZMatrix& y)
{
   err = std::max(err,
      ((x.real().array() - y.real().array()) / (1.0 + x.real().array().abs()))
         .abs()
         .maxCoeff());

   err = std::max(err,
      ((x.imag().array() - y.imag().array()) / (1.0 + x.imag().array().abs()))
         .abs()
         .maxCoeff());
}

template <typename TData>
inline void initInfluence(TData& y, const DecoupledZSparse& val,
   const DecoupledZSparse& bc)
{
   throw std::logic_error("Not yet implemented");
}

inline void initInfluence(DecoupledZMatrix& y, const DecoupledZSparse& val,
   const DecoupledZSparse& bc)
{
   assert(bc.real().rows() >= val.real().cols());

   // Real value
   int cols = val.real().cols();
   y.real().resize(val.real().rows(), 3 * cols);

   for (int i = 0; i < cols; i++)
   {
      y.real().col(3 * i + 1) = bc.real().row(i).transpose();
      y.real().col(3 * i + 2) = val.real().col(i);
   }

   // Imaginary value
   cols = val.imag().cols();
   y.imag().resize(val.imag().rows(), 3 * cols);

   for (int i = 0; i < cols; i++)
   {
      y.imag().col(3 * i + 1) = bc.real().row(i).transpose();
      y.imag().col(3 * i + 2) = val.imag().col(i);
   }
}

template <typename TData>
inline void computeSetInfluence(TData& reg, const TData& x)
{
   assert(reg.cols() == x.cols() * 3);

   for (int i = 0; i < x.cols(); i++)
   {
      reg.col(3 * i) = x.col(i);
   }
}

inline void computeSetInfluence(DecoupledZMatrix& reg,
   const DecoupledZMatrix& x)
{
   assert(reg.real().cols() == 3 * x.real().cols());
   assert(reg.imag().cols() == 3 * x.imag().cols());

   for (int i = 0; i < x.real().cols(); i++)
   {
      reg.real().col(3 * i) = x.real().col(i);
   }

   for (int i = 0; i < x.imag().cols(); i++)
   {
      reg.imag().col(3 * i) = x.imag().col(i);
   }
}

template <typename TData>
void computeInfluenceCorrection(TData& y, const TData& x)
{
   assert(x.cols() % 3 == 0);

   // compute influence matrix
   int nBC = x.cols() / 3;
   TData mat(nBC, nBC);
   for (int j = 0; j < nBC; j++)
   {
      const auto bc = x.col(3 * j + 1).transpose();
      for (int k = 0; k < nBC; k++)
      {
         mat(j, k) = (bc * x.col(3 * k)).value();
      }
   }
   mat = mat.inverse();

   TData bcVal(nBC, 1);
   for (int j = 0; j < y.cols(); j++)
   {
      for (int k = 0; k < nBC; k++)
      {
         bcVal(k, 0) = (x.col(3 * k + 1).transpose() * y.col(j)).value();
      }
      bcVal = mat * bcVal;
      for (int k = 0; k < nBC; k++)
      {
         y.col(j) -= bcVal(k) * x.col(3 * k);
      }
   }
}

inline void computeInfluenceCorrection(DecoupledZMatrix& y,
   const DecoupledZMatrix& x)
{
   assert(x.real().cols() % 3 == 0);
   assert(x.imag().cols() % 3 == 0);

   // compute influence matrix
   int nBC = x.real().cols() / 3;
   Matrix mat(nBC, nBC);
   for (int j = 0; j < nBC; j++)
   {
      const auto bc = x.real().col(3 * j + 1).transpose();
      for (int k = 0; k < nBC; k++)
      {
         mat(j, k) = bc * x.real().col(3 * k);
      }
   }
   mat = mat.inverse();

   Matrix bcVal(nBC, 1);
   for (int j = 0; j < y.real().cols(); j++)
   {
      for (int k = 0; k < nBC; k++)
      {
         auto bc = x.real().col(3 * k + 1).transpose();
         bcVal(k, 0) = bc * y.real().col(j);
      }
      bcVal = mat * bcVal;
      for (int k = 0; k < nBC; k++)
      {
         y.real().col(j) -= bcVal(k) * x.real().col(3 * k);
      }
   }

   for (int j = 0; j < y.imag().cols(); j++)
   {
      for (int k = 0; k < nBC; k++)
      {
         auto bc = x.imag().col(3 * k + 1).transpose();
         bcVal(k, 0) = bc * y.imag().col(j);
      }
      bcVal = mat * bcVal;
      for (int k = 0; k < nBC; k++)
      {
         y.imag().col(j) -= bcVal(k) * x.imag().col(3 * k);
      }
   }
}

} // namespace details

} // namespace Views
} // namespace PredictorCorrector
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_PREDICTORCORRECTOR_VIEWS_ITIMESTEPPER_HPP
