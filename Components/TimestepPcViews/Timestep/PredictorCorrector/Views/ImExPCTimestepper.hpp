/**
 * @file ImExPCTimestepper.hpp
 * @brief Implementation of a templated (coupled) equation timestepper for
 * Implicit-Explicit Predictor-Corrector schemes.
 */

#ifndef QUICC_TIMESTEP_PREDICTORCORRECTOR_VIEWS_IMEXPCTIMSTEPPER_HPP
#define QUICC_TIMESTEP_PREDICTORCORRECTOR_VIEWS_IMEXPCTIMSTEPPER_HPP

// System includes
//

// Project includes
//
#include "QuICC/Register/Error.hpp"
#include "QuICC/Register/Explicit.hpp"
#include "QuICC/Register/Implicit.hpp"
#include "QuICC/Register/Intermediate.hpp"
#include "QuICC/Register/Solution.hpp"
#include "Timestep//PredictorCorrector/Views/ITimestepper.hpp"
#include "Timestep/PredictorCorrector/IImExPCScheme.hpp"
#include "View/ViewDense.hpp"

namespace QuICC {

namespace Timestep {

namespace PredictorCorrector {

namespace Views {

/**
 * @brief Implementation of a templated (coupled) equation timestepper for
 * Implicit-Explicit Predictor-Corrector schemes
 */
template <typename TOperator, typename TData, typename TImpl>
class ImExPCTimestepper
    : public ITimestepper<TOperator, TData, TImpl>
{
public:
   /**
    * @brief Constructor
    *
    * @param timeId  Solver timing with respect to timestepping
    */
   ImExPCTimestepper(const std::size_t timeId);

   /**
    * @brief Destructor
    */
   virtual ~ImExPCTimestepper() = default;

   /**
    * @brief Set timestepper scheme
    */
   void setScheme(SharedIImExPCScheme spScheme);

   /**
    * @brief Number of substeps
    */
   int steps() const final;

   /**
    * @brief Implicit coefficient a for linear operator
    *
    * A = (T + a L)
    */
   MHDFloat aIm(const int step) const final;

   /**
    *
    * @brief Get current timestep fraction
    */
   virtual MHDFloat stepFraction() const;

   /**
    * @brief Prepare fields for implicit solve
    */
   bool preSolve();

   /**
    * @brief Work on fields after implicit solve
    *
    * @param step    Current substep
    */
   bool postSolve();

   /**
    * @brief Add RHS and solution data storage
    *
    * @param rows Number of rows of matrix
    * @param cols Number of columns required
    */
   virtual void addStorage(const int rows, const int cols);

   /**
    * @brief Add to RHS
    *
    * @param rhs RSH values
    * @param startRow Start row
    */
   void addRhs(const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& rhs, const std::size_t startRow);

   /**
    * @brief Extract solution
    *
    * @param sol solution values
    * @param startRow Start row
    */
   void getSolution(View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& sol, const std::size_t startRow);

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

protected:
   /**
    * @brief Timestepping scheme
    */
   SharedIImExPCScheme mspScheme;

private:
};

template <typename TOperator, typename TData, typename TImpl>
ImExPCTimestepper<TOperator, TData, TImpl>::ImExPCTimestepper(
   const std::size_t timeId) :
    ITimestepper<TOperator, TData, TImpl>(timeId)
{}

template <typename TOperator, typename TData, typename TImpl>
void ImExPCTimestepper<TOperator, TData, TImpl>::setScheme(
   SharedIImExPCScheme spScheme)
{
   this->mspScheme = spScheme;
}

template <typename TOperator, typename TData, typename TImpl>
int ImExPCTimestepper<TOperator, TData, TImpl>::steps() const
{
   return this->mspScheme->steps();
}

template <typename TOperator, typename TData, typename TImpl>
MHDFloat
ImExPCTimestepper<TOperator, TData, TImpl>::stepFraction() const
{
   return this->mspScheme->cEx(this->mStep);
}

template <typename TOperator, typename TData, typename TImpl>
MHDFloat ImExPCTimestepper<TOperator, TData, TImpl>::aIm(
   const int step) const
{
   return this->mspScheme->aIm(step);
}

template <typename TOperator, typename TData, typename TImpl>
void ImExPCTimestepper<TOperator, TData, TImpl>::updateSolutions()
{
   details::computeSet(this->reg(Register::Intermediate::id()),
         this->reg(Register::Solution::id()));

   if (this->mspScheme->useEmbedded())
   {
      details::computeSet(this->reg(Register::Error::id()),
            this->reg(Register::Solution::id()));
   }
}

template <typename TOperator, typename TData, typename TImpl>
void ImExPCTimestepper<TOperator, TData, TImpl>::addStorage(
   const int rows, const int cols)
{
   // Assert for non zero rows and columns
   assert(rows > 0);
   assert(cols > 0);

   // Add additional registers
   std::vector<std::size_t> ids = {
      Register::Solution::id(),
      Register::Rhs::id(),
      Register::Error::id(),
      Register::Explicit::id(),
      Register::Intermediate::id()
   };
   this->addRegister(rows, cols, ids);

   // Init storage for inhomogeneous boundary value
   this->initInhomogeneous(rows, cols);

   // Register for influence kernels
   ids = {Register::Influence::id()};
   this->addRegister(1, 1, ids);
}

template <typename TOperator, typename TData, typename TImpl>
bool ImExPCTimestepper<TOperator, TData, TImpl>::preSolve()
{
   const bool hasInfluence =
      this->hasLinearOperator(Tag::Operator::Influence::id());

   // Use influence matrix
   if (hasInfluence)
   {
      const bool isFirstPass = (this->mOpId == Tag::Operator::Lhs::id());
      if (isFirstPass)
      {
         this->mOpId = Tag::Operator::Influence::id();
         this->mId = 0.0;

         return true;
      }
      else
      {
         this->mOpId = Tag::Operator::Lhs::id();
      }
   }

   // First step
   if (this->mStep == 0)
   {
      // Reset error
      if (this->mspScheme->useEmbedded())
      {
         this->mError = 0.0;
      }

      // Build RHS
      MHDFloat aIm = (1.0 - this->mspScheme->aIm(this->mStep)) * this->mDt;
      MHDFloat aMass = 1.0;
      MHDFloat aN = this->mDt;
      if (this->mHasExplicit)
      {
         details::computeSet(this->reg(Register::Explicit::id()), -1.0,
               this->reg(Register::Rhs::id()));
         details::computeSet(this->reg(Register::Rhs::id()), aN,
               this->reg(Register::Explicit::id()));
      }
      details::computeAMXPY(this->reg(Register::Rhs::id()),
            this->mMassMatrix, aMass,
            this->reg(Register::Intermediate::id()));
      details::computeAMXPY(this->reg(Register::Rhs::id()),
            this->linearOperator(Tag::Operator::Rhs::id(), 0), aIm,
            this->reg(Register::Intermediate::id()));

      this->mId = this->mspScheme->aIm(this->mStep);

      // Include inhomogeneous boundary conditions
      this->addInhomogeneous();
   }
   else
   {
      MHDFloat aNold = -this->mspScheme->aIm(this->mStep) * this->mDt;
      MHDFloat aNnew = this->mspScheme->aIm(this->mStep) * this->mDt;
      if (this->mHasExplicit)
      {
         details::computeSet(this->reg(Register::Error::id()), aNold,
            this->reg(Register::Explicit::id()));
         details::computeSet(this->reg(Register::Explicit::id()), -1.0,
            this->reg(Register::Rhs::id()));
         details::computeSet(this->reg(Register::Rhs::id()), aNnew,
            this->reg(Register::Explicit::id()));
      }
      details::computeXPAY(this->reg(Register::Rhs::id()),
         this->reg(Register::Error::id()), 1.0);

      this->mId = this->mspScheme->aIm(this->mStep);
   }

   return true;
}

template <typename TOperator, typename TData, typename TImpl>
bool ImExPCTimestepper<TOperator, TData, TImpl>::postSolve()
{
   const bool hasInfluence =
      this->hasLinearOperator(Tag::Operator::Influence::id());

   if (hasInfluence)
   {
      const bool isFirstPass = (this->mOpId == Tag::Operator::Influence::id());

      if (isFirstPass)
      {
         // Apply quasi-inverse
         details::computeMV(this->reg(Register::Rhs::id()),
               this->mMassMatrix,
               this->reg(Register::Solution::id()));

         return true;
      }
      else
      {
         // Correct solution with green's function
         details::computeInfluenceCorrection(
               this->reg(Register::Solution::id()),
               this->reg(Register::Influence::id()));
      }
   }

   if (this->mStep == 0)
   {
      // Store predictor solution
      details::computeSet(this->reg(Register::Intermediate::id()),
            this->reg(Register::Solution::id()));

      this->mStep += 1;
   }
   else
   {
      // Build corrected solution
      details::computeSet(this->reg(Register::Error::id()),
            this->reg(Register::Solution::id()));
      details::computeXPAY(this->reg(Register::Solution::id()),
            this->reg(Register::Intermediate::id()), 1.0);
      details::computeSet(this->reg(Register::Intermediate::id()),
            this->reg(Register::Solution::id()));

      // Compute error
      if (this->mspScheme->useEmbedded())
      {
         details::computeErrorFromDiff(this->mError,
               this->reg(Register::Error::id()),
               this->reg(Register::Intermediate::id()));
      }

      this->mStep += 1;

      // Check if we are done
      if (this->mStep == this->steps())
      {
         this->mStep = 0;
      }
   }

   return false;
}

template <typename TOperator, typename TData, typename TImpl>
void ImExPCTimestepper<TOperator, TData, TImpl>::setSolution(const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& sol, const std::size_t start)
{
   details::computeSet(this->reg(Register::Solution::id()),
         sol);
}

template <typename TOperator, typename TData, typename TImpl>
void ImExPCTimestepper<TOperator, TData, TImpl>::addRhs(const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& rhs, const std::size_t start)
{
   details::computeAXPY(this->reg(Register::Rhs::id()), 1.0,
         rhs);
}

template <typename TOperator, typename TData, typename TImpl>
void ImExPCTimestepper<TOperator, TData, TImpl>::getSolution(View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& sol, const std::size_t start)
{
   details::computeSet(sol, this->reg(Register::Solution::id()));
}

} // namespace Views
} // namespace PredictorCorrector
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_PREDICTORCORRECTOR_VIEWS_IMEXPCTIMSTEPPER_HPP
