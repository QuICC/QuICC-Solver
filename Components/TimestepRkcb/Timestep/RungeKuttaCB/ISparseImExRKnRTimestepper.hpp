/**
 * @file ISparseImExRKnRTimestepper.hpp
 * @brief Implementation of base for the templated (coupled) equation
 * timestepper for Implicit-Explicit Runge-Kutta (nR) schemes.
 */

#ifndef QUICC_TIMESTEP_RUNGEKUTTACB_ISPARSEIMEXRKNRTIMESTEPPER_HPP
#define QUICC_TIMESTEP_RUNGEKUTTACB_ISPARSEIMEXRKNRTIMESTEPPER_HPP

// System includes
//
#include <map>
#include <memory>

// Project includes
//
#include "QuICC/Register/Error.hpp"
#include "QuICC/Register/Explicit.hpp"
#include "QuICC/Register/Implicit.hpp"
#include "QuICC/Register/Intermediate.hpp"
#include "QuICC/Timestep/ISparseTimestepper.hpp"
#include "Timestep/RungeKuttaCB/IImExRKCBScheme.hpp"

namespace QuICC {

namespace Timestep {

namespace RungeKuttaCB {

/**
 * @brief Implementation of a templated (coupled) equation timestepper for
 * Implicit-Explicit Runge-Kutta (2R) schemes
 */
template <typename TOperator, typename TData, template <typename> class TSolver>
class ISparseImExRKnRTimestepper
    : public ISparseTimestepper<TOperator, TData, TSolver>
{
public:
   /**
    * @brief Constructor
    *
    * @param start   Starting index (for example without m=0)
    * @param timeId  Solver timing with respect to timestepping
    */
   ISparseImExRKnRTimestepper(const int start, const std::size_t timeId);

   /**
    * @brief Destructor
    */
   virtual ~ISparseImExRKnRTimestepper() = default;

   /**
    * @brief Set timestepper scheme
    */
   void setScheme(SharedIImExRKCBScheme spScheme);

   /**
    * @brief Initialise solution after data was copied
    */
   virtual void initSolutions();

   /**
    * @brief Add RHS and solution data storage
    *
    * @param rows Number of rows of matrix
    * @param cols Number of columns required
    */
   virtual void addStorage(const int rows, const int cols);

   /**
    * @brief Get current timestep fraction
    */
   MHDFloat stepFraction() const;

protected:
   /**
    * @brief Number of substeps
    */
   virtual int steps() const;

   /**
    * @brief Implicit coefficient a for linear operator
    *
    * A = (T + a L)
    */
   virtual MHDFloat aIm(const int step) const;

   /**
    * @brief Timestepping scheme
    */
   SharedIImExRKCBScheme mspScheme;

private:
};

template <typename TOperator, typename TData, template <typename> class TSolver>
ISparseImExRKnRTimestepper<TOperator, TData,
   TSolver>::ISparseImExRKnRTimestepper(const int start,
   const std::size_t timeId) :
    ISparseTimestepper<TOperator, TData, TSolver>(start, timeId)
{}

template <typename TOperator, typename TData, template <typename> class TSolver>
void ISparseImExRKnRTimestepper<TOperator, TData, TSolver>::setScheme(
   SharedIImExRKCBScheme spScheme)
{
   this->mspScheme = spScheme;
}

template <typename TOperator, typename TData, template <typename> class TSolver>
int ISparseImExRKnRTimestepper<TOperator, TData, TSolver>::steps() const
{
   return this->mspScheme->steps();
}

template <typename TOperator, typename TData, template <typename> class TSolver>
MHDFloat
ISparseImExRKnRTimestepper<TOperator, TData, TSolver>::stepFraction() const
{
   return this->mspScheme->cEx(this->mStep);
}

template <typename TOperator, typename TData, template <typename> class TSolver>
MHDFloat ISparseImExRKnRTimestepper<TOperator, TData, TSolver>::aIm(
   const int step) const
{
   return this->mspScheme->aIm(step, step);
}

template <typename TOperator, typename TData, template <typename> class TSolver>
void ISparseImExRKnRTimestepper<TOperator, TData, TSolver>::initSolutions()
{
   for (size_t i = this->mZeroIdx; i < this->nSystem(); i++)
   {
      details::computeSet(this->reg(Register::Intermediate::id()).at(i),
         this->reg(Register::Solution::id()).at(i));

      if (this->mspScheme->useEmbedded())
      {
         details::computeSet(this->reg(Register::Error::id()).at(i),
            this->reg(Register::Solution::id()).at(i));
      }
   }
}

template <typename TOperator, typename TData, template <typename> class TSolver>
void ISparseImExRKnRTimestepper<TOperator, TData, TSolver>::addStorage(
   const int rows, const int cols)
{
   // Assert for non zero rows and columns
   assert(rows > 0);
   assert(cols > 0);

   ISparseTimestepper<TOperator, TData, TSolver>::addStorage(rows, cols);

   // Add additional registers
   std::vector<std::size_t> ids = {Register::Implicit::id(),
      Register::Explicit::id(), Register::Intermediate::id()};

   // Add storage for embedded scheme
   if (this->mspScheme->useEmbedded())
   {
      ids.push_back(Register::Error::id());
   }

   this->addRegister(rows, cols, ids);
}
} // namespace RungeKuttaCB
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_RUNGEKUTTACB_ISPARSEIMEXRKNRTIMESTEPPER_HPP
