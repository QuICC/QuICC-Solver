/**
 * @file SparseImExPCTimestepper.hpp
 * @brief Implementation of a templated (coupled) equation timestepper for Implicit-Explicit Predictor-Corrector schemes.
 */

#ifndef QUICC_TIMESTEP_SPARSEIMEXPCTIMESTEPPER_HPP
#define QUICC_TIMESTEP_SPARSEIMEXPCTIMESTEPPER_HPP

// System includes
//

// Project includes
//
#include "QuICC/Timestep/ISparseTimestepper.hpp"
#include "QuICC/Register/Explicit.hpp"
#include "QuICC/Register/Implicit.hpp"
#include "QuICC/Register/Explicit.hpp"
#include "QuICC/Register/Intermediate.hpp"
#include "QuICC/Register/Error.hpp"
#include "QuICC/Timestep/IImExPCScheme.hpp"

namespace QuICC {

namespace Timestep {

   /**
    * @brief Implementation of a templated (coupled) equation timestepper for Implicit-Explicit Predictor-Corrector schemes
    */
   template <typename TOperator,typename TData,template <typename> class TSolver> class SparseImExPCTimestepper: public ISparseTimestepper<TOperator,TData,TSolver>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          * @param timeId  Solver timing with respect to timestepping
          */
         SparseImExPCTimestepper(const int start, const std::size_t timeId);

         /**
          * @brief Destructor
          */
         virtual ~SparseImExPCTimestepper() = default;

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
          * @brief Get current timestep fraction
          */
         virtual MHDFloat stepFraction() const;

         /**
          * @brief Prepare fields for implicit solve
          */
         virtual bool preSolve();

         /**
          * @brief Work on fields after implicit solve
          *
          * @param step    Current substep
          */
         virtual bool postSolve();

         /**
          * @brief Add RHS and solution data storage
          *
          * @param rows Number of rows of matrix
          * @param cols Number of columns required
          */
         virtual void addStorage(const int rows, const int cols);

         /**
          * @brief Initialise solution after data was copied
          */
         virtual void initSolutions();

         /**
          * @brief Update solver after solution was updated
          */
         virtual void updateSolutions();

      protected:
         /**
          * @brief Timestepping scheme
          */
         SharedIImExPCScheme   mspScheme;

      private:
   };

   template <typename TOperator,typename TData,template <typename> class TSolver> SparseImExPCTimestepper<TOperator,TData,TSolver>::SparseImExPCTimestepper(const int start, const std::size_t timeId)
      : ISparseTimestepper<TOperator,TData,TSolver>(start, timeId)
   {
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseImExPCTimestepper<TOperator,TData,TSolver>::setScheme(SharedIImExPCScheme spScheme)
   {
      this->mspScheme = spScheme;
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> int SparseImExPCTimestepper<TOperator,TData,TSolver>::steps() const
   {
      return this->mspScheme->steps();
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> MHDFloat SparseImExPCTimestepper<TOperator,TData,TSolver>::stepFraction() const
   {
      return this->mspScheme->cEx(this->mStep);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> MHDFloat SparseImExPCTimestepper<TOperator,TData,TSolver>::aIm(const int step) const
   {
      return this->mspScheme->aIm(step);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseImExPCTimestepper<TOperator,TData,TSolver>::initSolutions()
   {
      for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
      {
         internal::computeSet(this->reg(Register::Intermediate::id()).at(i), this->reg(Register::Solution::id()).at(i));

         if(this->mspScheme->useEmbedded())
         {
            internal::computeSet(this->reg(Register::Error::id()).at(i), this->reg(Register::Solution::id()).at(i));
         }
      }
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseImExPCTimestepper<TOperator,TData,TSolver>::updateSolutions()
   {
      for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
      {
         internal::computeSet(this->reg(Register::Intermediate::id()).at(i), this->reg(Register::Solution::id()).at(i));

         if(this->mspScheme->useEmbedded())
         {
            internal::computeSet(this->reg(Register::Error::id()).at(i), this->reg(Register::Solution::id()).at(i));
         }
      }
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void SparseImExPCTimestepper<TOperator,TData,TSolver>::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      ISparseTimestepper<TOperator,TData,TSolver>::addStorage(rows,cols);

      // Add additional registers
      std::vector<std::size_t> ids = {Register::Error::id(), Register::Explicit::id(), Register::Intermediate::id()};
      this->addRegister(rows, cols, ids);

      // Register for influence kernels
      ids = {Register::Influence::id()};
      this->addRegister(1, 1, ids);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> bool SparseImExPCTimestepper<TOperator,TData,TSolver>::preSolve()
   {
      const bool hasInfluence = this->hasSolverMatrix(Tag::Operator::Influence::id());

      // Use influence matrix
      if(hasInfluence)
      {
         const bool isFirstPass = (this->mOpId == Tag::Operator::Lhs::id());
         if(isFirstPass)
         {
            this->mOpId = Tag::Operator::Influence::id();
            this->mId = 0.0;

            return true;
         } else
         {
            this->mOpId = Tag::Operator::Lhs::id();
         }
      }

      // First step
      if(this->mStep == 0)
      {
         // Reset error
         if(this->mspScheme->useEmbedded())
         {
            this->mError = 0.0;
         }

         // Build RHS
         MHDFloat aIm = (1.0 - this->mspScheme->aIm(this->mStep))*this->mDt;
         MHDFloat aMass = 1.0;
         MHDFloat aN = this->mDt;
         for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
         {
            if(this->mHasExplicit)
            {
               internal::computeSet(this->reg(Register::Explicit::id()).at(i), -1.0, this->reg(Register::Rhs::id()).at(i));
               internal::computeSet(this->reg(Register::Rhs::id()).at(i), aN, this->reg(Register::Explicit::id()).at(i));
            }
            internal::computeAMXPY(this->reg(Register::Rhs::id()).at(i), this->mMassMatrix.at(i), aMass, this->reg(Register::Intermediate::id()).at(i));
            internal::computeAMXPY(this->reg(Register::Rhs::id()).at(i), this->mRHSMatrix.at(i), aIm, this->reg(Register::Intermediate::id()).at(i));
         }

         this->mId = this->mspScheme->aIm(this->mStep);

         // Include inhomogeneous boundary conditions
         this->addInhomogeneous();
      }
      else
      {
         MHDFloat aNold = -this->mspScheme->aIm(this->mStep)*this->mDt;
         MHDFloat aNnew = this->mspScheme->aIm(this->mStep)*this->mDt;
         for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
         {
            if(this->mHasExplicit)
            {
               internal::computeSet(this->reg(Register::Error::id()).at(i), aNold, this->reg(Register::Explicit::id()).at(i));
               internal::computeSet(this->reg(Register::Explicit::id()).at(i), -1.0, this->reg(Register::Rhs::id()).at(i));
               internal::computeSet(this->reg(Register::Rhs::id()).at(i), aNnew, this->reg(Register::Explicit::id()).at(i));
            }
            internal::computeXPAY(this->reg(Register::Rhs::id()).at(i), this->reg(Register::Error::id()).at(i), 1.0);
         }

         this->mId = this->mspScheme->aIm(this->mStep);
      }

      return true;
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> bool SparseImExPCTimestepper<TOperator,TData,TSolver>::postSolve()
   {
      const bool hasInfluence = this->hasSolverMatrix(Tag::Operator::Influence::id());

      if(hasInfluence)
      {
         const bool isFirstPass = (this->mOpId == Tag::Operator::Influence::id());

         if(isFirstPass)
         {
            // Apply quasi-inverse
            for(std::size_t i = this->mZeroIdx; i < this->nSystem(); i++)
            {
               internal::computeMV(this->reg(Register::Rhs::id()).at(i), this->mMassMatrix.at(i), this->reg(Register::Solution::id()).at(i));
            }

            return true;
         }
         else
         {
            // Correct solution with green's function
            for(std::size_t i = this->mZeroIdx; i < this->nSystem(); i++)
            {
               internal::computeInfluenceCorrection(this->reg(Register::Solution::id()).at(i), this->reg(Register::Influence::id()).at(i));
            }
         }
      }

      if(this->mStep == 0)
      {
         // Store predictor solution
         for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
         {
            internal::computeSet(this->reg(Register::Intermediate::id()).at(i), this->reg(Register::Solution::id()).at(i));
         }

         this->mStep += 1;
      }
      else
      {
         // Build corrected solution
         for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
         {
            internal::computeSet(this->reg(Register::Error::id()).at(i), this->reg(Register::Solution::id()).at(i));
            internal::computeXPAY(this->reg(Register::Solution::id()).at(i), this->reg(Register::Intermediate::id()).at(i), 1.0);
            internal::computeSet(this->reg(Register::Intermediate::id()).at(i), this->reg(Register::Solution::id()).at(i));
         }

         // Compute error
         if(this->mspScheme->useEmbedded())
         {
            for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
            {
               internal::computeErrorFromDiff(this->mError, this->reg(Register::Error::id()).at(i), this->reg(Register::Intermediate::id()).at(i));
            }
         }

         this->mStep += 1;

         // Check if we are done
         if(this->mStep == this->steps())
         {
            this->mStep = 0;
         }
      }

      return false;
   }

} // Timestep
} // QuICC

#endif // QUICC_TIMESTEP_SPARSEIMEXPCTIMESTEPPER_HPP
