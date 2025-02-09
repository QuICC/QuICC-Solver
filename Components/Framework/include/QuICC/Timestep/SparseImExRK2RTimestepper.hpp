/**
 * @file SparseImExRK2RTimestepper.hpp
 * @brief Implementation of a templated (coupled) equation timestepper for Implicit-Explicit Runge-Kutta (2R) schemes.
 *
 * The implementation is based on Cavaglieri & Bewley, "Low-storage implicit/explicit Runge-Kutta schemes for the simulation of stiff high-dimensional ODE systems", JCP, 2015
 */

#ifndef QUICC_TIMESTEP_SPARSEIMEXRK2RTIMESTEPPER_HPP
#define QUICC_TIMESTEP_SPARSEIMEXRK2RTIMESTEPPER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Timestep/ISparseImExRKnRTimestepper.hpp"
#include "QuICC/Register/Solution.hpp"

namespace QuICC {

namespace Timestep {

   /**
    * @brief Implementation of a templated (coupled) equation timestepper for Implicit-Explicit Runge-Kutta (2R) schemes
    */
   template <typename TOperator,typename TData,template <typename> class TSolver> class SparseImExRK2RTimestepper: public ISparseImExRKnRTimestepper<TOperator,TData,TSolver>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          * @param timeId  Solver timing with respect to timestepping
          */
         SparseImExRK2RTimestepper(const int start, const std::size_t timeId);

         /**
          * @brief Destructor
          */
         virtual ~SparseImExRK2RTimestepper() = default;

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

      protected:

      private:
   };

   template <typename TOperator,typename TData,template <typename> class TSolver> SparseImExRK2RTimestepper<TOperator,TData,TSolver>::SparseImExRK2RTimestepper(const int start, const std::size_t timeId)
      : ISparseImExRKnRTimestepper<TOperator,TData,TSolver>(start, timeId)
   {
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> bool SparseImExRK2RTimestepper<TOperator,TData,TSolver>::preSolve()
   {
      if(this->mHasExplicit)
      {
         // Update explicit term with explicit (nonlinear) values
         for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
         {
            internal::computeSet(this->reg(Register::Explicit::id()).at(i), this->reg(Register::Rhs::id()).at(i));
         }
      }

      if(this->mHasExplicit && this->mStep > 0)
      {
         // Update intermediate solution
         MHDFloat bIm = this->mspScheme->bIm(this->mStep)*this->mDt;
         MHDFloat bEx = -this->mspScheme->bEx(this->mStep)*this->mDt;
         for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
         {
            internal::computeAMXPBYPZ(this->reg(Register::Intermediate::id()).at(i), this->mMassMatrix.at(i), bIm, this->reg(Register::Implicit::id()).at(i), bEx, this->reg(Register::Explicit::id()).at(i));
         }

         // Embedded lower order scheme solution
         if(this->mspScheme->useEmbedded())
         {
            bIm = this->mspScheme->bImErr(this->mStep)*this->mDt;
            bEx = -this->mspScheme->bExErr(this->mStep)*this->mDt;
            for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
            {
               internal::computeAMXPBYPZ(this->reg(Register::Error::id()).at(i), this->mMassMatrix.at(i), bIm, this->reg(Register::Implicit::id()).at(i), bEx, this->reg(Register::Explicit::id()).at(i));
            }
         }

         this->mStep += 1;
      }

      // First step
      if(this->mStep == 0)
      {
         // Reset error
         if(this->mspScheme->useEmbedded())
         {
            this->mError = 0.0;
         }

         // Build RHS for implicit term
         for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
         {
            internal::computeMV(this->reg(Register::Rhs::id()).at(i), this->mRHSMatrix.at(i), this->reg(Register::Intermediate::id()).at(i));
         }

         // Set ID for solver
         this->mId = this->mspScheme->aIm(this->mStep, this->mStep);

         this->mRegisterId = Register::Implicit::id();

         return true;

      }
      // Last step has no implicit solve
      else if(this->mStep == this->mspScheme->steps())
      {
         if(this->mRegisterId == Register::Solution::id())
         {
            // Compute error estimate using embedded scheme
            if(this->mspScheme->useEmbedded())
            {
               for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
               {
                  internal::computeSet(this->reg(Register::Rhs::id()).at(i), this->reg(Register::Error::id()).at(i));
               }

               // Set mass matrix ID for solver
               this->mId = 0.0;

               // Set explicit store register for solution
               this->mRegisterId = Register::Error::id();

               // Include inhomogeneous boundary conditions
               this->addInhomogeneous();

               return true;
            } else
            {
               for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
               {
                  internal::computeSet(this->reg(Register::Solution::id()).at(i), this->reg(Register::Intermediate::id()).at(i));
               }

               // Explicit nonlinear term at next step
               this->mHasExplicit = true;

               // Reset step to 0
               this->mStep = 0;

               // Reset register ID
               this->mRegisterId = Register::Implicit::id();

               return false;
            }
         } else if(this->mRegisterId == Register::Error::id())
         {
            for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
            {
               internal::computeSet(this->reg(Register::Solution::id()).at(i), this->reg(Register::Intermediate::id()).at(i));
            }

            for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
            {
               internal::computeError(this->mError, this->reg(Register::Intermediate::id()).at(i), this->reg(Register::Error::id()).at(i));

               internal::computeSet(this->reg(Register::Error::id()).at(i), this->reg(Register::Intermediate::id()).at(i));
            }

            // Explicit nonlinear term at next step
            this->mHasExplicit = true;

            // Reset step to 0
            this->mStep = 0;

            // Reset register ID
            this->mRegisterId = Register::Implicit::id();

            return false;
         } else
         {
            for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
            {
               internal::computeSet(this->reg(Register::Rhs::id()).at(i), this->reg(Register::Intermediate::id()).at(i));
            }

            // Set mass matrix ID for solver
            this->mId = 0.0;

            // Set explicit store register for solution
            this->mRegisterId = Register::Solution::id();

            // Include inhomogeneous boundary conditions
            this->addInhomogeneous();

            return true;
         }
      } else
      {
         if(this->mRegisterId == Register::Implicit::id())
         {
            // Build RHS for implicit term
            MHDFloat aIm = (this->mspScheme->aIm(this->mStep, this->mStep-1) - this->mspScheme->bIm(this->mStep-1))*this->mDt;
            MHDFloat aEx = -(this->mspScheme->aEx(this->mStep, this->mStep-1) - this->mspScheme->bEx(this->mStep-1))*this->mDt;
            for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
            {
               internal::computeAMXPYPBZ(this->reg(Register::Explicit::id()).at(i), this->mMassMatrix.at(i), aIm, this->reg(Register::Implicit::id()).at(i), this->reg(Register::Intermediate::id()).at(i), aEx);
               internal::computeSet(this->reg(Register::Rhs::id()).at(i), this->reg(Register::Explicit::id()).at(i));
            }

            // Set mass matrix ID for solver
            this->mId = 0.0;

            // Set explicit store register for solution
            this->mRegisterId = Register::Explicit::id();

            // Include inhomogeneous boundary conditions
            this->addInhomogeneous();

         } else if(this->mRegisterId == Register::Explicit::id())
         {
            // Update explicit term
            for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
            {
               internal::computeMV(this->reg(Register::Rhs::id()).at(i), this->mRHSMatrix.at(i), this->reg(Register::Explicit::id()).at(i));
            }

            // Set ID for solver
            this->mId = this->mspScheme->aIm(this->mStep, this->mStep);

            // Set implicit store register for solution
            this->mRegisterId = Register::Implicit::id();
         }

         return true;
      }
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> bool SparseImExRK2RTimestepper<TOperator,TData,TSolver>::postSolve()
   {
      if(this->mRegisterId == Register::Explicit::id())
      {
         // Update explicit term
         for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
         {
            internal::computeSet(this->reg(Register::Explicit::id()).at(i), this->reg(Register::Solution::id()).at(i));
         }

         // Loop back to presolve but without new nonlinear term
         this->mHasExplicit = false;

         return true;

      } else if(this->mRegisterId == Register::Implicit::id())
      {
         // Update implicit term
         for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
         {
            internal::computeSet(this->reg(Register::Implicit::id()).at(i), this->reg(Register::Solution::id()).at(i));
         }

      } else if(this->mRegisterId == Register::Solution::id())
      {
         // Update intermediate term
         for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
         {
            internal::computeSet(this->reg(Register::Intermediate::id()).at(i), this->reg(Register::Solution::id()).at(i));
         }

         // Loop back to presolve but without new nonlinear term
         this->mHasExplicit = false;

         return true;

      } else if(this->mRegisterId == Register::Error::id())
      {
         // Update error term
         for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
         {
            internal::computeSet(this->reg(Register::Error::id()).at(i), this->reg(Register::Solution::id()).at(i));
         }

         // Loop back to presolve but without new nonlinear term
         this->mHasExplicit = false;

         return true;
      }

      if(this->mStep == 0)
      {
         // Update intermediate solution
         MHDFloat bIm = this->mspScheme->bIm(this->mStep)*this->mDt;
         MHDFloat bEx = -this->mspScheme->bEx(this->mStep)*this->mDt;
         for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
         {
            internal::computeAMXPBYPMZ(this->reg(Register::Intermediate::id()).at(i), this->mMassMatrix.at(i), bIm, this->reg(Register::Implicit::id()).at(i), bEx, this->reg(Register::Explicit::id()).at(i));
         }

         // Embedded lower order scheme solution
         if(this->mspScheme->useEmbedded())
         {
            bIm = this->mspScheme->bImErr(this->mStep)*this->mDt;
            bEx = -this->mspScheme->bExErr(this->mStep)*this->mDt;
            for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
            {
               internal::computeAMXPBYPMZ(this->reg(Register::Error::id()).at(i), this->mMassMatrix.at(i), bIm, this->reg(Register::Implicit::id()).at(i), bEx, this->reg(Register::Explicit::id()).at(i));
            }
         }

         // Increase step counter
         this->mStep += 1;

         // Loop back to presolve but without new nonlinear term
         this->mHasExplicit = false;

         return true;

      }
      else
      {
         // Prepare solution for new nonlinear term
         MHDFloat aIm = this->mspScheme->aIm(this->mStep, this->mStep)*this->mDt;
         for(size_t i = this->mZeroIdx; i < this->nSystem(); i++)
         {
            internal::computeXPAY(this->reg(Register::Solution::id()).at(i), this->reg(Register::Explicit::id()).at(i), aIm);
         }

         // Next step will have nonlinear term
         this->mHasExplicit = true;

         return false;
      }
   }
}
}

#endif // QUICC_TIMESTEP_SPARSEIMEXRK2RTIMESTEPPER_HPP
