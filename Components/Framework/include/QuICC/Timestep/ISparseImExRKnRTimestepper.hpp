/** 
 * @file ISparseImExRKnRTimestepper.hpp
 * @brief Implementation of base for the templated (coupled) equation timestepper for Implicit-Explicit Runge-Kutta (nR) schemes.
 */

#ifndef QUICC_TIMESTEP_ISPARSEIMEXRKNRTIMESTEPPER_HPP
#define QUICC_TIMESTEP_ISPARSEIMEXRKNRTIMESTEPPER_HPP

// Configuration includes
//
#include <memory>
#include <map>

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/Timestep/IImExRKCBScheme.hpp"
#include "QuICC/SparseSolvers/SparseLinearSolver.hpp"

namespace QuICC {

namespace Timestep {

   namespace internal
   {
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
      template <typename TData> void computeSet(TData& z, const MHDFloat a, const TData& y);

      /**
       * @brief Compute z = a*y
       */
      void computeSet(DecoupledZMatrix& z, const MHDFloat a, const DecoupledZMatrix& y);

      /**
       * @brief Compute z = A*y
       */
      template <typename TOperator,typename TData> void computeMV(TData& z, const TOperator& mat, const TData& y);

      /**
       * @brief Compute z = A*y
       */
      void computeMV(DecoupledZMatrix& z, const SparseMatrix& mat, const DecoupledZMatrix& y);

      /**
       * @brief Compute z = a*x + b*y + z
       */
      template <typename TData> void computeAXPBYPZ(TData& z, const MHDFloat a, const TData& x, const MHDFloat b, const TData& y);

      /**
       * @brief Compute z = a*x + b*y + z
       */
      void computeAXPBYPZ(DecoupledZMatrix& z, const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b, const DecoupledZMatrix& y);

      /**
       * @brief Compute y = a*M*x + y
       */
      template <typename TData> void computeAMXPY(TData& y, const SparseMatrix& mat, const MHDFloat a, const TData& x);

      /**
       * @brief Compute y = a*M*x + y
       */
      void computeAMXPY(DecoupledZMatrix& y, const SparseMatrix& mat, const MHDFloat a, const DecoupledZMatrix& x);

      /**
       * @brief Compute z = a*M*x + y + b*z
       */
      template <typename TData> void computeAMXPYPBZ(TData& z, const SparseMatrix& mat, const MHDFloat a, const TData& x, const TData& y, const MHDFloat b);

      /**
       * @brief Compute z = a*M*x + y + b*z
       */
      void computeAMXPYPBZ(DecoupledZMatrix& z, const SparseMatrix& mat, const MHDFloat a, const DecoupledZMatrix& x, const DecoupledZMatrix& y, const MHDFloat b);

      /**
       * @brief Compute z = a*M*x + b*y + z
       */
      template <typename TData> void computeAMXPBYPZ(TData& z, const SparseMatrix& mat, const MHDFloat a, const TData& x, const MHDFloat b, const TData& y);

      /**
       * @brief Compute z = a*M*x + b*y + z
       */
      void computeAMXPBYPZ(DecoupledZMatrix& z, const SparseMatrix& mat, const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b, const DecoupledZMatrix& y);

      /**
       * @brief Compute z = a*M*x + b*y + M*z
       */
      template <typename TData> void computeAMXPBYPMZ(TData& z, const SparseMatrix& mat, const MHDFloat a, const TData& x, const MHDFloat b, const TData& y);

      /**
       * @brief Compute z = a*M*x + b*y + M*z
       */
      void computeAMXPBYPMZ(DecoupledZMatrix& z, const SparseMatrix& mat, const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b, const DecoupledZMatrix& y);

      /**
       * @brief Compute y = a*x + y
       */
      template <typename TData> void computeAXPY(TData& y, const MHDFloat a, const TData& x);

      /**
       * @brief Compute y = a*x + y
       */
      void computeAXPY(DecoupledZMatrix& y, const MHDFloat a, const DecoupledZMatrix& x);

      /**
       * @brief Compute y = x + a*y
       */
      template <typename TData> void computeXPAY(TData& y, const TData& x, const MHDFloat a);

      /**
       * @brief Compute y = x + a*y
       */
      void computeXPAY(DecoupledZMatrix& y, const DecoupledZMatrix& x, const MHDFloat a);

      /**
       * @brief Compute z = M*y
       */
      template <typename TOperator,typename TData> void computeMY(TData& z, const TOperator& mat, const TData& y);

      /**
       * @brief Compute z = M*y
       */
      void computeMY(DecoupledZMatrix& z, const SparseMatrix& mat, const DecoupledZMatrix& y);

      /**
       * @brief Compute error
       */
      template <typename TData> void computeError(MHDFloat& err, const TData& x, const TData& y);

      /**
       * @brief Compute error
       */
      void computeError(MHDFloat& err, const DecoupledZMatrix& x, const DecoupledZMatrix& y);

   }

   /**
    * @brief Implementation of a templated (coupled) equation timestepper for Implicit-Explicit Runge-Kutta (2R) schemes
    */
   template <typename TOperator,typename TData,template <typename> class TSolver> class ISparseImExRKnRTimestepper: public Solver::SparseLinearSolver<TOperator,TData,TSolver>
   {
      public:
         /**
          * @brief Enum for registers
          */
         enum RegisterId {
            // 
            IMPLICIT_REGISTER = 0,
            // 
            EXPLICIT_REGISTER,
            // 
            SOLUTION_REGISTER,
            // 
            ERROR_REGISTER,
         };

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
         virtual ~ISparseImExRKnRTimestepper();

         /**
          * @brief Set timestepper scheme
          */
         void setScheme(SharedIImExRKCBScheme spScheme);

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
          * @brief Update the LHS matrix with new timedependence
          */
         void updateTimeMatrix(const MHDFloat dt);

         /**
          * @brief Set RHS matrix at t_n
          *
          * @param idx Index of the matrix
          */
         TOperator& rRHSMatrix(const int idx);

         /**
          * @brief Add RHS and solution data storage
          * 
          * @param rows Number of rows of matrix
          * @param cols Number of columns required
          */
         virtual void addStorage(const int rows, const int cols);

         /**
          * @brief Build the scheme operators
          *
          * @param idx  Solver index
          * @param ops  Operators for the timestepper
          */
         virtual void buildOperators(const int idx, const std::map<std::size_t, DecoupledZSparse>& ops, const MHDFloat dt, const int size);

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
         MHDFloat stepFraction() const;

      protected:
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
         RegisterId  mRegisterId;

         /**
          * @brief RHS operator
          */
         std::vector<TOperator>   mRHSMatrix;

         /**
          * @brief Mass matrix operator
          */
         std::vector<SparseMatrix>   mMassMatrix;

         /**
          * @brief Storage for implicit solution piece
          */
         std::vector<TData>  mImSolution;

         /**
          * @brief Storage for explicit solution piece
          */
         std::vector<TData>  mExSolution;

         /**
          * @brief Storage for intermediate solution
          */
         std::vector<TData>  mIntSolution;

         /**
          * @brief Storage for lower order embedded scheme solution
          */
         std::vector<TData>  mErrSolution;

         /**
          * @brief Timestepping scheme
          */
         SharedIImExRKCBScheme   mspScheme;

      private:
   };

   template <typename TOperator,typename TData,template <typename> class TSolver> ISparseImExRKnRTimestepper<TOperator,TData,TSolver>::ISparseImExRKnRTimestepper(const int start, const std::size_t timeId)
      : Solver::SparseLinearSolver<TOperator,TData,TSolver>(start, timeId), mHasExplicit(true), mStep(0), mDt(-1.0), mError(-1.0), mRegisterId(IMPLICIT_REGISTER)
   {
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> ISparseImExRKnRTimestepper<TOperator,TData,TSolver>::~ISparseImExRKnRTimestepper()
   {
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> MHDFloat ISparseImExRKnRTimestepper<TOperator,TData,TSolver>::error() const
   {
      return this->mError;
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> bool ISparseImExRKnRTimestepper<TOperator,TData,TSolver>::finished()
   {
      return (this->mStep == 0);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> MHDFloat ISparseImExRKnRTimestepper<TOperator,TData,TSolver>::stepFraction() const
   {
      return this->mspScheme->cEx(this->mStep);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void  ISparseImExRKnRTimestepper<TOperator,TData,TSolver>::updateTimeMatrix(const MHDFloat dt)
   {
      // Update stored timestep
      MHDFloat oldDt = this->mDt;
      this->mDt = dt;

      for(int step = 0; step < this->mspScheme->steps(); ++step)
      {
         // Update is only required if aIm is not zero
         if(this->mspScheme->aIm(step,step) != 0.0)
         {
            // Loop over matrices within same step
            for(int i = 0; i < this->nSystem(); ++i)
            {
               // Get the number of nonzero elements in time dependence
               size_t nnz = this->mRHSMatrix.at(i).nonZeros();

               // Update LHS and RHS matrices
               for (size_t k=0; k< static_cast<size_t>(this->mRHSMatrix.at(i).outerSize()); ++k)
               {
                  typename TOperator::InnerIterator lhsIt(this->rLHSMatrix(this->mspScheme->aIm(step,step), i),k);
                  for(typename TOperator::InnerIterator timeIt(this->mRHSMatrix.at(i),k); timeIt; ++timeIt)
                  {
                     // Only keep going if nonzero elements are left
                     if(nnz > 0)
                     {
                        assert(lhsIt.col() == timeIt.col());
                        assert(lhsIt.row() <= timeIt.row());

                        // LHS matrix might have additional nonzero entries
                        while(lhsIt.row() < timeIt.row() && lhsIt)
                        {
                           ++lhsIt;
                        }

                        // Update LHS matrix
                        if(timeIt.row() == lhsIt.row())
                        {
                           // Update values
                           lhsIt.valueRef() += this->mspScheme->aIm(step,step)*(oldDt - this->mDt)*timeIt.value();

                           // Update nonzero counter
                           nnz--;
                        }

                        // Update LHS iterators and counters
                        ++lhsIt;

                     } else
                     {
                        break;
                     }
                  }
               }

               // Abort if some nonzero entries where not updated
               if(nnz != 0)
               {
                  throw std::logic_error("Update of timestepping matrices failed");
               }
            }
         }
      }
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void ISparseImExRKnRTimestepper<TOperator,TData,TSolver>::setScheme(SharedIImExRKCBScheme spScheme)
   {
      this->mspScheme = spScheme;
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void ISparseImExRKnRTimestepper<TOperator,TData,TSolver>::initMatrices(const int n)
   {
      // Initialise base matrices
      for(int i = 0; i < this->mspScheme->steps(); i++)
      {
         Solver::SparseLinearSolver<TOperator,TData,TSolver>::initMatrices(this->mspScheme->aIm(i,i), n);
      }

      // Do not reinitialise if work already done by other field
      if(this->mRHSMatrix.size() == 0)
      {
         // Reserve space for the RHS matrices
         this->mRHSMatrix.reserve(n);

         // Initialise storage for RHS matrices
         for(int i = 0; i < n; ++i)
         {
            // Create storage for LHS matrices
            this->mRHSMatrix.push_back(TOperator());
         }
      }

      // Do not reinitialise if work already done by other field
      if(this->mMassMatrix.size() == 0)
      {
         // Reserve space for the RHS matrices
         this->mMassMatrix.reserve(n);

         // Initialise storage for RHS matrices
         for(int i = 0; i < n; ++i)
         {
            // Create storage for LHS matrices
            this->mMassMatrix.push_back(SparseMatrix());
         }
      }
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void ISparseImExRKnRTimestepper<TOperator,TData,TSolver>::initSolutions()
   {
      for(size_t i = this->mZeroIdx; i < this->mRHSData.size(); i++)
      {
         internal::computeSet(this->mIntSolution.at(i), this->mSolution.at(i));

         if(this->mspScheme->useEmbedded())
         {
            internal::computeSet(this->mErrSolution.at(i), this->mSolution.at(i));
         }
      }
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> TOperator& ISparseImExRKnRTimestepper<TOperator,TData,TSolver>::rRHSMatrix(const int idx)
   {
      return this->mRHSMatrix.at(idx);
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void ISparseImExRKnRTimestepper<TOperator,TData,TSolver>::addStorage(const int rows, const int cols)
   {
      // Assert for non zero rows and columns
      assert(rows > 0);
      assert(cols > 0);

      Solver::SparseLinearSolver<TOperator,TData,TSolver>::addStorage(rows,cols);

      // Add RK storage
      this->mImSolution.push_back(TData(rows,cols));
      this->mExSolution.push_back(TData(rows,cols));
      this->mIntSolution.push_back(TData(rows,cols));

      // Initialize storage for embedded scheme
      if(this->mspScheme->useEmbedded())
      {
         this->mErrSolution.push_back(TData(rows,cols));
      }
   }

   template <typename TOperator,typename TData,template <typename> class TSolver> void ISparseImExRKnRTimestepper<TOperator,TData,TSolver>::buildOperators(const int idx, const std::map<std::size_t,DecoupledZSparse>& ops, const MHDFloat dt, const int size)
   {
      std::map<std::size_t,DecoupledZSparse>::const_iterator iOpA = ops.find(ModelOperator::ImplicitLinear::id());
      std::map<std::size_t,DecoupledZSparse>::const_iterator iOpB = ops.find(ModelOperator::Time::id());
      std::map<std::size_t,DecoupledZSparse>::const_iterator iOpC = ops.find(ModelOperator::Boundary::id());

      // Update timestep
      this->mDt = dt;

      // Set explicit matrix
      this->rRHSMatrix(idx).resize(size, size);
      Solver::internal::addOperators(this->rRHSMatrix(idx), 1.0, iOpA->second);

      // Set mass matrix
      this->mMassMatrix.at(idx).resize(size, size);
      Solver::internal::addOperators(this->mMassMatrix.at(idx), 1.0, iOpB->second);

      // Set implicit matrix
      for(int i = 0; i < this->mspScheme->steps(); ++i)
      {
         // Set LHS matrix
         this->rLHSMatrix(this->mspScheme->aIm(i,i), idx).resize(size, size);
         Solver::internal::addOperators(this->rLHSMatrix(this->mspScheme->aIm(i,i), idx), 1.0, iOpB->second);
         Solver::internal::addOperators(this->rLHSMatrix(this->mspScheme->aIm(i,i), idx), -this->mspScheme->aIm(i,i)*this->mDt, iOpA->second);
         Solver::internal::addOperators(this->rLHSMatrix(this->mspScheme->aIm(i,i), idx), 1.0, iOpC->second);
      }
   }

   namespace internal
   {
      template <typename TData> inline void computeSet(TData& y, const TData& x)
      {
         y = x;
      }

      inline void computeSet(DecoupledZMatrix& y, const DecoupledZMatrix& x)
      {
         y.real() = x.real();

         y.imag() = x.imag();
      }

      template <typename TData> inline void computeSet(TData& y, const MHDFloat a, const TData& x)
      {
         y = a*x;
      }

      inline void computeSet(DecoupledZMatrix& y, const MHDFloat a, const DecoupledZMatrix& x)
      {
         y.real() = a*x.real();

         y.imag() = a*x.imag();
      }

      template <typename TData> inline void computeAXPY(TData& y, const MHDFloat a, const TData& x)
      {
         if(a != 0.0)
         {
            y += a*x;
         }
      }

      inline void computeAXPY(DecoupledZMatrix& y, const MHDFloat a, const DecoupledZMatrix& x)
      {
         if(a != 0.0)
         {
            y.real() += a*x.real();

            y.imag() += a*x.imag();
         }
      }

      template <typename TData> inline void computeXPAY(TData& y, const TData& x, const MHDFloat a)
      {
         if(a != 0.0)
         {
            y = x + a*y;

         } else
         {
            computeSet<TData>(y, x);
         }
      }

      inline void computeXPAY(DecoupledZMatrix& y, const DecoupledZMatrix& x, const MHDFloat a)
      {
         if(a != 0.0)
         {
            y.real() = x.real() + a*y.real();

            y.imag() = x.imag() + a*y.imag();

         } else
         {
            computeSet(y, x);
         }
      }

      template <typename TData> inline void computeXPAYPBZ(TData& z, const TData& x, const MHDFloat a, const TData& y, const MHDFloat b)
      {
         if(a == 0.0)
         {
            z = x + b*z;

         } else if(b == 0.0)
         {
            z = x + a*y;

         } else
         {
            z = x + a*y + b*z;
         }
      }

      inline void computeXPAYPBZ(DecoupledZMatrix& z, const DecoupledZMatrix& x, const MHDFloat a, const DecoupledZMatrix& y, const MHDFloat b)
      {
         if(a == 0.0)
         {
            z.real() = x.real() + b*z.real();

            z.imag() = x.imag() + b*z.imag();

         } else if(b == 0.0)
         {
            z.real() = x.real() + a*y.real();

            z.imag() = x.imag() + a*y.imag();

         } else
         {
            z.real() = x.real() + a*y.real() + b*z.real();

            z.imag() = x.imag() + a*y.imag() + b*z.imag();
         }
      }

      template <typename TData> inline void computeAXPBYPZ(TData& z, const MHDFloat a, const TData& x, const MHDFloat b, const TData& y)
      {
         if(a == 0.0)
         {
            z += b*y;

         } else if(b == 0.0)
         {
            z += a*x;

         } else
         {
            z += a*x + b*y;
         }
      }

      inline void computeAXPBYPZ(DecoupledZMatrix& z, const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b, const DecoupledZMatrix& y)
      {
         if(a == 0.0)
         {
            z.real() += b*y.real();

            z.imag() += b*y.imag();

         } else if(b == 0.0)
         {
            z.real() += a*x.real();

            z.imag() += a*x.imag();

         } else
         {
            z.real() += a*x.real() + b*y.real();

            z.imag() += a*x.imag() + b*y.imag();
         }
      }

      template <typename TData> void computeAMXPY(TData& y, const SparseMatrix& mat, const MHDFloat a, const TData& x)
      {
         if(a != 0.0)
         {
            y += mat*(a*x);
         }
      }

      inline void computeAMXPY(DecoupledZMatrix& y, const SparseMatrix& mat, const MHDFloat a, const DecoupledZMatrix& x)
      {
         if(a != 0.0)
         {
            y.real() += mat*(a*x.real());

            y.imag() += mat*(a*x.imag());
         }
      }

      template <typename TData> void computeAMXPYPBZ(TData& z, const SparseMatrix& mat, const MHDFloat a, const TData& x, const TData& y, const MHDFloat b)
      {
         if(a == 0.0)
         {
            z = y + b*z;

         } else
         {
            z = mat*(a*x) + y + b*z;
         }
      }

      inline void computeAMXPYPBZ(DecoupledZMatrix& z, const SparseMatrix& mat, const MHDFloat a, const DecoupledZMatrix& x, const DecoupledZMatrix& y, const MHDFloat b)
      {
         if(a == 0.0)
         {
            z.real() = y.real() + b*z.real();

            z.imag() = y.imag() + b*z.imag();

         } else
         {
            z.real() = mat*(a*x.real()) + y.real() + b*z.real();

            z.imag() = mat*(a*x.imag()) + y.imag() + b*z.imag();
         }
      }

      template <typename TData> void computeAMXPBYPZ(TData& z, const SparseMatrix& mat, const MHDFloat a, const TData& x, const MHDFloat b, const TData& y)
      {
         if(a == 0.0)
         {
            z += b*y;

         } else if(b == 0.0)
         {
            z += mat*(a*x);

         } else
         {
            z += mat*(a*x) + b*y;
         }
      }

      inline void computeAMXPBYPZ(DecoupledZMatrix& z, const SparseMatrix& mat, const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b, const DecoupledZMatrix& y)
      {
         if(a == 0.0)
         {
            z.real() += b*y.real();

            z.imag() += b*y.imag();

         } else if(b == 0.0)
         {
            z.real() += mat*(a*x.real());

            z.imag() += mat*(a*x.imag());

         } else
         {
            z.real() += mat*(a*x.real()) + b*y.real();

            z.imag() += mat*(a*x.imag()) + b*y.imag();
         }
      }

      template <typename TData> void computeAMXPBYPMZ(TData& z, const SparseMatrix& mat, const MHDFloat a, const TData& x, const MHDFloat b, const TData& y)
      {
         if(a == 0.0)
         {
            z = b*y + mat*z;

         } else if(b == 0.0)
         {
            z = mat*(a*x + z);

         } else
         {
            z = mat*(a*x + z) + b*y;
         }
      }

      inline void computeAMXPBYPMZ(DecoupledZMatrix& z, const SparseMatrix& mat, const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b, const DecoupledZMatrix& y)
      {
         if(a == 0.0)
         {
            z.real() = b*y.real() + mat*z.real();

            z.imag() = b*y.imag() + mat*z.imag();

         } else if(b == 0.0)
         {
            z.real() = mat*(a*x.real() + z.real());

            z.imag() = mat*(a*x.imag() + z.imag());

         } else
         {
            z.real() = mat*(a*x.real() + z.real()) + b*y.real();

            z.imag() = mat*(a*x.imag() + z.imag()) + b*y.imag();
         }
      }

      template <typename TOperator,typename TData> inline void computeMV(TData& y, const TOperator& A, const TData& x)
      {
         y = A*x;
      }

      inline void computeMV(DecoupledZMatrix& y, const SparseMatrix& A, const DecoupledZMatrix& x)
      {
         y.real() = A*x.real();

         y.imag() = A*x.imag();
      }

      template <typename TData> inline void computeError(MHDFloat& err, const TData& x, const TData& y)
      {
         err = std::max(err, ((x.array() - y.array())/(1.0 + x.array().abs())).abs().maxCoeff());
      }

      inline void computeError(MHDFloat& err, const DecoupledZMatrix& x, const DecoupledZMatrix& y)
      {
         err = std::max(err, ((x.real().array() - y.real().array())/(1.0 + x.real().array().abs())).abs().maxCoeff());

         err = std::max(err, ((x.imag().array() - y.imag().array())/(1.0 + x.imag().array().abs())).abs().maxCoeff());
      }

   }
}
}

#endif // QUICC_TIMESTEP_ISPARSEIMEXRKNRTIMESTEPPER_HPP
