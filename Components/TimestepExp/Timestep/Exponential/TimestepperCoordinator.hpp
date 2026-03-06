/**
 * @file TimestepCoordinator.hpp
 * @brief Implementation of the exponential timestep coordinator
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_TIMESTEPCOORDINATOR_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_TIMESTEPCOORDINATOR_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Debug/DebuggerMacro.h"
#include "Timestep/Exponential/Tags.hpp"
#include "Timestep/Exponential/TimestepperInfo.hpp"
#ifdef QUICC_DEBUG
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/ModelOperator/Coordinator.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#endif
#include "Types/Typedefs.hpp"
#include "View/ViewDense.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

   template <template <class,class,typename> class TStepper, typename TImpl> class TimestepperCoordinator;

   /**
    * @brief Implementation of the base for a general sparse solver coordinator
    */
   template <template <class,class,typename> class TStepper> class TimestepperCoordinator<TStepper, base_t>
   {
      public:
         /**
          * @brief Constructor
          */
         TimestepperCoordinator();

         /**
          * @brief Destructor
          */
         virtual ~TimestepperCoordinator() = default;

         /**
          * @brief Initialise coordinator
          *
          * @param sizeInfo  vector of size information
          */
         template <typename TScheme>
         void init(const MHDFloat dt, const std::vector<TimestepperInfo>& info, std::shared_ptr<TScheme> spScheme);

         /**
          * @brief Clear the RHS data of all solvers
          */
         void clearSolvers();

         /**
          * @brief Get current time
          */
         MHDFloat stepFraction();

         /**
          * @brief Solve all the linear systems
          */
         void solveSystems();

         /**
          * @brief Get error measure
          */
         MHDFloat error() const;

         /**
          * @brief Check if step is done
          */
         bool finishedStep() const;

         /**
          * @brief Update timestep
          */
         void updateTimestep(const MHDFloat dt);

         /**
          * @brief Update time dependencies
          */
         void updateMatrices();

         /**
          * @brief Update RHS
          */
         void updateRhs(const TimestepperInfo& info, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& rhs);

         /**
          * @brief Update solution
          */
         void updateSolution(const TimestepperInfo& info, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& sol);

         /**
          * @brief Update Solution with corrections
          */
         void updateSolution(const TimestepperInfo& info, const std::vector<std::tuple<MHDComplex, int, int>>& corr);

         /**
          * @brief Get solution
          */
         void getSolution(View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& sol, const TimestepperInfo& info);

      protected:
         /**
          * @brief Flag to signal end of computation
          */
         bool mFinished;

         /**
          * @brief Computation error
          */
         MHDFloat mError;

      private:
         /**
          * @brief Current timestep
          */
         MHDFloat mDt;

   };

   template <template <class,class,typename> class TStepper> TimestepperCoordinator<TStepper, base_t>::TimestepperCoordinator()
      : mFinished(false), mError(-1.0)
   {
   }

   template <template <class,class,typename> class TStepper> MHDFloat TimestepperCoordinator<TStepper, base_t>::error() const
   {
      return this->mError;
   }

   template <template <class,class,typename> class TStepper> bool TimestepperCoordinator<TStepper, base_t>::finishedStep() const
   {
      return this->mFinished;
   }

   template <template <class, class, typename> class TStepper>
   void TimestepperCoordinator<TStepper, base_t>::updateTimestep(const MHDFloat dt)
   {
      this->mDt = dt;
   }

   template <template <class,class,typename> class TStepper> template <typename TScheme> void TimestepperCoordinator<TStepper,base_t>::init(const MHDFloat dt, const std::vector<TimestepperInfo>& infos, std::shared_ptr<TScheme> spScheme)
   {
      this->mDt = dt;
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper, base_t>::updateSolution(const TimestepperInfo& info, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& sol)
   {
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper, base_t>::updateSolution(const TimestepperInfo& info, const std::vector<std::tuple<MHDComplex, int, int>>& corr)
   {
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper, base_t>::updateRhs(const TimestepperInfo& info, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& rhs)
   {
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper, base_t>::getSolution(View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& sol, const TimestepperInfo& info)
   {
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper, base_t>::clearSolvers()
   {
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper, base_t>::updateMatrices()
   {
   }

   template <template <class,class,typename> class TStepper> MHDFloat TimestepperCoordinator<TStepper,  base_t>::stepFraction()
   {
      return -42.0;
   }

   template <template <class,class,typename> class TStepper> void TimestepperCoordinator<TStepper,base_t>::solveSystems()
   {
   }

} // Exponential
} // Timestep
} // QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_TIMESTEPCOORDINATOR_HPP
