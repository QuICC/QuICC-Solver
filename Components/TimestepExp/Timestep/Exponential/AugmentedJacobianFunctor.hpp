/**
 * @file AugmentedJacobianFunctor.hpp
 * @brief Functor for an matrix free action of augmented Jacobian
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_AUGMENTEDJACOBIANFUNCTOR_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_AUGMENTEDJACOBIANFUNCTOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Timestep/Interface.hpp"
#include "Timestep/Exponential/Functors/FunctorData.hpp"
#include "Types/Typedefs.hpp"
#include "Memory/MemoryResource.hpp"
#include "Timestep/Exponential/Functors/DoNothingFunctor.hpp"
#include "Timestep/Exponential/Functors/StepperWrapperFunctor.hpp"
#include "Timestep/Exponential/Functors/ApplyConstraintFunctor.hpp"
#include "Timestep/Exponential/Functors/TransferOutputFunctor.hpp"
#include "Timestep/Exponential/Functors/TransferCorrectionFunctor.hpp"
#include "Timestep/Exponential/Functors/InputFunctor.hpp"
#include "Timestep/Exponential/Functors/GetInputFunctor.hpp"
#include "Timestep/Exponential/Functors/OutputFunctor.hpp"
#include "Timestep/Exponential/Tags.hpp"

namespace QuICC {

namespace Pseudospectral {
   /// Forward declaration
   class Coordinator;
}

namespace Timestep {

namespace Exponential {

   /// Forward declaration
   template <typename T1, typename T2, typename T3> class EpirkTimestepper;

/**
 * @brief Functor for the matrix-free action of augmented Jacobian
 */
class AugmentedJacobianFunctor
{
   public:
      /// Typedef for Field ID to solver field ID
      typedef std::map<SpectralFieldId, std::size_t> IdMap;

      /// Typedef for the Timestepper functor
      using TsFunctor = Functors::StepperWrapperFunctor<EpirkTimestepper<SparseMatrix, Matrix, base_t>>;

      /**
       * @brief ctor
       */
      AugmentedJacobianFunctor(std::shared_ptr<Functors::FunctorData> spData, const MHDFloat dt, const std::size_t regId, const std::size_t regCol, const int fixedIt,  Pseudospectral::Coordinator* pPseudo, std::shared_ptr<IdMap> idMap, std::shared_ptr<Memory::memory_resource> mem);

      /**
       * @brief ctor
       */
      ~AugmentedJacobianFunctor() = default;

      /**
       * @brief Update dt
       */
      void updateTimestep(const MHDFloat dt);

      /**
       * @brief Apply augmented Jacobian
       */
      void operator()(Eigen::Ref<Matrix> out, Eigen::Ref<Matrix> in)  const;

      /**
       * @brief Update B of augmented matrix
       */
      void updateB(const Matrix &matB);

      /**
       * @brief Set data handle
       */
      void setMatrixHandle(Matrix& tmp);

      /**
       * @brief Set timestepper
       */
      void setStepper(std::shared_ptr<TsFunctor> pStepper);

   private:
      using OviewFunctor = Functors::TransferOutputFunctor<TsFunctor>;
      using OcorrFunctor = Functors::TransferCorrectionFunctor<TsFunctor>;
      using IbeforeFunctor = Functors::ApplyConstraintFunctor;
      using IviewFunctor = Functors::GetInputFunctor<TsFunctor>;

      /**
       * @brief Apply Jacobian
       */
      void applyJacobian() const;

      /**
       * @brief Iteration index
       */
      const int mcFixedIt;

      /**
       * @brief Epsilon for complex step
       */
      double mcEps;

      /**
       * @brief Size if A
       */
      int mAn;

      /**
       * @brief Size if B
       */
      int mBn;

      /**
       * @brief Total size
       */
      int mN;

      /**
       * @brief Functor data
       */
      std::shared_ptr<Functors::FunctorData> mspData;

      /**
       * @brief Timestep
       */
      MHDFloat mDt;

      /**
       * @brief Register ID
       */
      std::size_t mRegId;

      /**
       * @brief Register column
       */
      std::size_t mRegCol;

      /**
       * @brief Handle matrix for Jacobian action
       */
      Matrix* mpHandle;

      /**
       * @brief Pointer to Pseudospectral coordinator
       */
      Pseudospectral::Coordinator *mpPseudo;

      /**
       * @brief Matrix B
       */
      Matrix matB;

      /**
       * @brief Shared field ID to solver field id
       */
      std::shared_ptr<IdMap> mpIdMap;

      /**
       * @brief
       */
      std::shared_ptr<Memory::memory_resource> _mem;

      /**
       * @brief Nothing functor
       */
      std::shared_ptr<Functors::DoNothingFunctor> mpNFunc;

      /**
       * @brief Timestepper functor
       */
      std::shared_ptr<TsFunctor> mpTsFunc;

      /**
       * @brief Output view functor
       */
      std::shared_ptr<OviewFunctor> mpOviewFunc;

      /**
       * @brief Output correction functor
       */
      std::shared_ptr<OcorrFunctor> mpOcorrFunc;

      /**
       * @brief Output functor
       */
      std::shared_ptr<Functors::OutputFunctor<OviewFunctor, OcorrFunctor>> mpOutFunc;

      /**
       * @brief Input before functor
       */
      std::shared_ptr<IbeforeFunctor> mpIbefFunc;

      /**
       * @brief Iput view functor
       */
      std::shared_ptr<IviewFunctor> mpIviewFunc;

      /**
       * @brief Input functor
       */
      std::shared_ptr<Functors::InputFunctor<IviewFunctor>> mpInFunc;
};

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_AUGMENTEDJACOBIANFUNCTOR_HPP
