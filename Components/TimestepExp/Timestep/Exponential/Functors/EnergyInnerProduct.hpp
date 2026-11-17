/**
 * @file EnergyInnerProduct.hpp
 * @brief Functor for energy based innner product
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_ENERGYINNERPRODUCT_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_ENERGYINNERPRODUCT_HPP

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
#include "Timestep/Exponential/Functors/TransferOutputFunctor.hpp"
#include "Timestep/Exponential/Functors/OutputFunctor.hpp"
#include "Timestep/Exponential/Functors/EnergyFunctor.hpp"
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

namespace Functors {

/**
 * @brief Functor for energy based inner product
 */
class EnergyInnerProduct
{
   public:
      /// Typedef for Field ID to solver field ID
      typedef std::map<SpectralFieldId, std::size_t> IdMap;

      /// Typedef for the Timestepper functor
      using TsFunctor = Functors::StepperWrapperFunctor<EpirkTimestepper<SparseMatrix, Matrix, base_t>>;

      /**
       * @brief ctor
       */
      EnergyInnerProduct(const std::size_t regId, const std::size_t regCol);

      /**
       * @brief ctor
       */
      void configure(std::shared_ptr<Functors::FunctorData> spData, std::shared_ptr<IdMap> idMap, std::shared_ptr<Memory::memory_resource> mem);

      /**
       * @brief ctor
       */
      ~EnergyInnerProduct() = default;

      /**
       * @brief Compute inner product
       */
      double operator()(const Matrix& u, const int i, const Matrix& v, const int j, const int n)  const;

      /**
       * @brief Compute inner product for multiple vectors
       */
      Array operator()(const Matrix& u, const int i0, const int i1, const Matrix& v, const int j, const int n)  const;

      /**
       * @brief Compute norm based on innner product
       */
      double norm(const Matrix& u, const int i, const int n)  const;

      /**
       * @brief Set data handle
       */
      void setWorkspaceHandle(Matrix& tmp);

      /**
       * @brief Set timestepper
       */
      void setStepper(std::shared_ptr<TsFunctor> pStepper);

   private:
      using OviewFunctor = Functors::TransferOutputFunctor<TsFunctor>;

      /**
       * @brief Initialize
       */
      void init();

      std::vector<SpectralFieldId> qstIds(const SpectralFieldId& id) const;

      /**
       * @brief Size if A
       */
      int mAn;

      /**
       * @brief Functor data
       */
      std::shared_ptr<Functors::FunctorData> mspData;

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
       * @brief Output functor
       */
      std::shared_ptr<Functors::OutputFunctor<OviewFunctor, Functors::DoNothingFunctor>> mpOutFunc;

      /**
       * @brief Output functor
       */
      std::shared_ptr<Functors::EnergyFunctor> mpEFunc;

      /**
       * @brief Field IDs
       */
      std::vector<SpectralFieldId> mFieldIds;

      /**
       * @brief ScalarField
       */
      std::shared_ptr<Framework::Selector::ComplexScalarField> mspField;

};

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_ENERGYINNERPRODUCT_HPP
