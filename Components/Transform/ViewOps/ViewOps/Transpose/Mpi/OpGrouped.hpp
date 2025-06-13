/**
 * @file OpGrouped.hpp
 * @brief Transpose operations on Views
 */
#pragma once

// External includes
//

// Project includes
//
#include "Operator/Unary.hpp"
#include "Profiler/Interface.hpp"
#include "View/View.hpp"
#include "ViewOps/Transpose/Mpi/CommGrouped.hpp"
#include "ViewOps/Transpose/Mpi/Coordinates.hpp"
#include "ViewOps/Transpose/StructArray.hpp"
#include "ViewOps/Transpose/Tags.hpp"


namespace QuICC {
/// @brief namespace for Transpose type operations
namespace Transpose {
/// @brief namespace for Mpi backends
namespace Mpi {

using namespace QuICC::Operator;

/// @brief Transpose operator
/// @tparam Tout
/// @tparam Tin
template <class Tout, class Tin, class Perm>
class OpGrouped : public UnaryBaseOp<OpGrouped<Tout, Tin, Perm>, Tout, Tin>
{
   using ScalarType = typename Tin::value_type::ScalarType;

public:
   /// @brief Constructor
   OpGrouped(std::shared_ptr<Memory::memory_resource> mem)
   {
      _comm = std::make_unique<CommGrouped<ScalarType>>(mem);
   };
   /// @brief default constructor
   OpGrouped() = delete;
   /// @brief default dtor
   ~OpGrouped() = default;

private:
   /// @brief action implementation
   /// @param out output View
   /// @param in input View
   void applyImpl(Tout& out, const Tin& in);
   /// @brief give access to base class
   friend UnaryBaseOp<OpGrouped<Tout, Tin, Perm>, Tout, Tin>;
   /// @brief communicator object
   std::unique_ptr<CommGrouped<ScalarType>> _comm;
   /// @brief maximum group size
   std::uint64_t _groupSize;
};


template <class Tout, class Tin, class Perm>
void OpGrouped<Tout, Tin, Perm>::applyImpl(Tout& out, const Tin& in)
{
   assert(in.size() >= 1);
   assert(out.size() >= 1);
   assert(in.size() == out.size());

   Profiler::RegionFixture<4> fix("Transpose::Mpi::OpGrouped::applyImpl");

   auto comm = _comm.get();
   assert(comm != nullptr);
   if (!comm->isSetup())
   {
      // Set group size
      _groupSize = in.size();
      // Get absolute coordinates
      std::vector<point_t> cooOld =
         View::getCoo<typename Tin::value_type, p012_t>(in[0]);
      std::vector<point_t> cooNew =
         View::getCoo<typename Tout::value_type, Perm>(out[0]);
      assert(cooOld.size() == in[0].size());
      assert(cooNew.size() == out[0].size());
      // Setup
      comm->setComm(cooNew, cooOld, _groupSize);
   }

   assert(in.size() == _groupSize);

   // Collect pointers to data
   constexpr std::size_t maxGroupSize = 16;
   assert(_groupSize <= maxGroupSize);
   // We are using a stack array aka structArray to copy directly to
   // the GPU kernel the pointers, we are no using std::array because
   // of compatibility with CUDA
   structArray<ScalarType*, maxGroupSize> outData;
   structArray<const ScalarType*, maxGroupSize> inData;
   for (std::size_t i = 0; i < _groupSize; ++i)
   {
      outData[i] = out[i].data();
      inData[i] = in[i].data();
   }

   comm->exchange(outData, inData);
}

} // namespace Mpi
} // namespace Transpose
} // namespace QuICC
