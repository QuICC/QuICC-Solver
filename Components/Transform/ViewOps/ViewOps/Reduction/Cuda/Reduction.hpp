/**
 * @file Reduction.hpp
 * @brief Reduction operations on Views
 * Allows for a reduction operation.
 */
#pragma once

// External includes
//
#include <cstdint>
#include <memory>

// Project includes
//
#include "Memory/Memory.hpp"
#include "Memory/MemoryResource.hpp"
#include "Operator/Unary.hpp"
#include "Operator/Nary.hpp"
#include "Environment/Cfl.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {
/// @brief namespace for Reduction type operations
namespace Reduction {
/// @brief namespace for Cuda backends
namespace Cuda {

/// @brief Reduction operator
/// @tparam Tout n dimensional view
/// @tparam Tin  n-1 dimensional view
/// @tparam Dir axis to perform reduction on
template <class Tout, class Tin, std::uint32_t Dir>
class Op : public Operator::UnaryBaseOp<Op<Tout, Tin, Dir>, Tout, Tin>
{
public:
   /// @brief ctor passing memory resource
   /// @param mem memory resource
   Op(std::shared_ptr<Memory::memory_resource> mem);
   /// @brief Default constructor
   Op() = delete;
   /// @brief dtor
   ~Op() = default;

private:
   /// @brief action implementation that does not overwrite the input
   /// @param out differentiatied physical space coefficient
   /// @param in input modes
   void applyImpl(Tout& out, const Tin& in);
   /// @brief action implementation that might modify the input
   /// @param out differentiatied physical space coefficient
   /// @param in input modes
   // void applyImpl(Tout& out, Tin& in);
   /// @brief give access to base class
   friend Operator::UnaryBaseOp<Op<Tout, Tin, Dir>, Tout, Tin>;
   /// @brief memory resource
   /// needs shared ptr for memory pools
   /// note, this must call the dtor last
   /// otherwise we cannot dealloc data
   /// \todo consider removing shared ptr and using singleton
   std::shared_ptr<Memory::memory_resource> _mem;
   /// @brief index typedef
   using IndexType = typename Tin::IndexType;
   /// @brief layer index cache
   Memory::MemBlock<IndexType> _layerIndex;
   /// @brief layer width cache
   Memory::MemBlock<IndexType> _layerWidth;
   /// @brief max layer width cache
   std::uint32_t _N;
   /// @brief input offset cache
   Memory::MemBlock<IndexType> _offSetIn;
   /// @brief output offset cache
   Memory::MemBlock<IndexType> _offSetOut;
};
/// @tparam T scalar
template <class T = double> struct MagVelFunctor
{
   /// @brief non dimensional scaling for transport term
   T _scaling;

   /// @brief ctor
   /// @param scaling
   MagVelFunctor(T scaling) : _scaling(scaling){};

   /// @brief deleted default constructor
   MagVelFunctor() = delete;

   /// @brief dtor
   ~MagVelFunctor() = default;

   /// @brief Cross product, component wise
   /// @param uj
   /// @param uk
   /// @param vj
   /// @param vk
   /// @return i component of cross product
   __host__ __device__ T operator()(T ui, T uj, T uk, T vi, T vj, T vk)
   {
      return _scaling * (ui * vi + uj * vj + uk * vk);
   }
};
template <class Functor, class Tout, class Tin0, class Tin1, class Tin2, class Tin3, class Tin4, class Tin5>
class OpCfl : public QuICC::Operator::NaryBaseOp<OpCfl<Functor, Tout, Tin0, Tin1, Tin2, Tin3, Tin4, Tin5>, Tout, Tin0, Tin1, Tin2, Tin3, Tin4, Tin5>
{
private:
   /// @brief stored functor, i.e. struct with method
   /// Tout::ScalarType operator()(Targs::ScalarType var, ...)
	Functor _f;

public:
   /// @brief capture functor by value
   /// @param f functor, i.e. struct with method
   /// Tout::ScalarType operator()(Targs::ScalarType var, ...)
   OpCfl(Functor f) : _f(f){};
   /// @brief default constructor
   OpCfl() = delete;
   /// @brief dtor
   ~OpCfl() = default;

private:
   /// @brief action implementation
   /// @param out output View
   /// @param ...args input Views
   void applyImpl(Tout& out, const Tin0&, const Tin1&, const Tin2&, const Tin3&, const Tin4&, const Tin5&);
   /// @brief give access to base class
   friend QuICC::Operator::NaryBaseOp<OpCfl<Functor, Tout, Tin0, Tin1, Tin2, Tin3, Tin4, Tin5>, Tout, Tin0, Tin1, Tin2, Tin3, Tin4, Tin5>;
};

} // namespace Cuda
} // namespace Reduction
} // namespace QuICC
