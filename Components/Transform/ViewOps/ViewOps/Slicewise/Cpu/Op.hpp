/**
 * @file Pointwise.hpp
 * @brief Pointwise operations on Views
 * Allows for any user defined pointwise operation.
 * The operation is defined via a functor object.
 * Value semantic lets a (good) compiler easily inline and
 * remove the indirect call.
 */
#pragma once

// External includes
//
#include <type_traits>

// Project includes
//
#include "Operator/Nary.hpp"
#include "Profiler/Interface.hpp"
#include "View/Attributes.hpp"

namespace QuICC {
/// @brief namespace for Pointwise type operations
namespace Pointwise {
/// @brief namespace for cpu backends
namespace Cpu {

using namespace QuICC::Operator;

/// @brief Pointwise operator
/// @tparam Functor Nary scalar functor
/// @tparam Tout output View
/// @tparam ...Targs input Views
template <class Functor, class Tout, class... Targs>
class Op : public NaryBaseOp<Op<Functor, Tout, Targs...>, Tout, Targs...>
{
private:
   /// @brief stored functor, i.e. struct with method
   /// Tout::ScalarType operator()(Targs::ScalarType var, ...)
   Functor _f;

public:
   /// @brief capture functor by value
   /// @param f functor, i.e. struct with method
   /// Tout::ScalarType operator()(Targs::ScalarType var, ...)
   Op(Functor f) : _f(f){};
   /// @brief default constructor
   Op() = delete;
   /// @brief dtor
   ~Op() = default;

private:
   /// @brief action implementation
   /// @param out output View
   /// @param ...args input Views
   void applyImpl(Tout& out, const Targs&... args);
   /// @brief give access to base class
   friend NaryBaseOp<Op<Functor, Tout, Targs...>, Tout, Targs...>;
   /// @brief index typedef
   using IndexType = typename Tout::IndexType;
   /// @brief layer index cache
   std::vector<IndexType> _layerIndex;
   /// @brief layer width cache
   std::vector<IndexType> _layerWidth;
};

template <class Functor, class Tout, class ...Targs>
void Op<Functor, Tout, Targs...>::applyImpl(Tout& out, const Targs&... args)
{
   Profiler::RegionFixture<4> fix("Pointwise::Cpu::applyImpl");

   // check Tout and Targs.. match Functor op
   using res_t = std::invoke_result_t<Functor, typename Targs::ScalarType...>;
   static_assert(std::is_same_v<typename Tout::ScalarType, res_t>,
      "Mismatch in functor or arguments");
   // check same size
   assert(((out.size() == args.size()) && ... ));

   // implemented only for physical space
   static_assert(std::is_same_v<Tout::Attributes, View::DCCSC3D>());

   // cache populated layers
   auto pointers = out.pointers();
   if (_layerIndex.size() < 1)
   {
      for (IndexType k = 0; k < pointers.size() - 1; ++k)
      {
         IndexType nCols = pointers[k + 1] - pointers[k];
         assert(nCols <= out.dims()[1]);
         // check if layer is populated
         if (nCols > 0)
         {
            _layerIndex.push_back(k);
            _layerWidth.push_back(nCols);
         }
      }
   }

   // apply slicewise functor
   std::size_t offSet = 0;
   for (IndexType h = 0; h < _layerIndex.size(); ++h)
   {
      // layer index
      auto l = _layerIndex[h];

      // DCCSC3D
      // get slice dimension
      auto M = out.lds();
      auto N = _layerWidth[h];

      // check mem bounds
      assert(offSet + M * N <= out.size());

      // column major
      for (std::size_t n = 0; n < N; ++n)
      {
         for (std::size_t m = 0; m < M; ++m)
         {
            auto mnk = offSet + m + n*M;
            out[mnk] = _f(l, args[mnk]...);
         }
      }

      // update offset
      offSet += M * N;
   }
}

} // namespace Cpu
} // namespace Pointwise
} // namespace QuICC
