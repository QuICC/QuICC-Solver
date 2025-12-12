/**
 * @file NoTwoGridOp.hpp
 * @brief Slicewise operations on Views without storing grid
 * Allows for any user defined Slicewise operation with grid dependence.
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
#include "Memory/Memory.hpp"
#include "Memory/MemoryResource.hpp"
#include "Operator/Nary.hpp"
#include "Profiler/Interface.hpp"
#include "Types/Internal/Casts.hpp"
#include "View/Attributes.hpp"
#include "View/View.hpp"


namespace QuICC {
/// @brief namespace for Slicewise type operations
namespace Slicewise {
/// @brief namespace for cpu backends
namespace Cpu {

using namespace QuICC::Operator;

/// @brief Slicewise operator
/// @tparam Functor Nary scalar functor
/// @tparam Tout output View
/// @tparam ...Targs input Views
template <std::uint8_t Dir, class Functor, class Tout, class Tgrid1, class Tgrid2,
   class... Targs>
class NoTwoGridOp : public NaryBaseOp<NoTwoGridOp<Dir, Functor, Tout, Tgrid1, Tgrid2, Targs...>,
              Tout, Tgrid1, Tgrid2, Targs...>
{
private:
   /// @brief stored functor, i.e. struct with method
   /// Tout::ScalarType operator()(Targs::ScalarType var, ...)
   Functor _f;

public:
   /// @brief capture functor by value
   /// @param f functor, i.e. struct with method
   /// Tout::ScalarType operator()(Targs::ScalarType var, ...)
   NoTwoGridOp(Functor f) : _f(f){};
   /// @brief default constructor
   NoTwoGridOp() = delete;
   /// @brief dtor
   ~NoTwoGridOp() = default;

private:
   /// @brief give access to base class
   friend NaryBaseOp<NoTwoGridOp<Dir, Functor, Tout, Tgrid1, Tgrid2, Targs...>, Tout,
      Tgrid1, Tgrid2, Targs...>;
   /// @brief action implementation
   /// @param out output View
   /// @param ...args input Views
   void applyImpl(Tout& out, const Tgrid1&, const Tgrid2&, const Targs&... args);
   /// @brief specialized implementation for Phi-Theta slice
   void phiThetaImpl(Tout& out, const Tgrid1& grid1, const Tgrid2& grid2, const Targs&... args);
   /// @brief specialized implementation for Phi-R slice
   void phiRImpl(Tout& out, const Tgrid1& grid1, const Tgrid2& grid2, const Targs&... args);
   /// @brief index typedef
   using IndexType = typename Tout::IndexType;
   /// @brief layer index cache
   std::vector<IndexType> _layerIndex;
   /// @brief layer width cache
   std::vector<IndexType> _layerWidth;
};


template <std::uint8_t Dir, class Functor, class Tout,
   class Tgrid1, class Tgrid2, class... Targs>
void NoTwoGridOp<Dir, Functor, Tout, Tgrid1, Tgrid2, Targs...>::applyImpl(Tout& out,
   const Tgrid1& grid1, const Tgrid2& grid2, const Targs&... args)
{
   Profiler::RegionFixture<4> fix("SlicewiseNoTwoGrid::Cpu::applyImpl");

   // check Tout and Targs.. match Functor op
   using res_t =
      std::invoke_result_t<Functor, typename Tgrid1::ScalarType, typename Tgrid2::ScalarType, typename Targs::ScalarType...>;
   static_assert(std::is_same_v<typename Tout::ScalarType, res_t>,
      "Mismatch in functor or arguments");
   // check same size
   assert(((out.size() == args.size()) && ...));

   // implemented only for physical space
   static_assert(std::is_same_v<Tout, View::View<double, View::DCCSC3D>>);

   if constexpr (Dir == 1)
   {
      phiRImpl(out, grid1, grid2, args...);
   }
   else if constexpr (Dir == 2)
   {
      phiThetaImpl(out, grid1, grid2, args...);
   }
   else
   {
      throw std::logic_error("This slice direction is not implemented.");
   }
}

template <std::uint8_t Dir, class Functor, class Tout,
   class Tgrid1, class Tgrid2, class... Targs>
void NoTwoGridOp<Dir, Functor, Tout, Tgrid1, Tgrid2, Targs...>::phiRImpl(Tout& out,
   const Tgrid1& grid1, const Tgrid2& grid2, const Targs&... args)
{
   assert(Dir == 1);

   // apply slicewise functor
   auto indices = out.indices()[1];
   // DCCSC3D
   // column height
   auto M = out.lds();
   for (std::size_t col = 0; col < indices.size(); ++col)
   {
      // column Id index or Theta Idx
      auto thetaIdx = indices[col];

      // check mem bounds
      assert((col + 1) * M <= out.size());
      assert(thetaIdx < grid1.size());
      assert(thetaIdx < grid2.size());

      // column major
      for (std::size_t m = 0; m < M; ++m)
      {
         auto mnk = m + col * M;
         out[mnk] = _f(grid1[thetaIdx], grid2[thetaIdx], args[mnk]...);
      }
   }
}


template <std::uint8_t Dir, class Functor, class Tout,
   class Tgrid1, class Tgrid2, class... Targs>
void NoTwoGridOp<Dir, Functor, Tout, Tgrid1, Tgrid2, Targs...>::phiThetaImpl(Tout& out,
   const Tgrid1& grid1, const Tgrid2& grid2, const Targs&... args)
{
   assert(Dir == 2);

   // cache populated layers
   auto pointers = out.pointers()[1];
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
            auto mnk = offSet + m + n * M;
            out[mnk] = _f(grid1[l], grid2[l], args[mnk]...);
         }
      }

      // update offset
      offSet += M * N;
   }
}

} // namespace Cpu
} // namespace Slicewise
} // namespace QuICC
