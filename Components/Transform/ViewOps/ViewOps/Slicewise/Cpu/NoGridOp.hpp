/**
 * @file NoGridOp.hpp
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
#include <utility>

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
/// @tparam Ng1 number of grid data 1, meaning depends on DIR
/// @tparam Ng2 number of grid data 2, meaning depends on DIR
/// @tparam Ng3 number of grid data 3, meaning depends on DIR
/// @tparam ...Targs input Views
template <std::uint8_t Dir, class Functor, class Tout, std::uint8_t Ng1,
   std::uint8_t Ng2, std::uint8_t Ng3, class... Targs>
class NoGridOp : public NaryBaseOp<NoGridOp<Dir, Functor, Tout, Ng1, Ng2, Ng3, Targs...>,
              Tout, Targs...>
{
private:
   /// @brief stored functor, i.e. struct with method
   /// Tout::ScalarType operator()(Targs::ScalarType var, ...)
   Functor _f;

public:
   /// @brief capture functor by value
   /// @param f functor, i.e. struct with method
   /// Tout::ScalarType operator()(Targs::ScalarType var, ...)
   NoGridOp(Functor f) : _f(f){};
   /// @brief default constructor
   NoGridOp() = delete;
   /// @brief dtor
   ~NoGridOp() = default;

private:
   /// @brief give access to base class
   friend NaryBaseOp<NoGridOp<Dir, Functor, Tout, Ng1, Ng2, Ng3, Targs...>, Tout,
      Targs...>;
   /// @brief action implementation
   void applyImpl(Tout& out, const Targs&... args);
   /// @brief specialized implementation for Phi-Theta slice
   template <std::size_t... Is, std::size_t... Js>
   void phiThetaImpl(Tout& out, const std::tuple<Targs...>& args, std::index_sequence<Is...>, std::index_sequence<Js...>);
   /// @brief specialized implementation for Phi-R slice
   template <std::size_t... Is, std::size_t... Js>
   void phiRImpl(Tout& out, const std::tuple<Targs...>& args, std::index_sequence<Is...>, std::index_sequence<Js...>);
   /// @brief specialized implementation for Phi slice
   template <std::size_t... Is, std::size_t... Js, std::size_t... Ks>
   void phiImpl(Tout& out, const std::tuple<const Targs&...>& args, std::index_sequence<Is...>, std::index_sequence<Js...>, std::index_sequence<Ks...>);
   /// @brief index typedef
   using IndexType = typename Tout::IndexType;
   /// @brief layer index cache
   std::vector<IndexType> _layerIndex;
   /// @brief layer width cache
   std::vector<IndexType> _layerWidth;
};


template <std::uint8_t Dir, class Functor, class Tout,
   std::uint8_t Ng1, std::uint8_t Ng2, std::uint8_t Ng3, class... Targs>
void NoGridOp<Dir, Functor, Tout, Ng1, Ng2, Ng3, Targs...>::applyImpl(Tout& out,
   const Targs&... args)
{
   Profiler::RegionFixture<4> fix("SlicewiseNoGrid::Cpu::applyImpl");

   // check Tout and Targs.. match Functor op
   using res_t =
      std::invoke_result_t<Functor, typename Targs::ScalarType...>;
   static_assert(std::is_same_v<typename Tout::ScalarType, res_t>,
      "Mismatch in functor or arguments");

   // implemented only for physical space
   static_assert(std::is_same_v<Tout, View::View<double, View::DCCSC3D>>);

   if constexpr (Dir == 1)
   {
      constexpr std::size_t n = sizeof...(Targs) - Ng1;
      phiRImpl(out,
            std::tuple<Targs...>(args...),
            std::make_integer_sequence<std::size_t, Ng1>(),
            std::make_integer_sequence<std::size_t, n>());
   }
   else if constexpr (Dir == 2)
   {
      constexpr std::size_t n = sizeof...(Targs) - Ng1;
      phiThetaImpl(out,
            std::tuple<Targs...>(args...),
            std::make_integer_sequence<std::size_t, Ng1>(),
            std::make_integer_sequence<std::size_t, n>());

   }
   else if constexpr (Dir == 3)
   {
      constexpr std::size_t n = sizeof...(Targs) - Ng1 - Ng2;
      phiImpl(out,
            std::forward_as_tuple(args...),
            std::make_integer_sequence<std::size_t, Ng1>(),
            std::make_integer_sequence<std::size_t, Ng2>(),
            std::make_integer_sequence<std::size_t, n>());

   }
   else
   {
      throw std::logic_error("This slice direction is not implemented.");
   }
}

template <std::uint8_t Dir, class Functor, class Tout,
   std::uint8_t Ng1, std::uint8_t Ng2, std::uint8_t Ng3, class... Targs>
   template <std::size_t... Is, std::size_t... Js>
void NoGridOp<Dir, Functor, Tout, Ng1, Ng2, Ng3, Targs...>::phiRImpl(Tout& out,
      const std::tuple<Targs...>& args, std::index_sequence<Is...>, std::index_sequence<Js...>)
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
      assert(thetaIdx < std::get<0>(args).size());

      // column major
      for (std::size_t m = 0; m < M; ++m)
      {
         auto mnk = m + col * M;
         out[mnk] = _f(std::get<Is>(args)[thetaIdx]..., std::get<Js+Ng1>(args)[mnk]...);
      }
   }
}

template <std::uint8_t Dir, class Functor, class Tout,
   std::uint8_t Ng1, std::uint8_t Ng2, std::uint8_t Ng3, class... Targs>
   template <std::size_t... Is, std::size_t... Js>
void NoGridOp<Dir, Functor, Tout, Ng1, Ng2, Ng3, Targs...>::phiThetaImpl(Tout& out,
   const std::tuple<Targs...>& args, std::index_sequence<Is...>, std::index_sequence<Js...>)
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
            out[mnk] = _f(std::get<Is>(args)[l]..., std::get<Js+Ng1>(args)[mnk]...);
         }
      }

      // update offset
      offSet += M * N;
   }
}

template <std::uint8_t Dir, class Functor, class Tout,
   std::uint8_t Ng1, std::uint8_t Ng2, std::uint8_t Ng3, class... Targs>
   template <std::size_t... Is, std::size_t... Js, std::size_t... Ks>
void NoGridOp<Dir, Functor, Tout, Ng1, Ng2, Ng3, Targs...>::phiImpl(Tout& out,
   const std::tuple<const Targs&...>& args, std::index_sequence<Is...>, std::index_sequence<Js...>, std::index_sequence<Ks...>)
{
   assert(Dir == 3);

   assert(Ng1 > 0);
   assert(Ng2 > 0);

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
   auto indices = out.indices()[1];
   std::size_t offSet = 0;
   std::size_t colIdx = 0;
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

      auto mnk = offSet;
      for (std::size_t n = 0; n < N; ++n)
      {
         auto thetaIdx = indices[colIdx];
         for (std::size_t m = 0; m < M; ++m)
         {
            out[mnk] = _f(std::get<Is>(args)[l]..., std::get<Js+Ng1>(args)[thetaIdx]..., std::get<Ks+Ng1+Ng2>(args)[mnk]...);
            mnk++;
         }
         colIdx++;
      }

      // update offset
      offSet += M * N;
   }
}

} // namespace Cpu
} // namespace Slicewise
} // namespace QuICC
