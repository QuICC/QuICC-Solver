/**
 * @file FunctorHelpers.hpp
 * @brief Helpers for working with functors
 */

#ifndef QUICC_PHYSICAL_DETAILS_FUNCTORHELPERS_HPP
#define QUICC_PHYSICAL_DETAILS_FUNCTORHELPERS_HPP

// System includes
//
#include <cstdint>
#include <tuple>
#include <type_traits>

// Project includes
//
#include "Std/Cuda/Utility.hpp"

namespace QuICC {

namespace Physical {

namespace details {

template <bool Enable, typename T, typename... Tcache>
   struct CachedFunctor
   {
      /// Enable caching
      static constexpr bool enableCaching(){return Enable;};

      /// Cache variables
      std::tuple<Tcache...> _s;

      /// @brief ctor
      CachedFunctor(): _s(Tcache{0}...) {};

      /// @brief dtor
      ~CachedFunctor() = default;
   };

template <typename T, template <typename> class TFunctor>
struct AddTmplFunctor : public TFunctor<T>
{
   /// @brief ctor
   /// @param scaling
   template <typename... Ts>
   AddTmplFunctor(Ts... s) : TFunctor<T>(s...) {};

   /// @brief deleted default constructor
   AddTmplFunctor() = delete;

   /// @brief dtor
   ~AddTmplFunctor() = default;

   /// @brief add inherited operator
   template <typename... Args> QUICC_CUDA_HOSTDEV T operator()(Args... args)
   {
      return call_with_tuple(args_to_tuple(args...),
         std::make_integer_sequence<std::size_t, sizeof...(Args) - 1>());
   }

   private:
      using TFunctor<T>::operator();

      /// @brief Create tuple of arguments
      template <typename... Args> decltype(auto) args_to_tuple(Args... args)
      {
         return std::tuple<Args...>(args...);
      }

      /// @brief Call functor tuple of arguments
      template <std::size_t... Is, typename... Args>
      T call_with_tuple(const std::tuple<Args...>& tuple,
         std::index_sequence<Is...>)
      {
         return std::get<sizeof...(Args) - 1>(tuple) +
                TFunctor<T>::operator()(std::get<Is>(tuple)...);
      }
};

template <typename T, template <typename> class TFunctor>
struct SubTmplFunctor : public TFunctor<T>
{
   /// @brief ctor
   /// @param scaling
   template <typename... Ts>
   SubTmplFunctor(Ts... s) : TFunctor<T>(-s...) {};

   /// @brief deleted default constructor
   SubTmplFunctor() = delete;

   /// @brief dtor
   ~SubTmplFunctor() = default;

   /// @brief add inherited operator
   template <typename... Args> QUICC_CUDA_HOSTDEV T operator()(Args... args)
   {
      return call_with_tuple(args_to_tuple(args...),
         std::make_integer_sequence<std::size_t, sizeof...(Args) - 1>());
   }

   private:
      using TFunctor<T>::operator();

      /// @brief Create tuple of arguments
      template <typename... Args> decltype(auto) args_to_tuple(Args... args)
      {
         return std::tuple<Args...>(args...);
      }

      /// @brief Call functor tuple of arguments
      template <std::size_t... Is, typename... Args>
      T call_with_tuple(const std::tuple<Args...>& tuple,
         std::index_sequence<Is...>)
      {
         return std::get<sizeof...(Args) - 1>(tuple) +
                TFunctor<T>::operator()(std::get<Is>(tuple)...);
      }
};

} // namespace details
} // namespace Physical
} // namespace QuICC
#endif // QUICC_PHYSICAL_DETAILS_FUNCTORHELPERS_HPP
