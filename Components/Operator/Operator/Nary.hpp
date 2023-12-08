/**
 * @file Nary.hpp
 * @brief Operator Nary interface
 */
#pragma once

namespace QuICC {
/// @brief This namespace provides base methods and classes for generic
/// operators
namespace Operator {

/** @brief this is the abstract interface class of a generic Nary operator.
 *  Its purpose is to allow the usage of a common pointer to akin operators
 *  (i.e. operator with same input/output types) and to keep a consistent API.
 *  @tparam Tout output type, in pratice always a View
 *  @tparam Tin input type, in pratice always a View
 *  @tparam Top operator type, in pratice always a View
 */
template <class Tout, class... Targs> class NaryOp
{
public:
   /// @brief apply action of the Nary operator
   /// @param out result of the operation
   virtual void apply(Tout&, const Targs&...) = 0;

   /// @brief dtor
   virtual ~NaryOp() = default;
};

/** @brief CRTP base class of a generic Nary operator.
 *  @tparam Derived specialized operator class
 */
template <class Derived, class Tout, class... Targs>
class NaryBaseOp : public NaryOp<Tout, Targs...>
{
public:
   /// @brief apply action of the Nary operator by dispatching to derived class
   /// implementation
   void apply(Tout& out, const Targs&... args) final
   {
      auto derivedPtr = static_cast<Derived*>(this);
      derivedPtr->applyImpl(out, args...);
   }
};


} // namespace Operator
} // namespace QuICC
