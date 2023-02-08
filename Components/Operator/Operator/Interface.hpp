/**
 * @file Interface.hpp
 * @brief Operator interface
 */
#pragma once

namespace QuICC {
/// @brief This namespace provides base methods and classes for generic operators
namespace Operator {

/** @brief this is the abstract interface class of a generic unary operator.
 *  Its purpose is to allow the usage of a common pointer to akin operators
 *  (i.e. operator with same input/output types) and to keep a consistent API.
 *  @tparam Tout output type, in pratice always a View
 *  @tparam Tin input type, in pratice always a View
 */
template<class Tout, class Tin>
class InterfaceOp
{
public:
    /// @brief apply action of the unary operator
    /// @param out result of the operation
    /// @param in input of the operation with const qualifier
    virtual void apply(Tout& out, const Tin& in) = 0;

    /// @brief apply action of the unary operator
    /// @param out result of the operation
    /// @param in input of the operation
    virtual void apply(Tout& out, Tin& in) = 0;

    /// @brief dtor
    virtual ~InterfaceOp() = default;
};

/** @brief CRTP base class of a generic unary operator.
 *  @tparam Derived specialized operator class
 *  @tparam Tout output type, in pratice always a View
 *  @tparam Tin input type, in pratice always a View
 */
template <class Derived, class Tout, class Tin>
class BaseOp : public InterfaceOp<Tout, Tin>
{
public:
    /// @brief apply action of the unary operator by dispatching to derived class implementation
    /// @param out result of the operation
    /// @param in input of the operation
    void apply(Tout& out, Tin& in) final {
        auto derivedPtr = static_cast<Derived*>(this);
        derivedPtr->applyImpl(out, in);
    }

    /// @brief apply action of the unary operator by dispatching to derived class implementation
    /// @param out result of the operation
    /// @param in input of the operation with const qualifier
    void apply(Tout& out, const Tin& in) final {
        auto derivedPtr = static_cast<Derived*>(this);
        derivedPtr->applyImpl(out, in);
    }
};


} // namespace Operator
} // namespace QuICC
