/**
 * @file Binary.hpp
 * @brief Operator binary interface
 */
#pragma once

namespace QuICC {
/// @brief This namespace provides base methods and classes for generic operators
namespace Operator {

/** @brief this is the abstract interface class of a generic binary operator.
 *  Its purpose is to allow the usage of a common pointer to akin operators
 *  (i.e. operator with same input/output types) and to keep a consistent API.
 *  @tparam Tout output type, in pratice always a View
 *  @tparam Tin input type, in pratice always a View
 *  @tparam Top operator type, in pratice always a View
 */
template<class Tout, class Tin, class Top>
class BinaryOp
{
public:
    /// @brief apply action of the binary operator
    /// @param out result of the operation
    /// @param in input of the operation with const qualifier
    /// @param op this input represents the operator so it is always const
    virtual void apply(Tout& out, const Tin& in, const Top& op) = 0;

    /// @brief dtor
    virtual ~BinaryOp() = default;
};

/** @brief CRTP base class of a generic binary operator.
 *  @tparam Derived specialized operator class
 *  @tparam Tout output type, in pratice always a View
 *  @tparam Tin input type, in pratice always a View
 */
template <class Derived, class Tout, class Tin, class Top>
class BinaryBaseOp : public BinaryOp<Tout, Tin, Top>
{
public:
    /// @brief apply action of the binary operator by dispatching to derived class implementation
    /// @param out result of the operation
    /// @param in input of the operation with const qualifier
    /// @param op this input represents the operator so it is always const
    void apply(Tout& out, const Tin& in, const Top& op) final {
        auto derivedPtr = static_cast<Derived*>(this);
        derivedPtr->applyImpl(out, in, op);
    }
};


} // namespace Operator
} // namespace QuICC
