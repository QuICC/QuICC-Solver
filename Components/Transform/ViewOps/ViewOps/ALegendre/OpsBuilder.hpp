/**
 * @file OpsBuilder.hpp
 * @brief Mapping Generic ALegendre operator builders to specific Ops
 */

#pragma once

// External includes
//

// Project includes
//
#include "QuICC/Polynomial/ALegendre/Plm.hpp"
#include "QuICC/Polynomial/ALegendre/dPlm.hpp"
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"
#include "ViewOps/ALegendre/Builder.hpp"
#include "ViewOps/ALegendre/Tags.hpp"

namespace QuICC {
namespace Transform {
namespace ALegendre {


/// @brief This namespace hides implementation details
namespace details {

/// @brief Generic mapping
/// @tparam VOP operator view type
/// @tparam TAG kind
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class TAG, class DIR>
struct OpsBuilderMap {
    using type = void;
};

} // namespace details

/// @brief Convenience wrapper
/// (you cannot specialize type aliases)
/// @tparam VOP operator view type
/// @tparam TAG kind
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class TAG, class DIR>
using OpsBuilder = typename details::OpsBuilderMap<VOP, TAG, DIR>::type;


namespace details {

using namespace QuICC::Polynomial::ALegendre;

/// convenience wrapper for builder function
/// \todo cleanup underlaying builder to match JW
template<class Tview, class TPolyBuilder, class Tdata, std::int16_t LlDiff>
struct HelperBuilder
{
    void compute(Tview opView, const Internal::Array& igrid,
   const Internal::Array& iweights)
   {
        builder<Tview, TPolyBuilder, Tdata, LlDiff>(opView, igrid, iweights);
   }
};

/// @brief P Builder
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR>
struct OpsBuilderMap<VOP, P_t, DIR> {
    using type = HelperBuilder<VOP, ::QuICC::Polynomial::ALegendre::Plm, ::QuICC::Internal::Array::Scalar, 0>;
};

/// @brief D1 Builder
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR>
struct OpsBuilderMap<VOP, D1_t, DIR> {
    using type = HelperBuilder<VOP, ::QuICC::Polynomial::ALegendre::dPlm, ::QuICC::Internal::Array::Scalar, 0>;
};

/// @brief Ll Builder
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR>
struct OpsBuilderMap<VOP, Ll_t, DIR> {
    using type = HelperBuilder<VOP, ::QuICC::Polynomial::ALegendre::Plm, ::QuICC::Internal::Array::Scalar, 1>;
};

/// @brief LlD1 Builder
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR>
struct OpsBuilderMap<VOP, LlD1_t, DIR> {
    using type = HelperBuilder<VOP, ::QuICC::Polynomial::ALegendre::dPlm, ::QuICC::Internal::Array::Scalar, 1>;
};

/// @brief DivS1 Builder
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR>
struct OpsBuilderMap<VOP, DivS1_t, DIR> {
    using type = HelperBuilder<VOP, ::QuICC::Polynomial::ALegendre::sin_1Plm, ::QuICC::Internal::Array::Scalar, 0>;
};

/// @brief DivS1Dp Builder
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR>
struct OpsBuilderMap<VOP, DivS1Dp_t, DIR> {
    using type = OpsBuilder<VOP, DivS1_t, DIR>;
};

/// @brief LlDivS1 Builder
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR>
struct OpsBuilderMap<VOP, LlDivS1_t, DIR> {
    using type = HelperBuilder<VOP, ::QuICC::Polynomial::ALegendre::sin_1Plm, ::QuICC::Internal::Array::Scalar, 1>;
};

/// @brief LlDivS1Dp Builder
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR>
struct OpsBuilderMap<VOP, LlDivS1Dp_t, DIR> {
    using type = OpsBuilder<VOP, LlDivS1_t, DIR>;
};

/// @brief Ll2 Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP>
struct OpsBuilderMap<VOP, Ll2_t, fwd_t> {
    using type = HelperBuilder<VOP, ::QuICC::Polynomial::ALegendre::Plm, ::QuICC::Internal::Array::Scalar, 2>;
};

/// @brief DivLl Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP>
struct OpsBuilderMap<VOP, DivLl_t, fwd_t> {
    using type = HelperBuilder<VOP, ::QuICC::Polynomial::ALegendre::Plm, ::QuICC::Internal::Array::Scalar, -1>;
};

/// @brief DivLlD1 Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP>
struct OpsBuilderMap<VOP, DivLlD1_t, fwd_t> {
    using type = HelperBuilder<VOP, ::QuICC::Polynomial::ALegendre::dPlm, ::QuICC::Internal::Array::Scalar, -1>;
};

/// @brief DivLlDivS1 Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP>
struct OpsBuilderMap<VOP, DivLlDivS1_t, fwd_t> {
    using type = HelperBuilder<VOP, ::QuICC::Polynomial::ALegendre::sin_1Plm, ::QuICC::Internal::Array::Scalar, -1>;
};

/// @brief DivLlDivS1Dp Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP>
struct OpsBuilderMap<VOP, DivLlDivS1Dp_t, fwd_t> {
    using type = OpsBuilder<VOP, DivLlDivS1_t, fwd_t>;
};

} // namespace details

} // namespace ALegendre
} // namespace Transform
} // namespace QuICC
