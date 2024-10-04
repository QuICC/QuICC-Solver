/**
 * @file Types.hpp
 * @brief
 */
#pragma once

// External includes
//
#include <variant>
#include <memory>
#ifdef QUICC_HAS_CUDA_BACKEND
#include <cuda/std/complex>
#endif


// Project includes
//
#include "Operator/Nary.hpp"
#include "Operator/Unary.hpp"
#include "View/View.hpp"

namespace QuICC {
namespace Graph {

using QuICC::Operator::NaryOp;
using QuICC::Operator::UnaryOp;

/// @brief Fourier Phys
using R_DCCSC3D_t = View::View<double, View::DCCSC3D>;
/// @brief Fourier Mods or AL Phys or JW Phys/Mods
using C_DCCSC3D_t = QuICC::View::View<std::complex<double>, View::DCCSC3D>;
using C_DCCSC3DJIK_t = QuICC::View::View<std::complex<double>, View::DCCSC3DJIK>;
#ifdef QUICC_HAS_CUDA_BACKEND
using Ccuda_DCCSC3DJIK_t = QuICC::View::View<cuda::std::complex<double>, View::DCCSC3DJIK>;
#endif
/// @brief AL Mods
using C_S1CLCSC3D_t = QuICC::View::View<std::complex<double>, View::S1CLCSC3D>;
using C_S1CLCSC3DJIK_t = QuICC::View::View<std::complex<double>, View::S1CLCSC3DJIK>;

/// @brief encode data types
using varData_t = std::variant<
  R_DCCSC3D_t,
  C_DCCSC3D_t
  // C_S1CLCSC3D_t
>;

/// @brief encode shared pointers to op type
using varOp_t = std::variant<
    std::shared_ptr<NaryOp<C_DCCSC3D_t, C_DCCSC3D_t>>,
    std::shared_ptr<NaryOp<R_DCCSC3D_t, R_DCCSC3D_t>>,
    std::shared_ptr<NaryOp<C_DCCSC3D_t, C_DCCSC3D_t, C_DCCSC3D_t>>,
    #ifdef QUICC_HAS_CUDA_BACKEND
    std::shared_ptr<NaryOp<Ccuda_DCCSC3DJIK_t, Ccuda_DCCSC3DJIK_t, Ccuda_DCCSC3DJIK_t>>,
    #endif
    std::shared_ptr<NaryOp<R_DCCSC3D_t, R_DCCSC3D_t, R_DCCSC3D_t>>,
    std::shared_ptr<NaryOp<R_DCCSC3D_t, R_DCCSC3D_t, R_DCCSC3D_t, R_DCCSC3D_t, R_DCCSC3D_t>>,
    std::shared_ptr<NaryOp<
      R_DCCSC3D_t, R_DCCSC3D_t, R_DCCSC3D_t, R_DCCSC3D_t,
      R_DCCSC3D_t, R_DCCSC3D_t, R_DCCSC3D_t>>,
    std::shared_ptr<UnaryOp<C_DCCSC3D_t, C_DCCSC3D_t>>,
    std::shared_ptr<UnaryOp<C_DCCSC3D_t, C_DCCSC3DJIK_t>>,
    std::shared_ptr<UnaryOp<C_DCCSC3DJIK_t, C_DCCSC3DJIK_t>>,
    std::shared_ptr<UnaryOp<C_DCCSC3DJIK_t, C_DCCSC3D_t>>,
    std::shared_ptr<UnaryOp<R_DCCSC3D_t, C_DCCSC3D_t>>,
    std::shared_ptr<UnaryOp<C_DCCSC3D_t, R_DCCSC3D_t>>,
    std::shared_ptr<UnaryOp<R_DCCSC3D_t, R_DCCSC3D_t>>,
    std::shared_ptr<UnaryOp<C_S1CLCSC3D_t, C_DCCSC3D_t>>,
    std::shared_ptr<UnaryOp<C_S1CLCSC3DJIK_t, C_DCCSC3DJIK_t>>,
    std::shared_ptr<UnaryOp<C_DCCSC3D_t, C_S1CLCSC3D_t>>,
    std::shared_ptr<UnaryOp<C_DCCSC3DJIK_t, C_S1CLCSC3DJIK_t>>
>;


template<typename T, std::size_t N>
struct MemRefDescriptor {
  T *allocated;
  T *aligned;
  intptr_t offset;
  intptr_t sizes[N];
  intptr_t strides[N];
};

} // namespace Graph
} // namespace QuICC
