/**
 * @file FftTypes.hpp
 * @brief Types supported by FFT backends
 */
#pragma once

// External includes
//

// Project includes
//
#include "View/View.hpp"

namespace QuICC {
namespace Fft {

/// @brief Complex dense 2D tensor, input modes view type
using CmodsDense2D_t = View::View<std::complex<double>, View::dense2D>;
/// @brief Complex dense 2D tensor, output phys view type
using CphysDense2D_t = View::View<std::complex<double>, View::dense2D>;
/// @brief Real dense 2D tensor, output phys view type
using RphysDense2D_t = View::View<double, View::dense2D>;
/// @brief Complex compressed sparse layer 3D tensor, input modes view type
using CmodsDCCSC3D_t = View::View<std::complex<double>, View::DCCSC3D>;
/// @brief Real compressed sparse layer 3D tensor, output phys view type
using RphysDCCSC3D_t = View::View<double, View::DCCSC3D>;
/// @brief Complex compressed sparse layer 3D tensor, output phys view type
using CphysDCCSC3DInOrder_t = View::View<std::complex<double>, View::DCCSC3DInOrder>;
/// @brief Complex compressed sparse layer 3D tensor, input mods view type
using CmodsDCCSC3DInOrder_t = View::View<std::complex<double>, View::DCCSC3DInOrder>;

} // namespace Fft
} // namespace QuICC
