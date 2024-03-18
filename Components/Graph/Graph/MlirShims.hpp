#pragma once


#include <Quiccir-c/Utils.h>
#include <complex>

#include "ViewOps/ViewMemoryUtils.hpp"

namespace QuICC
{
namespace Graph
{
    using R_DCCSC3D_t = View::View<double, View::DCCSC3D>;
    using C_DCCSC3D_t = QuICC::View::View<std::complex<double>, View::DCCSC3D>;
    using view3_cd_t = ViewDescriptor<std::complex<double>, std::uint32_t, 3>;

} // namespace Graph
} // namespace QuICC


