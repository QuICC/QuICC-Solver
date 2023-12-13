# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *


class Quicc(CMakePackage, CudaPackage,  ROCmPackage):
    """Quasi-Inverse Convection Code"""

    homepage = "https://github.com/QuICC/QuICC"
    git = "https://github.com/QuICC/QuICC.git"
    #  url = "https://github.com/QuICC/QuICC/archive/refs/tags/v0.8.2.tar.gz"

    maintainers = ["danielganellari"]

    version("develop", branch="gpu_legendre_hip")

    variant("kokkos", default=True)
    variant("openmp", default=True)
    #  variant("python", default=True)
    variant("tests", default=True)
    variant("build_type",
            default="Release",
            description="CMake build type",
            values=("Debug", "Release", "RelWithDebInfo"),
            )

    depends_on("cmake@3.18:", type="build")
    depends_on("mpi")

    depends_on("kokkos")
    depends_on("kokkos+openmp", when="+openmp")

    depends_on("kokkos+cuda+cuda_lambda+wrapper", when="+cuda%gcc")
    depends_on("kokkos+cuda", when="+cuda")

    # rocm dependencies
    depends_on("kokkos+rocm", when="+rocm")

    depends_on("kokkos-kernels")

    depends_on("hdf5")
    depends_on("fftw")
    depends_on("openblas")
    depends_on("boost")
    depends_on("python")
    depends_on("py-numpy")

    def cmake_args(self):
        options = [
            self.define_from_variant("QUICC_TESTSUITE_TRANSFORM", "tests"),
            self.define_from_variant("QUICC_USE_KOKKOS", "kokkos"),
            self.define_from_variant("QUICC_USE_KOKKOS_HIP", "rocm"),
            self.define_from_variant("QUICC_USE_KOKKOS_CUDA", "cuda"),
        ]

        options.append(
            "-DPython_ROOT_DIR=%s" % self.spec["python"].prefix)

        if "+cuda%gcc" in self.spec:
            options.append(
                "-DCMAKE_CXX_COMPILER=%s" % self.spec["kokkos-nvcc-wrapper"].kokkos_cxx
            )
            options.append("-DQUICC_USE_KOKKOS_CUDA=ON")

        if "+cuda" in self.spec:
            cuda_arch = self.spec.variants["cuda_arch"].value
            if cuda_arch[0] != "none":
                options += ["-DCMAKE_CUDA_FLAGS=-arch=sm_{0}".format(cuda_arch[0])]

        if "+rocm" in self.spec:
            options.append(self.define(
                "CMAKE_CXX_COMPILER", self.spec["hip"].hipcc))
            archs = ",".join(self.spec.variants['amdgpu_target'].value)
            options.append("-DHIP_HCC_FLAGS=--amdgpu-target={0}".format(archs))
            options.append("-DCMAKE_CXX_FLAGS=--amdgpu-target={0} --offload-arch={0}".format(archs))
            options.append("-DQUICC_USE_KOKKOS_HIP=ON")

        return options
