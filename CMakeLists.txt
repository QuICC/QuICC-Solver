#
# Project information
#
cmake_minimum_required(VERSION 3.16.3)


# Add release with assert
include(cmake.d/AddRelWithAssert.cmake)
# Add release for Valgrind
include(cmake.d/AddRelForValgrind.cmake)
# Set default build type
set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
options are: Debug Release RelWithDebInfo RelWithAssert RelForValgrind MinSizeRel.")
# Toggle-able build type
set_property(CACHE CMAKE_BUILD_TYPE
  PROPERTY STRINGS "Debug;Release;RelWithDebInfo;RelWithAssert;RelForValgrind;MinSizeRel")

project(QuICC CXX C)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set CMake module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake.d")

include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
  enable_language(CUDA)
else()
  message(STATUS "No CUDA support detected")
endif()


include(MacroOutOfSourceBuild)
MACRO_ENSURE_OUT_OF_SOURCE_BUILD(
   "${CMAKE_PROJECT_NAME} requires an out of source build. Delete CMakeCache.txt and CMakeFiles/"
)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS OFF) # Use standard C++
# Set CUDA C++ standard
set(CMAKE_CUDA_STANDARD 14) # Should be changed to 17. Needs CMake 3.18
set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)


# Set default install path
include(GNUInstallDirs)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}"
    CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
endif()

# Set default indent style
set(QUICC_CMAKE_INDENT "  " CACHE STRING "Indent for cmake output.")
mark_as_advanced(QUICC_CMAKE_INDENT)

#
# Enable testing
#
enable_testing()

#
# Include function library
#
include(ModelFunctions)
include(ConfigureFunctions)
include(MessageQuICC)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%-------------------------- CONFIGURABLE -------------------------------%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

option(QUICC_USE_SYSTEM_QUICC "Use system installed QuICC" OFF)

# Setup paths
include(structure/external)
include(structure/global)

# Import framework
if(QUICC_USE_SYSTEM_QUICC)
  find_package(quicc REQUIRED)
else()
  # Create global configurator
  include(setup/global)

  # Framework setup
  include(setup/Framework)
endif()

# Scripts setup
include(setup/Scripts)

# Debug setup
include(setup/Debug)

# Automatic setup
include(setup/automatic)

# Components setup
include(setup/Components)

# Models setup
include(setup/Models)

# CUDA setup
include(setup/CUDA)
