set(QUICC_WIGXJPF_LIB "quicc_3rd_wigxjpf")
set(QUICC_WIGXJPF_PREFIX "${PROJECT_BINARY_DIR}/3rd_party")
set(QUICC_WIGXJPF_SOURCE_DIR "${QUICC_WIGXJPF_PREFIX}/wigxjpf")

include(ExternalProject)
find_program(MAKE_EXE NAMES gmake nmake make)
ExternalProject_Add(extern_wigxjpf
  URL               http://fy.chalmers.se/subatom/wigxjpf/wigxjpf-1.13.tar.gz
  URL_HASH          SHA256=90ab9bfd495978ad1fdcbb436e274d6f4586184ae290b99920e5c978d64b3e6a
  SOURCE_DIR        ${QUICC_WIGXJPF_SOURCE_DIR}
  BUILD_IN_SOURCE ON
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
  UPDATE_COMMAND ""
  BUILD_BYPRODUCTS ${QUICC_WIGXJPF_SOURCE_DIR}/lib/libwigxjpf.a
  DOWNLOAD_EXTRACT_TIMESTAMP OFF
)

add_library(${QUICC_WIGXJPF_LIB} STATIC IMPORTED GLOBAL)
add_dependencies(${QUICC_WIGXJPF_LIB} extern_wigxjpf)
set(QUICC_WIGXJPF_INCLUDE_DIRS ${QUICC_WIGXJPF_SOURCE_DIR}/inc)
# This is a workaround for the fact that included directories of an imported
# target should exist in the filesystem already at the configuration time.
# ref: https://gitlab.kitware.com/cmake/cmake/-/issues/15052
file(MAKE_DIRECTORY ${QUICC_WIGXJPF_INCLUDE_DIRS})
set_target_properties(
  ${QUICC_WIGXJPF_LIB}
  PROPERTIES IMPORTED_LOCATION ${QUICC_WIGXJPF_SOURCE_DIR}/lib/libwigxjpf.a
             INTERFACE_INCLUDE_DIRECTORIES ${QUICC_WIGXJPF_INCLUDE_DIRS}
             INTERFACE_SYSTEM_INCLUDE_DIRECTORIES ${QUICC_WIGXJPF_PREFIX})
#target_link_libraries(${QUICC_WIGXJPF_LIB} INTERFACE ???)

# Alias
add_library(${QUICC_NAMESPACE}ThirdParty::Wigxjpf ALIAS
  "${QUICC_WIGXJPF_LIB}")
