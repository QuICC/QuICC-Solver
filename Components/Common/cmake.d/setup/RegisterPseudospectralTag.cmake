set(idNS "PseudospectralTag")
string(TOUPPER ${idNS} IDNS)
set(idID "PseudospectralTag")
string(TOUPPER ${idID} IDID)
# glob all ID files
file(GLOB AllFiles include/QuICC/${idNS}/*.hpp)
# list of excluded files
set(excluded 
  "Typedefs.hpp"
  )

# generate header and id list
set(idList "")
set(idHeader "")
foreach(fullname ${AllFiles})
  get_filename_component(fname ${fullname} NAME)
  if(NOT fname IN_LIST excluded)
    set(idHeader "${idHeader}\n#include \"QuICC/${idNS}/${fname}\"")
    get_filename_component(cname ${fullname} NAME_WE)
    set(idList "${idList}\n      ${cname}::id();")
  endif()
endforeach(fullname)

# Configure registration file
configure_file(
  "include/QuICC/IdTools/registerAll.hpp.in"
  "${PROJECT_BINARY_DIR}/include/QuICC/${idNS}/registerAll.hpp"
  )

# Configure ID interface registration files
configure_file(
  "include/QuICC/IdTools/IId.hpp.in"
  "${PROJECT_BINARY_DIR}/include/QuICC/${idNS}/I${idID}.hpp"
  )
configure_file(
  "include/QuICC/IdTools/IRegisterId.hpp.in"
  "${PROJECT_BINARY_DIR}/include/QuICC/${idNS}/IRegisterId.hpp"
  )
configure_file(
  "include/QuICC/IdTools/ICreator.hpp.in"
  "${PROJECT_BINARY_DIR}/include/QuICC/${idNS}/ICreator.hpp"
  )
configure_file(
  "include/QuICC/IdTools/CreatorImpl.hpp.in"
  "${PROJECT_BINARY_DIR}/include/QuICC/${idNS}/CreatorImpl.hpp"
  )
configure_file(
  "include/QuICC/IdTools/Coordinator.hpp.in"
  "${PROJECT_BINARY_DIR}/include/QuICC/${idNS}/Coordinator.hpp"
  )

# unset all variables
unset(idNS)
unset(IDNS)
unset(idID)
unset(IDID)
unset(idARG)
unset(idSIG)
unset(idCOMMA)
unset(idNAN)
unset(idList)
unset(idHeader)
