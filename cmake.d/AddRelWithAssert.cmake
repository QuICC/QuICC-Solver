set(_languages "CXX;C")

foreach(_lang IN LISTS _languages)
  message(DEBUG "_lang: ${_lang}")
  set(CMAKE_${_lang}_FLAGS_RELWITHASSERT_INIT "-Wall -O2")
endforeach()
