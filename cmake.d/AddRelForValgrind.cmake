set(_languages "CXX;C")

foreach(_lang IN LISTS _languages)
  message(DEBUG "_lang: ${_lang}")
  set(CMAKE_${_lang}_FLAGS_RELFORVALGRIND_INIT "-Wall -O2 -DBOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS -g")
endforeach()
