# Install script for directory: /fasttmp/zhanhy41/Code/scorec_phasta/phasta

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/phastaIO/cmake_install.cmake")
  include("/fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/shapeFunction/cmake_install.cmake")
  include("/fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/phSolver/cmake_install.cmake")
  include("/fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/converterIO/cmake_install.cmake")
  include("/fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/AcuStat/cmake_install.cmake")
  include("/fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/M2N/cmake_install.cmake")
  include("/fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/M2NFixBnd/cmake_install.cmake")
  include("/fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/checkphasta/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

file(WRITE "/fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/${CMAKE_INSTALL_MANIFEST}" "")
foreach(file ${CMAKE_INSTALL_MANIFEST_FILES})
  file(APPEND "/fasttmp/zhanhy41/Code/scorec_phasta/phasta/build_ehsan/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
endforeach()
