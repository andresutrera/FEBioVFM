cmake_minimum_required(VERSION 3.1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(FEBIOPlugins LANGUAGES CXX)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_FILES_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if(NOT name)
    set(name febioplugin)
endif()


# Set a default build type if none was specified
set(default_build_type "Release")
 
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Find FEBio SDK or git repo automatically
# FEBio
# Find FEBio SDK or git repo automatically
if(WIN32)
	set(TEMP_PATHS ${CMAKE_SOURCE_DIR}/.. ${CMAKE_SOURCE_DIR}/../.. $ENV{HOMEPATH}/ $ENV{HOMEPATH}/source/repos $ENV{HOMEPATH}/*)
else()
    set(TEMP_PATHS $ENV{HOME}/ $ENV{HOME}/*)
endif()
    
find_path(FEBio_SDK FECore/stdafx.h
    PATHS ${TEMP_PATHS}
    PATH_SUFFIXES FEBio
    DOC "Path to the FEBio SDK, or git repo.")
    
if(NOT FEBio_SDK)
    if(WIN32)
        set(TEMP_PATHS $ENV{PROGRAMFILES}/* ${CMAKE_SOURCE_DIR}/.. $ENV{HOMEPATH}/* )
    elseif(APPLE)
        set(TEMP_PATHS /Applications/* ${CMAKE_SOURCE_DIR}/.. $ENV{HOME}/*)
    else()
        set(TEMP_PATHS ${CMAKE_SOURCE_DIR}/.. $ENV{HOME}/*)
    endif() 
    
    find_path(FEBio_SDK "include/FECore/stdafx.h"
        PATHS ${TEMP_PATHS}
        PATH_SUFFIXES sdk
        DOC "Path to the FEBio SDK, or git repo.")
endif()

if(NOT FEBio_SDK)
    set(FEBio_SDK "" CACHE PATH "Path to the FEBio SDK, or git repo.")
    message(FATAL_ERROR "Unable to find path to FEBio SDK or git repo automatically. Please set FEBio_SDK to the path to your FEBio SDK or git repo.")
endif()

# Only update the include and lib directories if the FEBio_SDK path has been changed. 
if(NOT OLD_SDK)
    set(NEWPATH TRUE)
else()
    #cmake_path(CONVERT ${OLD_SDK} TO_CMAKE_PATH_LIST STD_OLD_SDK)
    #cmake_path(CONVERT ${FEBio_SDK} TO_CMAKE_PATH_LIST STD_FEBIO_SDK)
    #string(COMPARE NOTEQUAL ${STD_FEBIO_SDK} ${STD_OLD_SDK} NEWPATH)
    string(COMPARE NOTEQUAL ${FEBio_SDK} ${OLD_SDK} NEWPATH)
endif()

if(NEWPATH)
    # Is this the SDK?
    string(REGEX MATCH "sdk" IS_SDK ${FEBio_SDK})

    set(LIB_SUFFIXES "")
    if(IS_SDK)
        set(FEBio_INC "${FEBio_SDK}/include" CACHE PATH "Path to FEBio include directory." FORCE)
        
        if(WIN32)
            list(APPEND LIB_SUFFIXES "vs2017/Release" "vs2017/Debug")
        else()
            list(APPEND LIB_SUFFIXES "lib")
        endif()
    else()
        set(FEBio_INC ${FEBio_SDK} CACHE PATH "Path to FEBio include directory." FORCE)
        
        if(WIN32)
            list(APPEND LIB_SUFFIXES "cmbuild/lib/Release" "cmbuild/lib/Debug" "cbuild/lib/Release" "cbuild/lib/Debug" "build/lib/Release" "build/lib/Debug")
        else()
            list(APPEND LIB_SUFFIXES "cbuild/lib" "cmbuild/lib" "build/lib" "cbuild/Release/lib" "cmbuild/Release/lib" "build/Release/lib" "cbuild/Debug/lib" "cmbuild/Debug/lib" "build/Debug/lib")
        endif()
    endif()

    mark_as_advanced(FEBio_INC)

    # Find lib path
    find_library(FECORE  
        NAMES FECore fecore fecore_gcc64 fecore_lnx64
        PATHS ${FEBio_SDK}
        PATH_SUFFIXES ${LIB_SUFFIXES}
        DOC "FEBio library path")

    if(FECORE)
        get_filename_component(FECORE_TEMP ${FECORE} DIRECTORY)
        set(FEBio_LIB_DIR ${FECORE_TEMP} CACHE PATH "Path to the FEBio lib directory." FORCE)
        mark_as_advanced(FEBio_LIB_DIR)
        unset(FECORE_TEMP)
        unset(FECORE CACHE)
    else()
        set(FEBio_LIB_DIR CACHE PATH "Path to the FEBio lib directory." FORCE)
        message(SEND_ERROR "Unable to find FEBio Library path automatically. Set FEBio_LIB_DIR.")
        unset(FECORE CACHE)
    endif()
endif()

set(OLD_SDK ${FEBio_SDK} CACHE INTERNAL "Old SDK path.")
mark_as_advanced(OLD_SDK)

#### Check if paths are valid and find libs ####
function(findLib libDir libOut libName)
    # Optional arguments ARGV3-ARVG6 can be used as alernative names for the library
    find_library(TEMP NAMES ${libName} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6}
        PATHS ${${libDir}} NO_DEFAULT_PATH)
    
    if(TEMP)
        set(${libOut} ${TEMP} PARENT_SCOPE)
        unset(TEMP CACHE)
    else()
        if(WIN32)
            message(SEND_ERROR "Could not find ${libName}.lib. Check ${libDir}.")
        elseif(APPLE)
            message(SEND_ERROR "Could not find lib${libName}.so, lib${libName}.a, or lib${libName}.dylib Check ${libDir}.")
        else()
            message(SEND_ERROR "Could not find lib${libName}.so, or lib${libName}.a. Check ${libDir}")
        endif()
        unset(TEMP CACHE)
    endif()
endfunction()


function(findLib2 libName libDir libOut)
    # Optional arguments ARGV3-ARVG6 can be used as alernative names for the library
    if(WIN32)
        find_library(TEMP NAMES ${libName}.lib ${ARGV3}.lib ${ARGV4}.lib ${ARGV5}.lib ${ARGV6}.lib
            PATHS ${${libDir}} NO_DEFAULT_PATH)
    else()
        find_library(TEMP NAMES lib${libName}.a lib${ARGV3}.a lib${ARGV4}.a lib${ARGV5}.a lib${ARGV6}.a
            PATHS ${${libDir}} NO_DEFAULT_PATH)
    endif()
    
    if(TEMP)
        set(${libOut} ${TEMP} PARENT_SCOPE)
        unset(TEMP CACHE)
    else()
        if(WIN32)
            message(SEND_ERROR "Could not find ${libName}.lib. Check ${libName}.")
        elseif(APPLE)
            message(SEND_ERROR "Could not find lib${libName}.so, lib${libName}.a, or lib${libName}.dylib Check ${libDir}.")
        else()
            message(SEND_ERROR "Could not find lib${libName}.so, or lib${libName}.a. Check ${libDir}")
        endif()
        unset(TEMP CACHE)
    endif()
endfunction()

set(MKL_LIB_DIR "/opt/intel/oneapi/mkl/latest/lib")
set(MKL_OMP_LIB "/opt/intel/oneapi/compiler/latest/lib/libiomp5.so")
if(1)
    set(MKL_LIB_NAMES mkl_intel_lp64 mkl_core mkl_intel_thread)
    foreach(namelib IN LISTS MKL_LIB_NAMES)
        findLib2(${namelib} MKL_LIB_DIR MKL_LIB)
        list(APPEND MKL_LIBS ${MKL_LIB})
        unset(MKL_LIB)
    endforeach(namelib)
    
    if(NOT EXISTS ${MKL_OMP_LIB})
        message(SEND_ERROR "Could not find MKL OMP library. Check MKL_OMP_LIB.")
    endif()
endif()


if (WIN32)
    set(FEBio_DEBUG_LIBS "")
    set(FEBio_RELEASE_LIBS "")

    set(FEBio_LIB_NAMES fecore febiolib febiomech febiomix xml febioplot febiorve)
    set(FEBio_DEBUG_DIR ${FEBio_LIB_DIR}/Debug)
    set(FEBio_RELEASE_DIR ${FEBio_LIB_DIR}/Release)

    foreach(name IN LISTS FEBio_LIB_NAMES)
        findLib(FEBio_DEBUG_DIR FEBio_LIB ${name})
        list(APPEND FEBio_DEBUG_LIBS ${FEBio_LIB})
        unset(FEBio_LIB)
    endforeach(name)

    foreach(name IN LISTS FEBio_LIB_NAMES)
        findLib(FEBio_RELEASE_DIR FEBio_LIB ${name})
        list(APPEND FEBio_RELEASE_LIBS ${FEBio_LIB})
        unset(FEBio_LIB)
    endforeach(name)

    unset(FEBio_LIB_NAMES)
else()
    set(FEBio_LIBS "")

    set(FEBio_LIB_NAMES fecore febiolib febiomech )
    set(FEBio_LIB_DIR "${FEBio_SDK}/lib")
    foreach(name IN LISTS FEBio_LIB_NAMES)
        findLib(FEBio_LIB_DIR FEBio_LIB ${name})
        list(APPEND FEBio_LIBS ${FEBio_LIB})
        unset(FEBio_LIB)
    endforeach(name)

    unset(FEBio_LIB_NAMES)
endif()






if(APPLE)
    # OpenMP
    find_package(OpenMP QUIET)
endif()

##### Set appropriate defines and includes #####

if(WIN32)
    add_definitions(-DWIN32 /MP /openmp /D_ALLOW_RANDOM_DISTRIBUTION_CONST_OPERATOR)
elseif(APPLE)
    add_definitions(-D__APPLE__)
    set(CMAKE_OSX_DEPLOYMENT_TARGET 10.13)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
    
    # Add openMP flags for macOS if found.
    if(${OpenMP_C_FOUND})
        add_compile_options(-Xpreprocessor -fopenmp)
        include_directories(${OpenMP_C_INCLUDE_DIRS})
    endif()
else()
	add_definitions(-DLINUX)
    add_compile_options(-static-libstdc++ -static-libgcc -w -Wall -fopenmp)
    
    set(CMAKE_BUILD_RPATH_USE_LINK_PATH FALSE)
    set(CMAKE_BUILD_RPATH $ORIGIN/../lib/)
endif()

# Extra compiler flags for intel compiler
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	set(CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS} -static-intel)
endif()

##### Add library #####

file(GLOB_RECURSE HDR CONFIGURE_DEPENDS "src/*.h" "src/*.hpp")
file(GLOB_RECURSE SRC CONFIGURE_DEPENDS "src/*.cpp")

add_library(${name} SHARED ${HDR} ${SRC})
set_property(TARGET ${name} PROPERTY AUTOGEN_BUILD_DIR ${CMAKE_BINARY_DIR}/CMakeFiles/AutoGen/${name}_autogen)
target_include_directories(${name} PRIVATE ${FEBio_INC} ${CMAKE_SOURCE_DIR}/src)

##### Link Libraries #####

if (WIN32)

    foreach(name IN LISTS FEBio_RELEASE_LIBS)
        target_link_libraries(${name} optimized ${name})
    endforeach()

    foreach(name IN LISTS FEBio_DEBUG_LIBS)
        target_link_libraries(${name} debug ${name})
    endforeach()
else()
    target_link_libraries(${name} PRIVATE -Wl,--start-group ${MKL_OMP_LIB} ${MKL_LIBS} -Wl,--end-group ${FEBio_LIBS} )
endif()
