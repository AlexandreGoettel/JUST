# General config
cmake_minimum_required(VERSION 3.4.1 FATAL_ERROR)
project(NuSolarFit
	VERSION 0.0.1
	DESCRIPTION "Preliminary tests for juno solar fits.")

SET(EXENAME ${PROJECT_NAME})

#-------------------------------------------------------------------------------
# ADD ROOT
#-------------------------------------------------------------------------------
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

# Make sure the ROOT variables are set correctly
message(STATUS "ROOT INCLUDE DIR: " ${ROOT_INCLUDE_DIRS})

find_package(ROOT 6.22 CONFIG REQUIRED Minuit)

# This allows to include ROOT headers in the project files
INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIRS})

#-------------------------------------------------------------------------------
# Set Variables
#-------------------------------------------------------------------------------
SET(EXENAME ${PROJECT_NAME})
SET(INCLUDEPATH ${CMAKE_CURRENT_SOURCE_DIR}/include)
SET(SRCPATH ${CMAKE_CURRENT_SOURCE_DIR}/src)
SET(SOURCE
	${SRCPATH}/DataReader.cpp
	${SRCPATH}/FitResults.cpp
	${SRCPATH}/Fitter.cpp
	${SRCPATH}/Main.cpp
	${SRCPATH}/OutputManager.cpp
	${SRCPATH}/Parser.cpp
	${SRCPATH}/ToyDataGenerator.cpp
)
SET(HEADERS
	${INCLUDEPATH}/DataReader.h
	${INCLUDEPATH}/FitResults.h
	${INCLUDEPATH}/Fitter.h
	${INCLUDEPATH}/OutputManager.h
	${INCLUDEPATH}/Parser.h
	${INCLUDEPATH}/ToyDataGenerator.h
)

#-------------------------------------------------------------------------------
# Add Executables here
#-------------------------------------------------------------------------------
add_executable(${EXENAME} ${SOURCE} ${HEADERS})

### Link libraries here ###
target_link_libraries(${EXENAME} ${ROOT_LIBRARIES})
### Make sure the compiler can find the include files
target_include_directories(${EXENAME} PRIVATE ${INCLUDEPATH})
### Install ###
INSTALL(TARGETS ${EXENAME} DESTINATION bin)
