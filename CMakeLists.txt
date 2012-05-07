project(matlab_io)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

###########################################################
# find and configure the isis core library
###########################################################
FIND_PACKAGE(ISIS REQUIRED)
INCLUDE_DIRECTORIES(${ISIS_INCLUDE_DIRS})

###############################################################################
# add global options
###############################################################################
ADD_DEFINITIONS(-D_ENABLE_LOG=1)
ADD_DEFINITIONS(-D_ENABLE_DEBUG=1)

###############################################################################
# add a target
###############################################################################
find_library(LIB_MWMCLMCR "mwmclmcr")
find_path(INCPATH_MAT mat.h PATH_SUFFIXES "include")

include_directories(${INCPATH_MAT})
add_library(isisImageFormat_matlab SHARED imageFormat_matlab.cpp)
target_link_libraries(isisImageFormat_matlab isis_core ${LIB_MWMCLMCR} ${ISIS_LIB_DEPENDS})
set_target_properties(isisImageFormat_matlab PROPERTIES INSTALL_RPATH_USE_LINK_PATH yes) #needed if -lmwmclmcr is not in the searchpath

if(WIN32)
	install(TARGETS isisImageFormat_matlab RUNTIME DESTINATION bin COMPONENT "IO_plugins")
else(WIN32)
	install(TARGETS isisImageFormat_matlab DESTINATION ${ISIS_PLUGIN_INFIX} COMPONENT "IO_plugins" )
endif(WIN32)
