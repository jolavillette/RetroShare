## ------------------------------------------------------------------------- ##
 # CMakeLists.txt
 # This file is part of RetroShare-SuperProject.
 #
 # Copyright (C) 2021-2022 Gioacchino Mazzurco <gio@retroshare.cc>
 # Copyright (C) 2026      David Bears <dbear4q@gmail.com>
 #
 # SPDX-License-Identifier: AGPL-3.0-or-later
 #
 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU Affero General Public License as
 # published by the Free Software Foundation, either version 3 of the
 # License, or (at your option) any later version.
 #
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU Affero General Public License for more details.
 #
 # You should have received a copy of the GNU Affero General Public License
 # along with this program.  If not, see <https://www.gnu.org/licenses/>.
## ------------------------------------------------------------------------- ##

# 3.7 : if(VERSION_GREATER_EQUAL)
cmake_minimum_required(VERSION 3.7...4.4)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
	list(PREPEND CMAKE_PROJECT_TOP_LEVEL_INCLUDES
		"${CMAKE_CURRENT_SOURCE_DIR}/mk/cmake/fetch_provider.cmake"
	)
else()
	include(mk/cmake/fetch_provider.cmake)
endif()

project(retroshare-superproject)

# Configuration Options
option(RS_GUI "Build RetroShare GUI" ON)
option(RS_SERVICE "Build RetroShare Service" ON)
option(RS_FRIENDSERVER "Build RetroShare Friendserver" ON)
option(RS_PLUGINS "Build RetroShare Plugins" OFF)
# Windows has no run-time search path: a program finds its DLLs in its own
# directory. Put the programs and the DLLs built here (a shared libretroshare)
# in one directory, so that what is built runs from the build tree.
if(WIN32 AND NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
endif()

# The engine looks for the plugins shipped with the program next to the running
# binary, in ../lib/retroshare/extensions6 (libretroshare, rsinit.cc). Build
# them there, so that the program run from the build tree finds them; they are
# installed to the same place relative to the installed binary.
set(RS_PLUGIN_DIR_FROM_EXE "../lib/retroshare/extensions6")
if(APPLE)
	set(rs_gui_exe_build_dir
		${PROJECT_BINARY_DIR}/retroshare-gui/retroshare-gui.app/Contents/MacOS
	)
elseif(CMAKE_RUNTIME_OUTPUT_DIRECTORY)
	set(rs_gui_exe_build_dir ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
else()
	set(rs_gui_exe_build_dir ${PROJECT_BINARY_DIR}/retroshare-gui)
endif()
get_filename_component(rs_plugin_binary_dir_default
	"${rs_gui_exe_build_dir}/${RS_PLUGIN_DIR_FROM_EXE}" ABSOLUTE
)
set(RS_PLUGIN_BINARY_DIR ${rs_plugin_binary_dir_default} CACHE PATH
	"directory to place plugin shared objects at build time (does not affect install location)"
)


# Default to Release build type if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	set(CMAKE_BUILD_TYPE "Release" CACHE
		STRING "Choose the type of build (Debug or Release)" FORCE
	)
endif()


### some misc global properties

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# export symbols for plugin/module linking
set(CMAKE_ENABLE_EXPORTS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(APPLE)
	# On macOS, allow undefined symbols in shared modules/libraries at link time
	# because they will be resolved at runtime by the host executable (retroshare-gui)
	set(CMAKE_SHARED_LINKER_FLAGS
		"${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup"
	)
	set(CMAKE_MODULE_LINKER_FLAGS
		"${CMAKE_MODULE_LINKER_FLAGS} -undefined dynamic_lookup"
	)
endif()


################################################################################
### dependencies

list(APPEND CMAKE_MODULE_PATH
	"${PROJECT_SOURCE_DIR}/mk/cmake"
	"${PROJECT_SOURCE_DIR}/mk/cmake/fetch_provider_packages"
)
set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)

# The dependencies that libretroshare can fetch are also submodules of this
# repository, at the same commits. When a submodule is checked out, build from
# it: nothing is downloaded. FETCHCONTENT_SOURCE_DIR_<NAME> is FetchContent's
# own switch for local sources; FETCH_sam3_SOURCE_DIR is its equivalent for
# libsam3, which is not built by CMake.
set(rs_dependency_submodules
	"RNP|supportlibs/librnp|src/libsexpp/CMakeLists.txt"
	"RESTBED|supportlibs/restbed|CMakeLists.txt"
	"UDP-DISCOVERY-CPP|supportlibs/udp-discovery-cpp|CMakeLists.txt"
	"RAPIDJSON|supportlibs/rapidjson|include/rapidjson/rapidjson.h"
	"BITDHT|supportlibs/libbitdht|CMakeLists.txt"
	"OPENPGPSDK|supportlibs/openpgpsdk|CMakeLists.txt"
	"JNI.HPP|supportlibs/jni.hpp|CMakeLists.txt"
	"CMARK|supportlibs/cmark|CMakeLists.txt"
)
foreach(rs_dependency_submodule ${rs_dependency_submodules})
	string(REPLACE "|" ";" rs_dependency_submodule "${rs_dependency_submodule}")
	list(GET rs_dependency_submodule 0 rs_dependency_name)
	list(GET rs_dependency_submodule 1 rs_dependency_path)
	list(GET rs_dependency_submodule 2 rs_dependency_marker)
	set(rs_dependency_dir "${PROJECT_SOURCE_DIR}/${rs_dependency_path}")
	if(NOT FETCHCONTENT_SOURCE_DIR_${rs_dependency_name}
		AND EXISTS "${rs_dependency_dir}/${rs_dependency_marker}")
		set(FETCHCONTENT_SOURCE_DIR_${rs_dependency_name} "${rs_dependency_dir}"
			CACHE PATH "local sources of ${rs_dependency_name} (submodule)"
		)
	endif()
endforeach()
if(NOT FETCH_sam3_SOURCE_DIR
	AND EXISTS "${PROJECT_SOURCE_DIR}/supportlibs/libsam3/Makefile")
	set(FETCH_sam3_SOURCE_DIR "${PROJECT_SOURCE_DIR}/supportlibs/libsam3"
		CACHE PATH "local sources of libsam3 (submodule)"
	)
endif()

# avoid using system libretroshare unless the user specifically asks for it
set(FETCH_PROVIDER_RetroShare_MODE FETCH_ONLY CACHE STRING
	"control how the Fetch dependency provider handles RetroShare"
)

if(RS_GXSWIKIPOS)
	add_subdirectory(supportlibs/pegmarkdown)
endif()

################################################################################
### add subdirectories

if(RS_GUI)
	add_subdirectory(retroshare-gui)

	# I guess it's not possible to use plugins without the GUI?
	if(RS_PLUGINS)
		add_subdirectory(plugins)
	endif()
endif()

if(RS_SERVICE)
	add_subdirectory(retroshare-service)
endif()

if(RS_FRIENDSERVER)
	add_subdirectory(retroshare-friendserver)
endif()

if(RS_WEBUI)
	add_subdirectory(retroshare-webui)
endif()
