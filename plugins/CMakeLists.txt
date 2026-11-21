## ------------------------------------------------------------------------- ##
 # plugins/CMakeLists.txt
 # This file is part of RetroShare.
 #
 # Copyright (C) 2026      David Bears <dbear4q@gmail.com>
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

cmake_minimum_required(VERSION 3.0...4.4)

# Per-plugin toggles: allow building a subset of plugins. Useful when one plugin
# is mid-port (e.g. to Qt6) and would otherwise break the whole build: disable
# it with -DRS_PLUGIN_<NAME>=OFF and keep compiling the rest.
option(RS_PLUGIN_FEEDREADER "Build the FeedReader plugin" ON)
option(RS_PLUGIN_VOIP       "Build the VOIP plugin"       ON)

if(RS_GUI)
	# don't accidentally find an installed version of RetroShare-GUI
	include(RedirectFindPackage)
	redirect_find_package(RetroShare-GUI)

	# A plugin is loaded into the running program and uses the same
	# libretroshare as it. That only holds when libretroshare is a shared
	# library: linked statically, every plugin would carry its own copy of the
	# whole library, with its own global state.
	find_package(RetroShare REQUIRED)
	get_target_property(retroshare_type retroshare::retroshare TYPE)
	if(NOT retroshare_type STREQUAL "SHARED_LIBRARY")
		message(FATAL_ERROR
			"RetroShare plugins need a shared libretroshare, but it is a "
			"${retroshare_type}. Configure with -DRS_LIBRETROSHARE_SHARED=ON, "
			"or without plugins (-DRS_PLUGINS=OFF)."
		)
	endif()

	# Install the plugins where the engine looks for them: RS_PLUGIN_DIR_FROM_EXE
	# (set by the super-project) relative to the installed program.
	include(GNUInstallDirs)
	if(NOT DEFINED RS_PLUGIN_DIR_FROM_EXE)
		set(RS_PLUGIN_DIR_FROM_EXE "../lib/retroshare/extensions6")
	endif()
	if(APPLE)
		set(rs_gui_exe_install_dir
			"${CMAKE_INSTALL_BINDIR}/retroshare-gui.app/Contents/MacOS"
		)
	else()
		set(rs_gui_exe_install_dir "${CMAKE_INSTALL_BINDIR}")
	endif()
	get_filename_component(RS_PLUGIN_INSTALL_DIR
		"/${rs_gui_exe_install_dir}/${RS_PLUGIN_DIR_FROM_EXE}" ABSOLUTE
	)
	string(REGEX REPLACE "^/" "" RS_PLUGIN_INSTALL_DIR "${RS_PLUGIN_INSTALL_DIR}")

	if(RS_PLUGIN_FEEDREADER)
		add_subdirectory(FeedReader)
	endif()
	if(RS_PLUGIN_VOIP)
		add_subdirectory(VOIP)
	endif()
endif()
