## ------------------------------------------------------------------------- ##
 # mk/cmake/fetch_provider_packages/FetchRetroShare.cmake
 # This file is part of RetroShare-GUI.
 #
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

# 3.11: FetchContent
# 3.14: FetchContent_MakeAvailable()
cmake_minimum_required(VERSION 3.14...4.4)

set(FETCH_RetroShare_GIT_TAG "" CACHE STRING
	"commit of libretroshare to fetch instead of the one the repository pins"
)
mark_as_advanced(FETCH_RetroShare_GIT_TAG)

# This program lives in the RetroShare repository, next to the libretroshare
# submodule, and that submodule is the libretroshare to build with: its sources
# when it is checked out, else the commit the repository pins, fetched from the
# submodule's address. Never a branch: the sources built would depend on the
# day of the build.
get_filename_component(rs_repository_dir
	"${CMAKE_CURRENT_LIST_DIR}/../../../.." ABSOLUTE
)
if(EXISTS "${rs_repository_dir}/libretroshare/CMakeLists.txt")
	set(rs_libretroshare_source SOURCE_DIR "${rs_repository_dir}/libretroshare")
else()
	find_program(RS_GIT_EXECUTABLE git)
	mark_as_advanced(RS_GIT_EXECUTABLE)
	set(rs_libretroshare_url "")
	if(RS_GIT_EXECUTABLE)
		if(NOT FETCH_RetroShare_GIT_TAG)
			execute_process(
				COMMAND "${RS_GIT_EXECUTABLE}" rev-parse HEAD:libretroshare
				WORKING_DIRECTORY "${rs_repository_dir}"
				OUTPUT_VARIABLE FETCH_RetroShare_GIT_TAG
				OUTPUT_STRIP_TRAILING_WHITESPACE
				ERROR_QUIET
			)
		endif()
		execute_process(
			COMMAND "${RS_GIT_EXECUTABLE}" config -f .gitmodules
				submodule.libretroshare.url
			WORKING_DIRECTORY "${rs_repository_dir}"
			OUTPUT_VARIABLE rs_libretroshare_url
			OUTPUT_STRIP_TRAILING_WHITESPACE
			ERROR_QUIET
		)
	endif()
	if(NOT FETCH_RetroShare_GIT_TAG OR NOT rs_libretroshare_url)
		message(FATAL_ERROR
			"libretroshare: no sources in ${rs_repository_dir}/libretroshare, "
			"and the commit the repository pins cannot be read (not a git "
			"checkout?). Check out the libretroshare submodule, or set "
			"FETCHCONTENT_SOURCE_DIR_RETROSHARE to its sources."
		)
	endif()
	set(rs_libretroshare_source
		GIT_REPOSITORY "${rs_libretroshare_url}"
		GIT_TAG ${FETCH_RetroShare_GIT_TAG}
		GIT_PROGRESS TRUE
	)
endif()

set(RS_USE_I2P_SAM3 TRUE)

include(FetchContent)
FetchContent_Declare(RetroShare
	${rs_libretroshare_source}
)
FetchContent_MakeAvailable(RetroShare)
set(RetroShare_FOUND TRUE)
