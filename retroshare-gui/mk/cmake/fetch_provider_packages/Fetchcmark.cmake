## ------------------------------------------------------------------------- ##
 # mk/cmake/fetch_provider_packages/Fetchcmark.cmake
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

set(FETCH_cmark_GIT_TAG "" CACHE STRING
	"git tag to use for fetching cmark"
)
mark_as_advanced(FETCH_cmark_GIT_TAG)
if(NOT FETCH_cmark_GIT_TAG)
	# A commit, not a branch: everybody builds the same sources, and git does
	# not go back to the network once they are there.
	# RetroShare's cmark submodule.
	set(FETCH_cmark_GIT_TAG
		b9c7a496ba7dd9c3495bae2ff2855899e47b245d
	)
endif()

# of cmark only the library is wanted: no tests, no shared library
set(CMARK_TESTS OFF)
set(CMARK_SHARED OFF)
set(CMARK_STATIC ON)

include(FetchContent)
FetchContent_Declare(cmark
	GIT_REPOSITORY "https://github.com/commonmark/cmark.git"
	GIT_TAG ${FETCH_cmark_GIT_TAG}
	GIT_PROGRESS TRUE
	# keeps cmark's program and its install rules out of our build and packages
	EXCLUDE_FROM_ALL
)
# The pinned cmark asks for a CMake version older than 3.5, which CMake 4
# refuses. CMAKE_POLICY_VERSION_MINIMUM is CMake's way to build such a project
# without editing it; set it for this dependency only.
set(fetch_cmark_policy_version_minimum ${CMAKE_POLICY_VERSION_MINIMUM})
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
FetchContent_MakeAvailable(cmark)
set(CMAKE_POLICY_VERSION_MINIMUM ${fetch_cmark_policy_version_minimum})
unset(fetch_cmark_policy_version_minimum)
set(cmark_FOUND TRUE)

if(NOT TARGET cmark::cmark)
	if(TARGET libcmark_static)
		# The pinned cmark: its static library is the target libcmark_static
		# ("cmark" is the program there), and the target does not carry what a
		# user of the library needs.
		target_include_directories(libcmark_static INTERFACE
			"$<BUILD_INTERFACE:${cmark_SOURCE_DIR}/src>"
			"$<BUILD_INTERFACE:${cmark_BINARY_DIR}/src>"
		)
		target_compile_definitions(libcmark_static INTERFACE CMARK_STATIC_DEFINE)
		add_library(cmark::cmark ALIAS libcmark_static)
	else()
		# later versions of cmark name the library target cmark
		add_library(cmark::cmark ALIAS cmark)
	endif()
endif()
