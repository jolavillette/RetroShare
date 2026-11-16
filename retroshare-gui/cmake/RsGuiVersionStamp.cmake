################################################################################
# retroshare-gui/cmake/RsGuiVersionStamp.cmake                                 #
# Copyright (C) 2026, Retroshare team <contact@retroshare.cc>                  #
#                                                                              #
# This program is free software: you can redistribute it and/or modify         #
# it under the terms of the GNU Affero General Public License as               #
# published by the Free Software Foundation, either version 3 of the           #
# License, or (at your option) any later version.                              #
#                                                                              #
# This program is distributed in the hope that it will be useful,              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                #
# GNU Affero General Public License for more details.                          #
#                                                                              #
# You should have received a copy of the GNU Affero General Public License     #
# along with this program.  If not, see <https://www.gnu.org/licenses/>.       #
################################################################################

# Build time version stamping, run with `cmake -P` at every build.
#
# Writes rsguiversion_build.h, which holds RS_GUI_VERSION, the `git describe`
# of the RetroShare super-project. Only rshare.cpp (through rsguiversion.h) and
# the Windows resource file include that header, so a new commit rebuilds them
# instead of the whole GUI, and the version always describes the commit that is
# being built.
#
# The file is written through configure_file(), which leaves it untouched when
# its content does not change: a build without a new commit rebuilds nothing.
#
# Expected variables:
#   RS_STAMP_SOURCE_DIR  retroshare-gui source directory
#   RS_STAMP_TEMPLATE    rsguiversion_build.h.in
#   RS_STAMP_OUTPUT      rsguiversion_build.h to write
#   GIT_EXECUTABLE       git, may be empty

set(RS_GUI_VERSION_SET OFF)
set(RS_GUI_VERSION "")
if(GIT_EXECUTABLE)
	execute_process(
		COMMAND "${GIT_EXECUTABLE}" describe --tags --always --dirty
		WORKING_DIRECTORY "${RS_STAMP_SOURCE_DIR}"
		OUTPUT_VARIABLE RS_GUI_VERSION
		OUTPUT_STRIP_TRAILING_WHITESPACE
		ERROR_QUIET
		RESULT_VARIABLE V_GIT_RESULT )
	if(NOT V_GIT_RESULT EQUAL 0)
		set(RS_GUI_VERSION "")
	endif()
endif()

if(NOT RS_GUI_VERSION STREQUAL "")
	# Drop the leading 'v' to match the usual version style
	string(REGEX REPLACE "^v" "" RS_GUI_VERSION "${RS_GUI_VERSION}")
	set(RS_GUI_VERSION_SET ON)
endif()

configure_file("${RS_STAMP_TEMPLATE}" "${RS_STAMP_OUTPUT}" @ONLY)
