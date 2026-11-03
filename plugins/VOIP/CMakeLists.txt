## ------------------------------------------------------------------------- ##
 # CMakeLists.txt
 # This file is part of the VOIP plugin for RetroShare.
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

cmake_minimum_required(VERSION 3.0...4.4)

project(VOIPPlugin)

# --- Dependencies ---
# Version-agnostic Qt detection: probe Qt6 first, fall back to Qt5 (mirrors
# retroshare-gui/CMakeLists.txt). All Qt usage goes through ${QT_VERSION_MAJOR}
# so the plugin builds against whichever Qt the GUI uses.
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.22)
	find_package(Qt6 QUIET COMPONENTS Core)
endif()
if(Qt6_FOUND)
	set(QT_VERSION_MAJOR 6)
else()
	find_package(Qt5 COMPONENTS Core REQUIRED)
	set(QT_VERSION_MAJOR 5)
endif()
find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Core Widgets Multimedia Xml Network REQUIRED)
if(QT_VERSION_MAJOR EQUAL 6)
	# Core5Compat provides QRegExp/QTextCodec etc. referenced via shared
	# retroshare-gui headers; matches the GUI's own link.
	find_package(Qt6 COMPONENTS Core5Compat REQUIRED)
endif()
find_package(PkgConfig REQUIRED)
pkg_check_modules(AVCODEC REQUIRED libavcodec)
pkg_check_modules(AVUTIL REQUIRED libavutil)
pkg_check_modules(SPEEX REQUIRED speex)
pkg_check_modules(SPEEXDSP REQUIRED speexdsp)
find_package(RetroShare REQUIRED)
find_package(RetroShare-GUI REQUIRED)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# --- Sources ---
set(VOIP_SOURCES
	VOIPPlugin.cpp
	gui/VOIPConfigPanel.cpp
	services/p3VOIP.cc
	services/rsVOIPItems.cc
	gui/AudioStats.cpp
	gui/AudioWizard.cpp
	gui/SpeexProcessor.cpp
	gui/audiodevicehelper.cpp
	gui/VideoProcessor.cpp
	gui/QVideoDevice.cpp
	gui/VOIPChatWidgetHolder.cpp
	gui/VOIPGUIHandler.cpp
	gui/VOIPNotify.cpp
	gui/VOIPToasterItem.cpp
	gui/VOIPToasterNotify.cpp
)

# macOS + Qt6: the camera is driven straight from AVFoundation (Objective-C++),
# both for the access request and for frame capture. Qt6's QPermission API never
# contacts the OS for the camera on macOS, and its QMediaCaptureSession/QVideoSink
# pipeline delivers no frames there. See gui/MacCamera{Permission,Capture}.{h,mm}.
# Qt5 on macOS keeps the plain Qt camera pipeline.
if(APPLE AND QT_VERSION_MAJOR EQUAL 6)
	list(APPEND VOIP_SOURCES gui/MacCameraPermission.mm gui/MacCameraCapture.mm)
	# Both files are written for ARC (__bridge_transfer, strong properties, no
	# release calls). Without -fobjc-arc nothing they allocate is ever freed: every
	# camera start leaked the capture session, its input and output and the frame
	# grabber.
	set_source_files_properties(gui/MacCameraPermission.mm gui/MacCameraCapture.mm
		PROPERTIES COMPILE_OPTIONS "-fobjc-arc")
endif()

set(VOIP_FORMS
    gui/AudioStats.ui
    gui/AudioWizard.ui
    gui/VOIPConfigPanel.ui
    gui/VOIPToasterItem.ui
)

set(VOIP_RESOURCES
    gui/VOIP_images.qrc
    lang/VOIP_lang.qrc
    qss/VOIP_qss.qrc
)

# --- Translations (.ts -> .qm) ---
file(GLOB TS_FILES "lang/*.ts")
if(Qt${QT_VERSION_MAJOR}LinguistTools_FOUND)
    qt_add_translation(QM_FILES ${TS_FILES})
else()
    find_program(LRELEASE_EXECUTABLE NAMES lrelease-qt${QT_VERSION_MAJOR} lrelease-qt5 lrelease)
    foreach(ts_file ${TS_FILES})
        get_filename_component(ts_name ${ts_file} NAME_WE)
        set(qm_src_file "${CMAKE_CURRENT_SOURCE_DIR}/lang/${ts_name}.qm")
        set(qm_bin_file "${CMAKE_CURRENT_BINARY_DIR}/${ts_name}.qm")

        if(EXISTS ${qm_src_file})
            list(APPEND QM_FILES ${qm_src_file})
        elseif(LRELEASE_EXECUTABLE)
            add_custom_command(
                OUTPUT ${qm_bin_file}
                COMMAND ${LRELEASE_EXECUTABLE} ${ts_file} -qm ${qm_bin_file}
                DEPENDS ${ts_file}
                VERBATIM
            )
            list(APPEND QM_FILES ${qm_bin_file})
        endif()
    endforeach()
endif()

# --- Target ---
add_library(VOIP MODULE
    ${VOIP_SOURCES}
    ${VOIP_FORMS}
    ${VOIP_RESOURCES}
    ${QM_FILES}
)

target_include_directories(VOIP PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/gui
	${CMAKE_CURRENT_SOURCE_DIR}/services
	${CMAKE_CURRENT_SOURCE_DIR}/interface
	${AVCODEC_INCLUDE_DIRS}
	${AVUTIL_INCLUDE_DIRS}
	${SPEEX_INCLUDE_DIRS}
	${SPEEXDSP_INCLUDE_DIRS}
)

target_link_libraries(VOIP PRIVATE
	Qt${QT_VERSION_MAJOR}::Core
	Qt${QT_VERSION_MAJOR}::Widgets
	Qt${QT_VERSION_MAJOR}::Multimedia
	Qt${QT_VERSION_MAJOR}::Xml
	Qt${QT_VERSION_MAJOR}::Network
	${AVCODEC_LIBRARIES}
	${AVUTIL_LIBRARIES}
	${SPEEX_LIBRARIES}
	${SPEEXDSP_LIBRARIES}
	retroshare::retroshare
	retroshare::retroshare-gui
)
if(QT_VERSION_MAJOR EQUAL 6)
	target_link_libraries(VOIP PRIVATE Qt6::Core5Compat)
endif()

# macOS + Qt6: MacCameraPermission.mm / MacCameraCapture.mm call into AVFoundation.
if(APPLE AND QT_VERSION_MAJOR EQUAL 6)
	target_link_libraries(VOIP PRIVATE
		"-framework AVFoundation" "-framework Foundation"
		"-framework CoreVideo" "-framework CoreMedia")
endif()

set_target_properties(VOIP PROPERTIES
    PREFIX "lib"
		AUTOUIC_SEARCH_PATHS "${CMAKE_CURRENT_SOURCE_DIR}/gui"
)

if(RS_PLUGIN_BINARY_DIR)
	set_target_properties(VOIP PROPERTIES
		LIBRARY_OUTPUT_DIRECTORY ${RS_PLUGIN_BINARY_DIR}
		RUNTIME_OUTPUT_DIRECTORY ${RS_PLUGIN_BINARY_DIR}
	)
endif()

# FFmpeg/libav headers (libavcodec/libavutil, included by gui/VideoProcessor.*)
# use the C99 fixed-width constant macros (INT64_C, UINT64_C, ...). In a C++
# translation unit these are only exposed by <stdint.h> when
# __STDC_CONSTANT_MACROS is defined before it is included, so it must be defined
# for the whole target to compile the libav headers.
target_compile_definitions(VOIP PRIVATE __STDC_CONSTANT_MACROS)
if(WIN32)
    target_compile_definitions(VOIP PRIVATE WINDOWS_SYS)
endif()

target_compile_options(VOIP PRIVATE "-Wno-deprecated-declarations")
target_compile_definitions(VOIP PRIVATE RS_NO_WARN_DEPRECATED)

if(RS_DEVELOPMENT_BUILD)
    target_compile_options(VOIP PRIVATE "-g" "-O0")
endif()

# --- Installation ---
install(TARGETS VOIP DESTINATION ${RS_PLUGIN_INSTALL_DIR})
