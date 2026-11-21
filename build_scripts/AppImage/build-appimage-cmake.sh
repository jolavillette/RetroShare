#!/usr/bin/env bash
# build_scripts/AppImage/build-appimage-cmake.sh
#
# Make an AppImage from a CMake build of the super-project.
#
#     build_scripts/AppImage/build-appimage-cmake.sh <build dir> [<output dir>]
#
# The build must have been configured with
#     -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_INSTALL_LIBDIR=lib -DRS_LIBRETROSHARE_SHARED=ON
# (add -DRS_PLUGINS=ON for the plugins). `cmake --install` then lays out exactly
# what an AppImage needs: the programs in usr/bin, libretroshare in usr/lib, the
# plugins in usr/lib/retroshare/extensions6, the data in usr/share/retroshare,
# the desktop file and the icons. linuxdeploy adds Qt and the system libraries.
#
# Tools: linuxdeploy and linuxdeploy-plugin-qt AppImages, looked for in
# $TOOLS (default: build_scripts/AppImage/tools) and downloaded there when
# missing. `file` must be installed (appimagetool needs it). Set
# LDAI_RUNTIME_FILE to an AppImage runtime to work without network.
set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }

[ $# -ge 1 ] || die "usage: $0 <build dir> [<output dir>]"
BUILD="$(cd "$1" && pwd)"
OUT="$(mkdir -p "${2:-$BUILD}" && cd "${2:-$BUILD}" && pwd)"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
TOOLS="${TOOLS:-$HERE/tools}"
APPDIR="$OUT/AppDir"

cache() { sed -n "s/^$1:[A-Z]*=//p" "$BUILD/CMakeCache.txt" | head -1; }
[ -f "$BUILD/CMakeCache.txt" ] || die "$BUILD is not a CMake build directory"
[ "$(cache CMAKE_INSTALL_PREFIX)" = /usr ] || die "configure with -DCMAKE_INSTALL_PREFIX=/usr"
[ "$(cache CMAKE_INSTALL_LIBDIR)" = lib ] || die "configure with -DCMAKE_INSTALL_LIBDIR=lib"
case "$(cache RS_LIBRETROSHARE_SHARED)" in ON|TRUE|1|YES) ;; *) die "configure with -DRS_LIBRETROSHARE_SHARED=ON" ;; esac
command -v file >/dev/null || die "the 'file' command is needed by appimagetool"

LD="$TOOLS/linuxdeploy-x86_64.AppImage"
LD_QT="$TOOLS/linuxdeploy-plugin-qt-x86_64.AppImage"
mkdir -p "$TOOLS"
fetch() { [ -x "$1" ] || { echo ">>> downloading $(basename "$1")"; wget -q -O "$1" "$2" && chmod +x "$1"; }; }
fetch "$LD" https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage
fetch "$LD_QT" https://github.com/linuxdeploy/linuxdeploy-plugin-qt/releases/download/continuous/linuxdeploy-plugin-qt-x86_64.AppImage

# linuxdeploy-plugin-qt finds Qt through qmake: it must be the qmake of the Qt
# the GUI is linked against.
EXE_QT="$(objdump -p "$BUILD/retroshare-gui/retroshare-gui" | grep -oE 'libQt[56]Core' | head -1 | tr -dc 56)"
[ -n "$EXE_QT" ] || die "cannot read the Qt major of $BUILD/retroshare-gui/retroshare-gui"
if [ -z "${QMAKE:-}" ]; then
	case "$EXE_QT" in
		6) QMAKE="$(command -v qmake6 || command -v qmake-qt6 || true)" ;;
		5) QMAKE="$(command -v qmake-qt5 || command -v qmake5 || command -v qmake || true)" ;;
	esac
fi
[ -n "${QMAKE:-}" ] || die "no qmake for Qt$EXE_QT, set QMAKE"
[ "$("$QMAKE" -query QT_VERSION | cut -d. -f1)" = "$EXE_QT" ] || die "$QMAKE is not a Qt$EXE_QT qmake"
export QMAKE APPIMAGE_EXTRACT_AND_RUN=1

echo ">>> install into $APPDIR"
rm -rf "$APPDIR"
DESTDIR="$APPDIR" cmake --install "$BUILD" > "$OUT/appimage-install.log"
# what is only needed to build against RetroShare stays out of the image
rm -rf "$APPDIR/usr/include" "$APPDIR/usr/lib/cmake" "$APPDIR/usr/lib/libretroshare.so" "$APPDIR"/usr/lib/*.a

DESKTOP="$APPDIR/usr/share/applications/retroshare-gui.desktop"
[ -f "$DESKTOP" ] || die "no desktop file installed"
PLUGINS="$APPDIR/usr/lib/retroshare/extensions6"
# the plugins stay where RetroShare looks for them; linuxdeploy only pulls in
# what they link
DEPS_ONLY=()
[ -d "$PLUGINS" ] && DEPS_ONLY=(--deploy-deps-only "$PLUGINS")

echo ">>> linuxdeploy: Qt and system libraries"
"$LD" --appdir "$APPDIR" --executable "$APPDIR/usr/bin/retroshare-gui" \
	--desktop-file "$DESKTOP" "${DEPS_ONLY[@]}" --plugin qt > "$OUT/appimage-deploy.log" 2>&1 \
	|| die "linuxdeploy failed, see $OUT/appimage-deploy.log"

# Qt 6 Multimedia: keep the FFmpeg backend only. linuxdeploy-plugin-qt copies
# the GStreamer backend and libraries but none of the GStreamer plugins, and
# Qt then crashes at the first QMediaPlayer.
QTMM="$APPDIR/usr/plugins/multimedia"
if [ -f "$QTMM/libffmpegmediaplugin.so" ] && [ -f "$QTMM/libgstreamermediaplugin.so" ]; then
	rm -f "$QTMM/libgstreamermediaplugin.so" "$APPDIR"/usr/lib/libgst*-1.0.so* "$APPDIR"/usr/lib/liborc-0.4.so*
fi

echo ">>> linuxdeploy: AppImage"
export LINUXDEPLOY_OUTPUT_VERSION="${LINUXDEPLOY_OUTPUT_VERSION:-$(git -C "$ROOT" describe --tags --always 2>/dev/null || echo dev)}"
( cd "$OUT" && "$LD" --appdir "$APPDIR" "${DEPS_ONLY[@]}" --output appimage > "$OUT/appimage-output.log" 2>&1 ) \
	|| die "AppImage generation failed, see $OUT/appimage-output.log"

IMAGE="$(ls -1t "$OUT"/*.AppImage | head -1)"
GLIBC="$(find "$APPDIR" -type f \( -name '*.so*' -o -path '*/bin/*' \) | while read -r f; do objdump -T "$f" 2>/dev/null | grep -oE 'GLIBC_2\.[0-9]+(\.[0-9]+)?'; done | sort -uV | tail -1)"
echo ">>> $IMAGE ($(du -h "$IMAGE" | cut -f1)), needs ${GLIBC:-an unknown glibc} or later"
