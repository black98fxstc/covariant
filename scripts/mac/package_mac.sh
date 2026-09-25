#!/bin/bash
set -e

# Usage: ./package_mac.sh [path-to-build-folder]
BUILD_DIR="${1:-build}"
BIN_PATH="${BUILD_DIR}/leonard"
APP_NAME="Leonard.app"
ZIP_NAME="Leonard-macOS.zip"

if [ ! -f "$BIN_PATH" ]; then
    echo "Error: leonard binary not found at $BIN_PATH"
    exit 1
fi

echo "==> Compiling AppleScript droplet into $APP_NAME..."
osacompile -o "$APP_NAME" scripts/mac/droplet.applescript

echo "==> Embedding leonard binary into $APP_NAME/Contents/MacOS/..."
cp "$BIN_PATH" "$APP_NAME/Contents/MacOS/leonard"
chmod +x "$APP_NAME/Contents/MacOS/leonard"

echo "==> Configuring Info.plist for .wsp file association..."
PLIST="$APP_NAME/Contents/Info.plist"

# Add Document Types for .wsp files
/usr/libexec/PlistBuddy -c "Add :CFBundleDocumentTypes list" "$PLIST" 2>/dev/null || true
/usr/libexec/PlistBuddy -c "Add :CFBundleDocumentTypes:0:CFBundleTypeExtensions array" "$PLIST" 2>/dev/null || true
/usr/libexec/PlistBuddy -c "Add :CFBundleDocumentTypes:0:CFBundleTypeExtensions:0 string wsp" "$PLIST" 2>/dev/null || true
/usr/libexec/PlistBuddy -c "Add :CFBundleDocumentTypes:0:CFBundleTypeName string FlowJo Workspace" "$PLIST" 2>/dev/null || true
/usr/libexec/PlistBuddy -c "Add :CFBundleDocumentTypes:0:CFBundleTypeRole string Viewer" "$PLIST" 2>/dev/null || true

# Set clean Bundle ID and Version
/usr/libexec/PlistBuddy -c "Set :CFBundleIdentifier org.covariant.leonard" "$PLIST" 2>/dev/null || \
/usr/libexec/PlistBuddy -c "Add :CFBundleIdentifier string org.covariant.leonard" "$PLIST" 2>/dev/null || true

/usr/libexec/PlistBuddy -c "Set :CFBundleShortVersionString 0.1.0" "$PLIST" 2>/dev/null || \
/usr/libexec/PlistBuddy -c "Add :CFBundleShortVersionString string 0.1.0" "$PLIST" 2>/dev/null || true

echo "==> Packaging into $ZIP_NAME..."
zip -r -y "$ZIP_NAME" "$APP_NAME"

echo "==> Successfully created $ZIP_NAME!"
