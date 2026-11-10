#!/bin/bash
set -euo pipefail

# Signs Leonard.app with a Developer ID Application certificate, notarizes it
# with the App Store Connect API, staples the ticket, and re-creates
# Leonard-macOS.zip. Run after package_mac.sh.
#
# Required environment (stored as GitHub Actions secrets):
#   MACOS_CERTIFICATE_P12        base64-encoded .p12 (cert + private key)
#   MACOS_CERTIFICATE_PASSWORD   .p12 export password
#   MACOS_NOTARY_KEY_P8          App Store Connect API key (.p8 file contents)
#   MACOS_NOTARY_KEY_ID          API key ID
#   MACOS_NOTARY_ISSUER_ID       API key issuer ID

APP_NAME="Leonard.app"
ZIP_NAME="Leonard-macOS.zip"
BIN_PATH="$APP_NAME/Contents/MacOS/leonard"
ENTITLEMENTS="scripts/mac/entitlements.plist"

for var in MACOS_CERTIFICATE_P12 MACOS_CERTIFICATE_PASSWORD MACOS_NOTARY_KEY_P8 MACOS_NOTARY_KEY_ID MACOS_NOTARY_ISSUER_ID; do
    if [ -z "${!var:-}" ]; then
        echo "Error: $var is not set"
        exit 1
    fi
done

if [ ! -d "$APP_NAME" ]; then
    echo "Error: $APP_NAME not found; run package_mac.sh first"
    exit 1
fi

WORK_DIR="$(mktemp -d)"
KEYCHAIN="build-$(uuidgen).keychain"
KEYCHAIN_PASSWORD="$(uuidgen)"

cleanup() {
    security delete-keychain "$KEYCHAIN" 2>/dev/null || true
    rm -rf "$WORK_DIR"
}
trap cleanup EXIT

echo "$MACOS_CERTIFICATE_P12" | base64 --decode > "$WORK_DIR/certificate.p12"
printf '%s\n' "$MACOS_NOTARY_KEY_P8" > "$WORK_DIR/AuthKey.p8"

echo "==> Creating temporary keychain..."
security create-keychain -p "$KEYCHAIN_PASSWORD" "$KEYCHAIN"
security set-keychain-settings -lut 21600 "$KEYCHAIN"
security unlock-keychain -p "$KEYCHAIN_PASSWORD" "$KEYCHAIN"
security import "$WORK_DIR/certificate.p12" -k "$KEYCHAIN" \
    -P "$MACOS_CERTIFICATE_PASSWORD" -T /usr/bin/codesign
# Allow codesign to use the key without a GUI prompt on headless runners.
security set-key-partition-list -S apple-tool:,apple:,codesign: -s -k "$KEYCHAIN_PASSWORD" "$KEYCHAIN"
security list-keychains -d user -s "$KEYCHAIN"

IDENTITY="$(security find-identity -v -p codesigning "$KEYCHAIN" | grep 'Developer ID Application' | head -1 | sed 's/.*"\(.*\)"/\1/' || true)"
if [ -z "$IDENTITY" ]; then
    echo "Error: no Developer ID Application identity found in keychain"
    exit 1
fi
echo "==> Signing with identity: $IDENTITY"

echo "==> Signing embedded binary..."
codesign --force --options runtime --timestamp --sign "$IDENTITY" "$BIN_PATH"

echo "==> Signing app bundle (with entitlements)..."
codesign --force --options runtime --timestamp --entitlements "$ENTITLEMENTS" \
    --sign "$IDENTITY" "$APP_NAME"

echo "==> Verifying signature..."
codesign --verify --deep --strict --verbose=2 "$APP_NAME"

echo "==> Creating zip for notarization..."
rm -f "$ZIP_NAME"
ditto -c -k --sequesterRsrc --keepParent "$APP_NAME" "$ZIP_NAME"

echo "==> Submitting to Apple notary service (this can take several minutes)..."
xcrun notarytool submit "$ZIP_NAME" \
    --key "$WORK_DIR/AuthKey.p8" \
    --key-id "$MACOS_NOTARY_KEY_ID" \
    --issuer "$MACOS_NOTARY_ISSUER_ID" \
    --wait

echo "==> Stapling notarization ticket..."
xcrun stapler staple "$APP_NAME"

echo "==> Re-creating distributable zip with stapled app..."
rm -f "$ZIP_NAME"
ditto -c -k --sequesterRsrc --keepParent "$APP_NAME" "$ZIP_NAME"

echo "==> Verifying Gatekeeper assessment..."
spctl -a -vv "$APP_NAME"

echo "==> Successfully created signed and notarized $ZIP_NAME!"
