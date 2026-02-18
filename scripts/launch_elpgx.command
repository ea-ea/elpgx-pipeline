#!/bin/bash
# ElPGx Launcher for macOS
# Double-click this file to launch the ElPGx GUI

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
APP_PATH="$SCRIPT_DIR/../dist/ElPGx.app"

if [ -d "$APP_PATH" ]; then
    open "$APP_PATH"
else
    echo "ElPGx.app not found at $APP_PATH"
    echo "Please build the application first using:"
    echo "  pyinstaller build/elpgx.spec --clean"
    read -p "Press Enter to close..."
fi
