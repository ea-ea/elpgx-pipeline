#!/bin/bash
#
# ElPGx Build Script
# Builds standalone executables for macOS, Windows, and Linux
#
# Usage:
#   ./scripts/build_executable.sh          # Build for current platform
#   ./scripts/build_executable.sh --clean  # Clean and rebuild
#
# Requirements:
#   - Python 3.10+
#   - PyInstaller 6.0+
#   - Required Python packages
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="$PROJECT_ROOT/build"
DIST_DIR="$PROJECT_ROOT/dist"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}"
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║              ElPGx Build Script v1.0.0                            ║"
echo "║          Building Standalone Desktop Application                  ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

# Check Python version
PYTHON_VERSION=$(python3 --version 2>&1 | cut -d' ' -f2 | cut -d'.' -f1,2)
echo -e "${BLUE}[INFO]${NC} Python version: $PYTHON_VERSION"

if [[ $(echo "$PYTHON_VERSION < 3.10" | bc -l) -eq 1 ]]; then
    echo -e "${RED}[ERROR]${NC} Python 3.10 or higher is required"
    exit 1
fi

# Clean previous builds if requested
if [[ "$1" == "--clean" ]]; then
    echo -e "${YELLOW}[INFO]${NC} Cleaning previous builds..."
    rm -rf "$BUILD_DIR"/*.build
    rm -rf "$BUILD_DIR"/*.dist
    rm -rf "$DIST_DIR"
fi

# Create build directory
mkdir -p "$BUILD_DIR"
mkdir -p "$DIST_DIR"

# Install/upgrade PyInstaller
echo -e "${BLUE}[INFO]${NC} Checking PyInstaller..."
pip3 install --upgrade pyinstaller >/dev/null 2>&1

# Install project dependencies
echo -e "${BLUE}[INFO]${NC} Installing project dependencies..."
pip3 install -e "$PROJECT_ROOT[build]" >/dev/null 2>&1

# Create assets directory if missing
mkdir -p "$PROJECT_ROOT/assets"

# Build using spec file
echo -e "${BLUE}[INFO]${NC} Building standalone executable..."
cd "$PROJECT_ROOT"

pyinstaller build/elpgx.spec \
    --clean \
    --noconfirm \
    --distpath "$DIST_DIR" \
    --workpath "$BUILD_DIR"

# Check build success
if [[ $? -eq 0 ]]; then
    echo -e "${GREEN}"
    echo "╔═══════════════════════════════════════════════════════════════════╗"
    echo "║                    BUILD SUCCESSFUL!                              ║"
    echo "╚═══════════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"

    echo -e "${GREEN}[SUCCESS]${NC} Executables created in: $DIST_DIR"
    echo ""
    echo "Files created:"
    ls -la "$DIST_DIR/ElPGx/"

    # Platform-specific instructions
    if [[ "$OSTYPE" == "darwin"* ]]; then
        echo ""
        echo -e "${BLUE}macOS Usage:${NC}"
        echo "  GUI:  open $DIST_DIR/ElPGx.app"
        echo "  CLI:  $DIST_DIR/ElPGx/elpgx --help"
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        echo ""
        echo -e "${BLUE}Linux Usage:${NC}"
        echo "  GUI:  $DIST_DIR/ElPGx/ElPGx"
        echo "  CLI:  $DIST_DIR/ElPGx/elpgx --help"
    elif [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "win32" ]]; then
        echo ""
        echo -e "${BLUE}Windows Usage:${NC}"
        echo "  GUI:  $DIST_DIR\\ElPGx\\ElPGx.exe"
        echo "  CLI:  $DIST_DIR\\ElPGx\\elpgx.exe --help"
    fi
else
    echo -e "${RED}[ERROR]${NC} Build failed!"
    exit 1
fi
