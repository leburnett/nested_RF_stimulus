#!/usr/bin/env bash
# Launch the Nested RF Stimulus Ephys Dashboard
# Usage: dash-nested   (start the dashboard)
#
# Shell alias (add to ~/.zshrc):
#   alias dash-nested='/Users/burnettl/Documents/GitHub/nested_RF_stimulus/src/dashboard.sh'

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PIXI_DIR="$SCRIPT_DIR/dashboard"
PORT=8051

# Kill any existing process on the port
existing_pid=$(lsof -ti :$PORT 2>/dev/null || true)
if [ -n "$existing_pid" ]; then
    echo "Killing existing process on port $PORT (PID: $existing_pid)..."
    kill $existing_pid 2>/dev/null
    sleep 0.5
fi

echo "Starting dashboard on http://localhost:$PORT ..."

# Open browser once server is ready (background)
(
    for i in $(seq 1 30); do
        sleep 1
        if curl -s -o /dev/null http://localhost:$PORT/ 2>/dev/null; then
            open "http://localhost:$PORT/"
            break
        fi
    done
) &

pixi run -e default --manifest-path "$PIXI_DIR/pixi.toml" dashboard
