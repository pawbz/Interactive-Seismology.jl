#!/usr/bin/env bash
# Build and serve Interactive Seismology directly at IP:port, without nginx.
#
# Required: PUBLIC_HOST=203.0.113.10 ./serve_remote.sh start
# Optional: PUBLIC_LIVE_URL=https://notebooks.example.edu ./serve_remote.sh start

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SITE_DIR="${SITE_DIR:-$REPO_DIR/_site}"
RUNTIME_DIR="${RUNTIME_DIR:-$REPO_DIR/.remote-server}"
SITE_HOST="${SITE_HOST:-0.0.0.0}"
SITE_PORT="${SITE_PORT:-8000}"
LIVE_HOST="${LIVE_HOST:-0.0.0.0}"
LIVE_PORT="${LIVE_PORT:-1234}"
PUBLIC_HOST="${PUBLIC_HOST:-}"
PUBLIC_SCHEME="${PUBLIC_SCHEME:-http}"
JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$HOME/.julia}"

SITE_PID_FILE="$RUNTIME_DIR/site.pid"
LIVE_PID_FILE="$RUNTIME_DIR/live-notebooks.pid"
SITE_LOG="$RUNTIME_DIR/site.log"
LIVE_LOG="$RUNTIME_DIR/live-notebooks.log"

usage() {
    cat <<'EOF'
Usage: PUBLIC_HOST=<public-IP-or-hostname> ./serve_remote.sh <command>

Commands:
  start       Build the site and start the static and live services in the background.
  stop        Stop services started by this script.
  restart     Stop, rebuild, and start the services.
  status      Report whether the two services are running.
  foreground  Build the site, then run both services until Ctrl+C.

Environment overrides:
  PUBLIC_HOST       Required unless PUBLIC_LIVE_URL is supplied; e.g. 203.0.113.10
  PUBLIC_LIVE_URL   URL embedded in live notebooks; defaults to http://PUBLIC_HOST:1234
  SITE_PORT         Static site port (default: 8000)
  LIVE_PORT         Pluto live-widget port (default: 1234)
  SITE_HOST         Static bind address (default: 0.0.0.0)
  LIVE_HOST         Pluto bind address (default: 0.0.0.0)
  JULIA_DEPOT_PATH  Julia package depot (default: $HOME/.julia)

Open SITE_PORT and LIVE_PORT in the VM firewall. Browse only the static site URL,
for example http://203.0.113.10:8000/; port 1234 is the live-widget API.
EOF
}

fail() {
    echo "Error: $*" >&2
    exit 1
}

validate_port() {
    local port="$1"
    [[ "$port" =~ ^[0-9]+$ ]] && ((port >= 1 && port <= 65535)) || fail "invalid port: $port"
}

pid_is_running() {
    local pid_file="$1"
    local expected_command="$2"
    [[ -f "$pid_file" ]] || return 1
    local pid
    pid="$(<"$pid_file")"
    [[ "$pid" =~ ^[0-9]+$ ]] &&
        kill -0 "$pid" 2>/dev/null &&
        ps -p "$pid" -o command= 2>/dev/null | grep -Fq "$expected_command"
}

stop_service() {
    local name="$1"
    local pid_file="$2"
    local expected_command="$3"
    if ! [[ -f "$pid_file" ]]; then
        return
    fi

    local pid
    pid="$(<"$pid_file")"
    if pid_is_running "$pid_file" "$expected_command"; then
        echo "Stopping $name (PID $pid)..."
        kill "$pid"
    elif [[ "$pid" =~ ^[0-9]+$ ]] && kill -0 "$pid" 2>/dev/null; then
        fail "refusing to stop PID $pid: it is not the expected $name process"
    fi
    rm -f "$pid_file"
}

build_site() {
    local public_live_url="${PUBLIC_LIVE_URL:-}"
    if [[ -z "$public_live_url" ]]; then
        [[ -n "$PUBLIC_HOST" ]] || fail "set PUBLIC_HOST or PUBLIC_LIVE_URL before starting"
        public_live_url="${PUBLIC_SCHEME}://${PUBLIC_HOST}:${LIVE_PORT}"
    fi

    echo "Building the site with live widgets at $public_live_url ..."
    JULIA_DEPOT_PATH="$JULIA_DEPOT_PATH" \
    PLUTO_SLIDER_SERVER_URL="$public_live_url" \
        julia --project=pluto-deployment-environment -e 'import Pkg; Pkg.instantiate(); include("generate.jl")'
}

start_background() {
    mkdir -p "$RUNTIME_DIR"
    if pid_is_running "$SITE_PID_FILE" "http.server" || pid_is_running "$LIVE_PID_FILE" "serve_live_notebooks.jl"; then
        fail "a service is already running; use ./serve_remote.sh status or restart"
    fi

    rm -f "$SITE_PID_FILE" "$LIVE_PID_FILE"
    echo "Starting static site at http://${PUBLIC_HOST:-$SITE_HOST}:$SITE_PORT/ ..."
    nohup python3 -m http.server "$SITE_PORT" --bind "$SITE_HOST" --directory "$SITE_DIR" \
        >"$SITE_LOG" 2>&1 &
    echo $! >"$SITE_PID_FILE"

    echo "Starting live-widget API on $LIVE_HOST:$LIVE_PORT ..."
    nohup env JULIA_DEPOT_PATH="$JULIA_DEPOT_PATH" \
        julia --project=pluto-deployment-environment serve_live_notebooks.jl \
        --host "$LIVE_HOST" --port "$LIVE_PORT" >"$LIVE_LOG" 2>&1 &
    echo $! >"$LIVE_PID_FILE"

    sleep 0.2
    if ! pid_is_running "$SITE_PID_FILE" "http.server" || ! pid_is_running "$LIVE_PID_FILE" "serve_live_notebooks.jl"; then
        stop_service "static site" "$SITE_PID_FILE" "http.server"
        stop_service "live-widget API" "$LIVE_PID_FILE" "serve_live_notebooks.jl"
        fail "a service failed to start; inspect $SITE_LOG and $LIVE_LOG"
    fi

    echo "Started. Logs: $SITE_LOG and $LIVE_LOG"
    echo "Open: http://${PUBLIC_HOST:-$SITE_HOST}:$SITE_PORT/"
}

run_foreground() {
    mkdir -p "$RUNTIME_DIR"
    python3 -m http.server "$SITE_PORT" --bind "$SITE_HOST" --directory "$SITE_DIR" \
        >"$SITE_LOG" 2>&1 &
    local site_pid=$!
    trap 'kill "$site_pid" 2>/dev/null || true' EXIT INT TERM

    echo "Static site: http://${PUBLIC_HOST:-$SITE_HOST}:$SITE_PORT/"
    echo "Press Ctrl+C to stop both services."
    JULIA_DEPOT_PATH="$JULIA_DEPOT_PATH" \
        julia --project=pluto-deployment-environment serve_live_notebooks.jl \
        --host "$LIVE_HOST" --port "$LIVE_PORT"
}

main() {
    local command="${1:-}"
    validate_port "$SITE_PORT"
    validate_port "$LIVE_PORT"
    [[ "$SITE_PORT" != "$LIVE_PORT" ]] || fail "SITE_PORT and LIVE_PORT must be different"
    cd "$REPO_DIR"

    case "$command" in
        start)
            build_site
            start_background
            ;;
        stop)
            stop_service "static site" "$SITE_PID_FILE" "http.server"
            stop_service "live-widget API" "$LIVE_PID_FILE" "serve_live_notebooks.jl"
            ;;
        restart)
            stop_service "static site" "$SITE_PID_FILE" "http.server"
            stop_service "live-widget API" "$LIVE_PID_FILE" "serve_live_notebooks.jl"
            build_site
            start_background
            ;;
        status)
            pid_is_running "$SITE_PID_FILE" "http.server" && echo "static site: running (PID $(<"$SITE_PID_FILE"))" || echo "static site: stopped"
            pid_is_running "$LIVE_PID_FILE" "serve_live_notebooks.jl" && echo "live-widget API: running (PID $(<"$LIVE_PID_FILE"))" || echo "live-widget API: stopped"
            ;;
        foreground)
            build_site
            run_foreground
            ;;
        -h|--help|help|"")
            usage
            ;;
        *)
            usage >&2
            fail "unknown command: $command"
            ;;
    esac
}

main "$@"
