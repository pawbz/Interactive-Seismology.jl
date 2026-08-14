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
LIVE_PORT="${LIVE_PORT:-4321}"
PUBLIC_HOST="${PUBLIC_HOST:-}"
PUBLIC_SCHEME="${PUBLIC_SCHEME:-http}"
JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$HOME/.julia}"

SITE_PID_FILE="$RUNTIME_DIR/site.pid"
LIVE_PID_FILE="$RUNTIME_DIR/live-notebooks.pid"
SITE_LOG="$RUNTIME_DIR/site.log"
LIVE_LOG="$RUNTIME_DIR/live-notebooks.log"
LIVE_NOTEBOOKS_FILE="$RUNTIME_DIR/live-notebooks.list"

usage() {
    cat <<'EOF'
Usage: PUBLIC_HOST=<public-IP-or-hostname> ./serve_remote.sh <command>

Commands:
  start       Build the site and start the static and live services in the background.
  stop        Stop services started by this script.
  restart     Stop, rebuild, and start the services.
  status      Report service health plus host and live-widget memory use.
  foreground  Build the site, then run both services until Ctrl+C.

Environment overrides:
  PUBLIC_HOST       Required unless PUBLIC_LIVE_URL is supplied; e.g. 203.0.113.10
  PUBLIC_LIVE_URL   URL embedded in live notebooks; defaults to http://PUBLIC_HOST:4321
  SITE_PORT         Static site port (default: 8000)
  LIVE_PORT         Pluto live-widget port (default: 4321)
  SITE_HOST         Static bind address (default: 0.0.0.0)
  LIVE_HOST         Pluto bind address (default: 0.0.0.0)
  JULIA_DEPOT_PATH  Julia package depot (default: $HOME/.julia)

Open SITE_PORT and LIVE_PORT in the VM firewall. Browse only the static site URL,
for example http://203.0.113.10:8000/; port 4321 is the live-widget API.
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

format_memory() {
    local kib="${1:-0}"
    awk -v kib="$kib" 'BEGIN {
        if (kib >= 1048576) printf "%.2f GiB", kib / 1048576
        else if (kib >= 1024) printf "%.1f MiB", kib / 1024
        else printf "%d KiB", kib
    }'
}

process_rss_kib() {
    ps -o rss= -p "$1" 2>/dev/null | tr -d '[:space:]' || true
}

process_tree_pids() {
    local parent="$1"
    local child
    printf '%s\n' "$parent"
    command -v pgrep >/dev/null 2>&1 || return 0
    while IFS= read -r child; do
        [[ "$child" =~ ^[0-9]+$ ]] || continue
        process_tree_pids "$child"
    done < <(pgrep -P "$parent" 2>/dev/null || true)
}

process_cwd() {
    local pid="$1"
    if [[ -e "/proc/$pid/cwd" ]]; then
        readlink -f "/proc/$pid/cwd" 2>/dev/null || true
    fi
    return 0
}

show_host_memory() {
    [[ -r /proc/meminfo ]] || return 0
    local total_kib available_kib used_kib
    read -r total_kib available_kib < <(awk '
        /^MemTotal:/ { total = $2 }
        /^MemAvailable:/ { available = $2 }
        END { print total, available }
    ' /proc/meminfo)
    [[ "$total_kib" =~ ^[0-9]+$ && "$available_kib" =~ ^[0-9]+$ ]] || return 0
    used_kib=$((total_kib - available_kib))
    echo "host RAM: $(format_memory "$used_kib") used / $(format_memory "$total_kib") total ($(format_memory "$available_kib") available)"
}

show_live_memory() {
    local live_pid
    live_pid="$(<"$LIVE_PID_FILE")"

    local -a pids=()
    local pid rss total_kib=0
    while IFS= read -r pid; do
        [[ "$pid" =~ ^[0-9]+$ ]] && pids+=("$pid")
    done < <(process_tree_pids "$live_pid")

    for pid in "${pids[@]}"; do
        rss="$(process_rss_kib "$pid")"
        [[ "$rss" =~ ^[0-9]+$ ]] && total_kib=$((total_kib + rss))
    done

    echo "live-widget RAM (approximate RSS, process tree): $(format_memory "$total_kib") across ${#pids[@]} process(es)"
    [[ -f "$LIVE_NOTEBOOKS_FILE" ]] || {
        echo "  Per-notebook RAM will be available after the next live-widget restart."
        return
    }

    # Pluto workers use the notebook folder as their current directory. This is
    # exact for folders with one live notebook; shared folders stay grouped.
    local -a groups=()
    local notebook notebook_dir group cwd
    declare -A group_notebooks=() group_rss=() group_processes=()
    while IFS= read -r notebook; do
        [[ -n "$notebook" ]] || continue
        notebook_dir="$(dirname "$notebook")"
        group="$REPO_DIR/src/$notebook_dir"
        if [[ -n "${group_notebooks[$group]+x}" ]]; then
            group_notebooks["$group"]+="|$notebook"
        else
            groups+=("$group")
            group_notebooks["$group"]="$notebook"
        fi
    done < "$LIVE_NOTEBOOKS_FILE"

    local unassigned_kib=0 unassigned_processes=0
    for pid in "${pids[@]}"; do
        rss="$(process_rss_kib "$pid")"
        rss="${rss:-0}"
        cwd="$(process_cwd "$pid")"
        if [[ -n "$cwd" && -n "${group_notebooks[$cwd]+x}" ]]; then
            group_rss["$cwd"]=$(( ${group_rss[$cwd]:-0} + rss ))
            group_processes["$cwd"]=$(( ${group_processes[$cwd]:-0} + 1 ))
        else
            unassigned_kib=$((unassigned_kib + rss))
            unassigned_processes=$((unassigned_processes + 1))
        fi
    done

    echo "live notebook workers:"
    local label members
    for group in "${groups[@]}"; do
        members="${group_notebooks[$group]}"
        if [[ "$members" == *"|"* ]]; then
            label="$(basename "$group")/ (shared by ${members//|/, })"
        else
            label="$members"
        fi
        printf '  %s: %s (%s process(es))\n' \
            "$label" \
            "$(format_memory "${group_rss[$group]:-0}")" \
            "${group_processes[$group]:-0}"
    done
    printf '  server/other live processes: %s (%s process(es))\n' \
        "$(format_memory "$unassigned_kib")" "$unassigned_processes"
}

stop_service() {
    local name="$1"
    local pid_file="$2"
    local expected_command="$3"
    if ! [[ -f "$pid_file" ]]; then
        return 0
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
    if [[ "$pid_file" == "$LIVE_PID_FILE" ]]; then
        rm -f "$LIVE_NOTEBOOKS_FILE"
    fi
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

    rm -f "$SITE_PID_FILE" "$LIVE_PID_FILE" "$LIVE_NOTEBOOKS_FILE"
    echo "Starting static site at http://${PUBLIC_HOST:-$SITE_HOST}:$SITE_PORT/ ..."
    nohup python3 -m http.server "$SITE_PORT" --bind "$SITE_HOST" --directory "$SITE_DIR" \
        >"$SITE_LOG" 2>&1 &
    echo $! >"$SITE_PID_FILE"

    echo "Starting live-widget API on $LIVE_HOST:$LIVE_PORT ..."
    nohup env JULIA_DEPOT_PATH="$JULIA_DEPOT_PATH" \
        LIVE_NOTEBOOKS_FILE="$LIVE_NOTEBOOKS_FILE" \
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
    LIVE_NOTEBOOKS_FILE="$LIVE_NOTEBOOKS_FILE" \
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
            show_host_memory
            if pid_is_running "$LIVE_PID_FILE" "serve_live_notebooks.jl"; then
                echo "live-widget API: running (PID $(<"$LIVE_PID_FILE"))"
                show_live_memory
            else
                echo "live-widget API: stopped"
            fi
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
