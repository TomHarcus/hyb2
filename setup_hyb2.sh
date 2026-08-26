#!/usr/bin/env bash
set -e

# optional arg for if user uses apptainer (cluster only)
SIF="${1:-/exports/eddie/scratch/$USER/hyb2_latest.sif}"

# detect the shell type
case "$(basename "${SHELL:-bash}")" in
    zsh) RC="$HOME/.zshrc" ;;
    *)   RC="$HOME/.bashrc" ;;
esac

# pick which method is installed (prefer docker if both available)
if command -v docker >/dev/null 2>&1; then
    METHOD=docker
    if command -v apptainer >/dev/null 2>&1; then
        echo "both docker and apptainer found. using docker"
    fi
elif command -v apptainer >/dev/null 2>&1; then
    METHOD=apptainer
else
    echo "error: need docker or apptainer installed" >&2
    exit 1
fi

# write the hyb2() wrapper to its own file (overwrites each time)
if [ "$METHOD" = docker ]; then
    # quoted text so every $ stays literal
    cat > "$HOME/.hyb2.sh" <<'EOF'
hyb2() {
    local tty=""; [ -t 1 ] && tty="-t"
    docker run --rm $tty --user "$(id -u):$(id -g)" -e HOME=/tmp \
        -v "$PWD:/data" -w /data ghcr.io/tomharcus/hyb2:latest hyb2 "$@"
}
EOF
else
    # unqouted text, only SIF is expanded
    cat > "$HOME/.hyb2.sh" <<EOF
hyb2() {
    apptainer run \${TMPDIR:+--bind "\$TMPDIR"} --bind /exports/eddie/scratch/\$USER "$SIF" hyb2 "\$@"
}
EOF
fi

# source it from the rc file
LINE='[ -f ~/.hyb2.sh ] && source ~/.hyb2.sh'
grep -qxF "$LINE" "$RC" 2>/dev/null || echo "$LINE" >> "$RC"

echo "installed hyb2() into ~/.hyb2.sh (method: $METHOD)"
[ "$METHOD" = apptainer ] && echo "using SIF: $SIF"
echo "sourced from: $RC"
echo "open new terminal, or run: source $RC"



