#!/bin/sh
# runs on every apptainer run / docker run before the pipeline

# if tmpdir is unset or points at a tmpfs, fall back to real disk under $home
if [ -z "$TMPDIR" ] || [ "$(stat -f -c %T "$TMPDIR" 2>/dev/null)" = "tmpfs" ]; then
    export TMPDIR="${HOME:-/tmp}/scratch_tmp"
    mkdir -p "$TMPDIR"
fi

# hand off to whatever command was requested (hyb2-py, hyb2-fold, ...) or the cmd default
exec "$@"