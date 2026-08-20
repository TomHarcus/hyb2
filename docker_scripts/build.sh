#!/usr/bin/env bash
set -euo pipefail

IMAGE=hyb2
TAG=latest
ARCHIVE=hyb2.tar.gz

case "${1:-help}" in
  build)
    # build docker on Fedora. preflight in the Dockerfile fails here if broken
    ROOT="$(cd "$(dirname "$0")/.." && pwd)"
    docker build -t "$IMAGE:$TAG" "$ROOT"
    ;;

  smoke)
    # no data needed. proves entrypoint, env, and every tool resolve inside the image
    docker run --rm "$IMAGE:$TAG" hyb2-py --help
    docker run --rm "$IMAGE:$TAG" bash -lc \
      'RNAfold --version && bowtie2 --version | head -1 \
       && Rscript -e "suppressMessages(library(DESeq2))" && echo "tools OK"'
    ;;

  run)
    data="$(cd "$2" && pwd)"; shift 2

    # tqdm progress ui
    TTY=""; [ -t 1 ] && TTY="-t"

    # varna interactive gui
    GUI=""
    if [ -n "${DISPLAY:-}" ] && [ -d /tmp/.X11-unix ]; then
        xhost +local: >/dev/null 2>&1 || true
        GUI="-e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix"
    fi

    docker run --rm $TTY $GUI \
    --user "$(id -u):$(id -g)" -e HOME=/tmp \
    -v "$data:/data" -w /data "$IMAGE:$TAG" "$@"
    # TMPDIR is handled by entrypoint.sh (falls back to real disk under $HOME);
    # for a 50GB run, add:  -e TMPDIR=/data/tmp   and mkdir it under $data first
    ;;

  save)
    # freeze the image to a tarball to scp to Eddie 
    docker save "$IMAGE:$TAG" | gzip > "$ARCHIVE"
    echo "wrote $ARCHIVE ($(du -h "$ARCHIVE" | cut -f1))"
    ;;

  eddie)
    cat <<'EOF'
# on Eddie, after scp of hyb2.tar.gz 
module load apptainer/1.4.4
gunzip -f hyb2.tar.gz
apptainer build hyb2.sif docker-archive://hyb2.tar     # convert once

# run.  --bind "$TMPDIR" makes node-local scratch writable inside (big sort spills there)
# inside an SGE job $TMPDIR is already /local/$JOB_ID (real disk); interactively, set it.
apptainer run \
    --bind "$TMPDIR" \
    --bind /exports/eddie/scratch/<you>/mydata:/data \
    hyb2.sif hyb2-py -i /data/reads.sam -d /data/ref.fasta -o run -a MyRNA -x 3900 -l 300 -r cplfold

EOF
    ;;

  *)
    echo "usage: ./build.sh {build | smoke | run <DATA_DIR> <args> | save | eddie}"
    ;;
esac