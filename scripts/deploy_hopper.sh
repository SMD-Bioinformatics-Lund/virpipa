#!/bin/bash

set -euo pipefail

SOURCE_REPO="${VIRPIPA_DEPLOY_REPO:-/fs1/jonas/bare/virpipa}"
DEPLOY_ROOT="${VIRPIPA_DEPLOY_ROOT:-/fs1/pipelines/virpipa}"
DRY_RUN=false
FORCE=false

usage() {
    echo "Usage: $0 [--dry-run] [--force] <branch|tag|commit>" >&2
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true ;;
        --force) FORCE=true ;;
        -h|--help) usage; exit 0 ;;
        --*) die "Unknown option: $1" ;;
        *)
            [[ -z "${SOURCE_REF:-}" ]] || die "Only one Git ref may be specified"
            SOURCE_REF="$1"
            ;;
    esac
    shift
done

[[ -n "${SOURCE_REF:-}" ]] || { usage; exit 1; }
[[ -d "$SOURCE_REPO" ]] || die "Bare repository not found: $SOURCE_REPO"

COMMIT=$(git --git-dir="$SOURCE_REPO" rev-parse --verify "${SOURCE_REF}^{commit}" 2>/dev/null) || die "Git ref not found: $SOURCE_REF"
SHORT_COMMIT=$(git --git-dir="$SOURCE_REPO" rev-parse --short=12 "$COMMIT")
DESCRIBE=$(git --git-dir="$SOURCE_REPO" describe --tags --always "$COMMIT" 2>/dev/null || echo "$SHORT_COMMIT")
SAFE_DESCRIBE=$(printf '%s' "$DESCRIBE" | sed -e 's#[^A-Za-z0-9._-]#_#g' -e 's#^\.*##')
[[ -n "$SAFE_DESCRIBE" ]] || SAFE_DESCRIBE="$SHORT_COMMIT"
COMMIT_TIME=$(git --git-dir="$SOURCE_REPO" show -s --format=%ci "$COMMIT")
DEPLOY_TIME=$(date -u +%Y-%m-%dT%H:%M:%SZ)
RELEASE_ID="$(date -u +%Y%m%dT%H%M%SZ)-${SAFE_DESCRIBE}"
RELEASES_DIR="$DEPLOY_ROOT/releases"
RELEASE_DIR="$RELEASES_DIR/$RELEASE_ID"

check_current_release() {
    local current_link="$DEPLOY_ROOT/current"
    local current_dir checksum_status unexpected_status newer_status

    [[ -e "$current_link" || -L "$current_link" ]] || return 0
    [[ -L "$current_link" ]] || die "$current_link exists but is not a symlink"
    current_dir=$(readlink -f "$current_link") || die "Cannot resolve $current_link"
    [[ -d "$current_dir" ]] || die "Current release does not exist: $current_dir"
    [[ -f "$current_dir/VERSION" && -f "$current_dir/SHA256SUMS" ]] || die "Current release lacks VERSION or SHA256SUMS: $current_dir"

    checksum_status=$(cd "$current_dir" && { sha256sum -c SHA256SUMS 2>&1 || true; } | awk '$NF != "OK" { print }')
    unexpected_status=$(cd "$current_dir" && {
        find . -type f -o -type l
        echo ./SHA256SUMS
    } | sed 's#^\./##' | sort -u | comm -13 <(awk '{print $2}' SHA256SUMS | sed -e 's#^\*##' -e 's#^\./##' | { cat; echo SHA256SUMS; } | sort -u) -)
    newer_status=$(find "$current_dir" \( -type f -o -type l \) -newer "$current_dir/VERSION" ! -name VERSION ! -name SHA256SUMS -printf '%P\n' | sort)

    if [[ -n "$checksum_status" || -n "$unexpected_status" || -n "$newer_status" ]]; then
        echo "Current release drift detected: $current_dir" >&2
        [[ -z "$checksum_status" ]] || { echo "Checksum differences:" >&2; echo "$checksum_status" >&2; }
        [[ -z "$unexpected_status" ]] || { echo "Unexpected files:" >&2; echo "$unexpected_status" >&2; }
        [[ -z "$newer_status" ]] || { echo "Files newer than VERSION:" >&2; echo "$newer_status" >&2; }
        $FORCE || die "Refusing to deploy over detected drift; review it or rerun with --force"
        echo "WARNING: continuing because --force was supplied; the modified release will be preserved" >&2
    fi
}

if $DRY_RUN; then
    check_current_release
    echo "Would deploy $SOURCE_REF ($COMMIT)"
    echo "Source: $SOURCE_REPO"
    echo "Release: $RELEASE_DIR"
    echo "Would atomically update: $DEPLOY_ROOT/current"
    exit 0
fi

mkdir -p "$RELEASES_DIR"
exec 9>"$DEPLOY_ROOT/.deploy.lock"
flock -n 9 || die "Another VirPipa deployment is in progress"

check_current_release
[[ ! -e "$RELEASE_DIR" ]] || die "Release already exists: $RELEASE_DIR"

STAGING_DIR=$(mktemp -d "$RELEASES_DIR/.staging.${RELEASE_ID}.XXXXXX")
cleanup() {
    [[ -z "${STAGING_DIR:-}" || ! -d "$STAGING_DIR" ]] || rm -r "$STAGING_DIR"
    [[ -z "${VALIDATION_DIR:-}" || ! -d "$VALIDATION_DIR" ]] || rm -r "$VALIDATION_DIR"
    [[ -z "${LINK_TMP:-}" || ! -L "$LINK_TMP" ]] || rm "$LINK_TMP"
}
trap cleanup EXIT

git --git-dir="$SOURCE_REPO" archive "$COMMIT" | tar -x -C "$STAGING_DIR"
for required in main.nf nextflow.config workflows/hcvpipe.nf scripts; do
    [[ -e "$STAGING_DIR/$required" ]] || die "Candidate release is missing $required"
done

cat > "$STAGING_DIR/VERSION" <<EOF
ref=$SOURCE_REF
commit=$COMMIT
describe=$DESCRIBE
commit_time=$COMMIT_TIME
deployed_at=$DEPLOY_TIME
deployed_by=$(id -un)
deployed_from=$(hostname -f 2>/dev/null || hostname)
EOF

if ! command -v nextflow >/dev/null 2>&1; then
    set +u
    source /etc/profile
    module load Java/23.0.2 nextflow/26.04.3 singularity/3.8.0
    set -u
fi
VALIDATION_DIR=$(mktemp -d /tmp/virpipa-deploy-validation.XXXXXX)
(cd "$VALIDATION_DIR" && nextflow config "$STAGING_DIR" -profile slurm,hpc,apptainer >/dev/null)
rm -r "$VALIDATION_DIR"
VALIDATION_DIR=''

(cd "$STAGING_DIR" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
mv "$STAGING_DIR" "$RELEASE_DIR"
STAGING_DIR=''

LINK_TMP="$DEPLOY_ROOT/.current.${RELEASE_ID}"
ln -s "releases/$RELEASE_ID" "$LINK_TMP"
mv -Tf "$LINK_TMP" "$DEPLOY_ROOT/current"
LINK_TMP=''

printf '%s\t%s\t%s\t%s\t%s\n' "$DEPLOY_TIME" "$(id -un)" "$SOURCE_REF" "$COMMIT" "$RELEASE_ID" >> "$DEPLOY_ROOT/deployments.log"
echo "Deployed $SOURCE_REF ($COMMIT)"
echo "Current release: $RELEASE_DIR"
