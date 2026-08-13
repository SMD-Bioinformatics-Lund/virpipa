#!/bin/bash

set -euo pipefail

PROJECT_DIR=$(cd "$(dirname "$0")/.." && pwd -P)
TEST_ROOT=$(mktemp -d /tmp/virpipa-deploy-test.XXXXXX)
trap 'rm -r "$TEST_ROOT"' EXIT

SOURCE="$TEST_ROOT/source"
BARE="$TEST_ROOT/bare"
TARGET="$TEST_ROOT/target"
BIN="$TEST_ROOT/bin"
mkdir -p "$SOURCE/scripts" "$SOURCE/workflows" "$BIN"

git init -q "$SOURCE"
git -C "$SOURCE" config user.name 'VirPipa deploy test'
git -C "$SOURCE" config user.email 'virpipa-deploy-test@example.invalid'
printf 'nextflow.enable.dsl=2\n' > "$SOURCE/main.nf"
printf 'profiles { slurm {} hpc {} apptainer {} }\n' > "$SOURCE/nextflow.config"
printf 'workflow HCVPIPE {}\n' > "$SOURCE/workflows/hcvpipe.nf"
printf '#!/bin/bash\n' > "$SOURCE/scripts/helper.sh"
git -C "$SOURCE" add .
git -C "$SOURCE" commit -qm initial
git -C "$SOURCE" tag release/v-test
git init -q --bare "$BARE"
git -C "$SOURCE" push -q "$BARE" HEAD:master --tags

cat > "$BIN/nextflow" <<'EOF'
#!/bin/bash
[[ "$1" == config && "$2" == -profile && "$3" == slurm,hpc,apptainer ]]
EOF
chmod +x "$BIN/nextflow"

deploy() {
    PATH="$BIN:$PATH" VIRPIPA_DEPLOY_REPO="$BARE" VIRPIPA_DEPLOY_ROOT="$TARGET" \
        "$PROJECT_DIR/scripts/deploy_hopper.sh" "$@"
}

deploy --dry-run release/v-test | grep -q 'Would deploy release/v-test'
[[ ! -e "$TARGET/current" ]]

deploy release/v-test
FIRST=$(readlink -f "$TARGET/current")
[[ -f "$FIRST/VERSION" && -f "$FIRST/SHA256SUMS" ]]
[[ "$(basename "$FIRST")" == *-release_v-test ]]
grep -q '^ref=release/v-test$' "$FIRST/VERSION"
grep -q '^describe=release/v-test$' "$FIRST/VERSION"
(cd "$FIRST" && sha256sum -c SHA256SUMS >/dev/null)

printf '# hotfix\n' >> "$FIRST/main.nf"
if deploy master >"$TEST_ROOT/drift.out" 2>&1; then
    echo 'Expected drift detection to abort deployment' >&2
    exit 1
fi
grep -q 'Current release drift detected' "$TEST_ROOT/drift.out"

sleep 1
deploy --force master
SECOND=$(readlink -f "$TARGET/current")
[[ "$FIRST" != "$SECOND" ]]
grep -q '# hotfix' "$FIRST/main.nf"
! grep -q '# hotfix' "$SECOND/main.nf"

touch "$SECOND/untracked-hotfix.txt"
if deploy master >"$TEST_ROOT/unexpected.out" 2>&1; then
    echo 'Expected unexpected-file detection to abort deployment' >&2
    exit 1
fi
grep -q 'Unexpected files' "$TEST_ROOT/unexpected.out"
rm "$SECOND/untracked-hotfix.txt"

exec 8>"$TARGET/.deploy.lock"
flock -n 8
if deploy master >"$TEST_ROOT/lock.out" 2>&1; then
    echo 'Expected deployment lock to reject a concurrent deployment' >&2
    exit 1
fi
grep -q 'Another VirPipa deployment is in progress' "$TEST_ROOT/lock.out"
flock -u 8

if deploy does-not-exist >"$TEST_ROOT/ref.out" 2>&1; then
    echo 'Expected invalid ref to fail' >&2
    exit 1
fi
grep -q 'Git ref not found' "$TEST_ROOT/ref.out"

echo 'Deployment tests passed'
