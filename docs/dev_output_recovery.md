# Recover completed development runs

The November 11, 2026 campaign completed 240 analyses across 40 run directories
under `/fs1/results_dev/virpipa/hcv`. Its markers were published to production
and pointed at `/access/virpipa/hcv`, rather than those development results.
The analyses do not need to be rerun to correct the handoff.

## Evidence and prepared markers

The local recovery bundle is `/tmp/virpipa-dev-recovery-20261111/`. It contains:

- `production-markers/`: byte-verified copies of the 40 November 11 production
  `.sqlimport.done` / `.sqlimport.error` files, retaining their timestamps.
- `importer-evidence.txt`: relevant lines from Lennart's
  `/data/bnf/scripts/middleman/middleman.done.log` and `middleman.error.log`,
  plus the runner's status handling and verified dev environment paths.
- `manifest.json`: each run's log outcome, original marker checksum, verified
  dev target and aggregate QC record count (240 records total).
- `dev-markers/`: 40 corrected `<run>.sqlimport` files, each containing exactly
  `--run-dir /fs1/results_dev/virpipa/hcv/<run>/` followed by a newline.
- `run-dirs.txt` and `SHA256SUMS`: manual import targets and bundle checksums.

Preserve this temporary bundle in the campaign archive before `/tmp` is cleaned.
The original production markers were left in place and unchanged.

The importer logs record 36 accepted tasks and these four failures:

```text
241004_A01932_0484_BHJLNLDRX5+240910_A01932_0474_AHFNKHDRX5
250328_A01932_0567_AHVLT7DRX5+250128_A01932_0538_AHN7V3DRX5
260217_A00681_1211_AHJMLTDRX7
260605_A01932_0704_BHFTGCDMX2+260512_A00681_1243_BHMH3CDRX7
```

These outcomes come from log entries at 17:00–17:05, not marker suffixes alone.
Middleman accepts both zero exit status and recognized `already loaded` /
`DuplicateKeyError` messages. The logs do not retain importer stdout, exit codes
or failure stderr, so they establish task outcomes, not row changes or failure
causes. No database rollback or production re-import was performed.

## Manual development import

Copy the bundle from the laptop to Lennart, then verify it:

```bash
scp -r /tmp/virpipa-dev-recovery-20261111 lennart:/tmp/
```

On Lennart, use the existing dev checkout, environment and config. These paths
were checked during recovery; the dev results root is also visible there.
The config names `data_dev/virtitta.sqlite3` and the result root
`/fs1/results_dev/virpipa/hcv`. Recheck these before importing.

```bash
bundle=/tmp/virpipa-dev-recovery-20261111
(cd "$bundle" && sha256sum -c SHA256SUMS) || exit 1
export VIRTITTA_SRC=/data/bnf/dev/jonas/hcv/virtitta-dev-setup
export VIRTITTA_DEV_CONFIG="$VIRTITTA_SRC/virtitta.dev.toml"
export MAMBA_ROOT_PREFIX=/data/bnf/dev/jonas/hcv/micromamba
export VIRTITTA_MAMBA="$MAMBA_ROOT_PREFIX/micromamba"
export VIRTITTA_DEV_ENV="$MAMBA_ROOT_PREFIX/envs/virtitta"
export PYTHONPATH="$VIRTITTA_SRC"

while IFS= read -r run_dir; do
  test -d "$run_dir" || exit 1
  "$VIRTITTA_MAMBA" run -p "$VIRTITTA_DEV_ENV" python -m virtitta.cli import-run \
    --config "$VIRTITTA_DEV_CONFIG" --run-dir "$run_dir" || exit 1
done < "$bundle/run-dirs.txt"
```

## Development marker publication

Alternatively, after confirming that a development importer watches
`/fs1/results_dev/cron/virtitta` and uses the intended dev config, copy the bundle
to Hopper and publish only the corrected markers. This directory did not exist
when inspected during recovery. Creating this directory alone
does not configure a consumer. The production wrapper on Lennart uses production
Docker and must not be used for this dev import.

```bash
scp -r /tmp/virpipa-dev-recovery-20261111 rs-fe1:/tmp/
```

On Hopper:

```bash
bundle=/tmp/virpipa-dev-recovery-20261111
(cd "$bundle" && sha256sum -c SHA256SUMS) || exit 1
while IFS= read -r run_dir; do
  test -d "$run_dir" || exit 1
done < "$bundle/run-dirs.txt"
mkdir -p /fs1/results_dev/cron/virtitta
cp --no-clobber "$bundle"/dev-markers/*.sqlimport /fs1/results_dev/cron/virtitta/
```

These are handoff commands for operator execution; recovery preparation did not
publish markers or import into either database. Future runs should use the
[`HCV-dev` launcher entry](deployment.md#launcher-integration) with `--dev`.
