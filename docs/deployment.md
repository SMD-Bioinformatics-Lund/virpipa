# Hopper deployment

Production VirPipa releases are installed from the Hopper bare repository, not
from the development checkout:

- Git source: `/fs1/jonas/bare/virpipa/`
- development and test checkout: `/fs1/jonas/src/virpipa/`
- production root: `/fs1/pipelines/virpipa/`

## Deploy a revision

Load the Hopper Nextflow environment or let the deploy utility load the standard
production modules, then name the branch, tag, or commit explicitly:

```bash
bash scripts/deploy_hopper.sh --dry-run v1.1
bash scripts/deploy_hopper.sh v1.1
```

The utility extracts the complete tracked tree into
`/fs1/pipelines/virpipa/releases/<timestamp>-<git-describe>/`, validates the production
Nextflow configuration, writes `VERSION` and `SHA256SUMS`, and atomically points
`current` at the new release. Characters unsuitable for a directory name are
replaced with `_`; `VERSION` retains the exact description and full commit.
Deployments are locked so that only one can run at a time. Existing releases are retained.

Before promotion, the current release is checked for failed checksums,
unexpected files, and files newer than its deployment metadata. Deployment
stops if direct changes are found. After reviewing and preserving any hotfix,
`--force` permits creation of a new release without modifying the old one.

Inspect the active version with:

```bash
readlink -f /fs1/pipelines/virpipa/current
cat /fs1/pipelines/virpipa/current/VERSION
```

## Launcher integration

The live `[HCV]` entry in `/fs2/sw/bnf-scripts/pipeline_files.config` should use
`/fs1/pipelines/virpipa/current`. `HCV-dev` should continue to use the development
checkout.

Until the launcher canonicalizes the pipeline path, do not deploy while an HCV
coordinator is queued or running. The symlink works, but a queued coordinator
uses the release selected when it starts and a switch during an analysis could
make later `${projectDir}` references cross release boundaries.

The preferred launcher enhancement is to split the configured pipeline value
into its first path token and remaining arguments, resolve only the path with
Perl's `abs_path()`, and write the resolved release path into the generated
`sbatch` file. Verify this using `start_nextflow_analysis.pl --nostart` before
enabling it. The deploy utility deliberately does not modify the shared launcher.

## Manual rollback

First verify the selected release using its checksums. Then take the deployment
lock and atomically repoint `current`:

```bash
release=/fs1/pipelines/virpipa/releases/<release>
(cd "$release" && sha256sum -c SHA256SUMS)
(
  flock -n 9 || exit 1
  ln -s "releases/$(basename "$release")" /fs1/pipelines/virpipa/.current.rollback
  mv -Tf /fs1/pipelines/virpipa/.current.rollback /fs1/pipelines/virpipa/current
) 9>/fs1/pipelines/virpipa/.deploy.lock
```

## Container caches

Production processes execute prebuilt SIF images directly from
`/fs1/resources/containers/`; these are shared runtime images, not a Singularity
download cache. Singularity's normal per-user cache may remain under
`$HOME/.singularity/cache`. Do not place runtime caches inside a release.
