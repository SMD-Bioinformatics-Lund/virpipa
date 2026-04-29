# AGENTS.md - VirPipa Repository Guidance

## General coding guidelines

### 1. Think Before Coding

**Don't assume. Don't hide confusion. Surface tradeoffs.**

Before implementing:
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them - don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

### 2. Simplicity First

**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

### 3. Surgical Changes

**Touch only what you must. Clean up only your own mess.**

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it - don't delete it.

When your changes create orphans:
- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: Every changed line should trace directly to the user's request.

### 4. Goal-Driven Execution

**Define success criteria. Loop until verified.**

Transform tasks into verifiable goals:
- "Add validation" → "Write tests for invalid inputs, then make them pass"
- "Fix the bug" → "Write a test that reproduces it, then make it pass"
- "Refactor X" → "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:
```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.

## Purpose

VirPipa is an HCV probe-capture assembly pipeline built with Nextflow DSL2. The standing engineering goal is output parity with the original bash pipeline in `scripts/hcvpipe.sh`, while using Nextflow-native structure where it improves the implementation without changing emitted results.

## Always-On Rules

- Treat `master` as the active base branch unless the user explicitly asks for another branch or a new feature branch.
- Never add container images to git.
- The repo has two execution environments: a local laptop for light tests and Hopper for real runs.
- On Hopper, use `runme.sh` / `runme_test.sh` rather than ad hoc `nextflow run` commands because the environment is site-specific.
- If you commit changes that must run on Hopper, push to the Hopper bare repo before using the wrapper scripts.
- Prefer bash in Nextflow `script:` blocks unless Groovy is clearly the better tool.
- Do not change pipeline behavior speculatively. Keep changes surgical and tie each change directly to the requested behavior.

## Development Policy

- Build the pipeline one step at a time and verify each new or changed module against the known bash behavior.
- Check module parameters and emitted filenames before coding.
- Match the bash pipeline outputs, especially in the published `results/` tree. Nextflow internals may differ; published outputs should not.
- Checksums are not the only parity signal. Some formats embed paths, timestamps, or command lines, so compare payloads when needed.
- Prefer metadata-preserving tuple channels over bare paths.
- When a workflow decision depends on `sample_id` vs `run_name`, be explicit. Batch behavior has already broken on that distinction before.

## Environment Summary

### Laptop

- Use the `skrotis` conda/mamba environment for local work.
- The laptop is suitable for tiny tests and fixture-backed module checks, not full real-sample runs.
- The mounted HPC filesystem is available under `/mnt/fs1`.
- Repo-local containers live under `assets/containers/` when using local container-backed tests.

### Hopper

- Hopper has no internet.
- Hopper paths use `/fs1/...`.
- Site wrappers and environment modules matter; use the provided run scripts.
- The Hopper bare repo remote is the publication point for code that the wrappers pull from.

## Canonical Commands

### Config / validation

```bash
nextflow config -profile hpc,apptainer
nextflow config -profile slurm,hpc,apptainer
```

### Common local module checks

```bash
nextflow run test_module.nf -profile local,tiny --module subsample --subsample_reads 25
nextflow run test_module.nf -profile local_containers,tiny --module hostile
nextflow run test_module.nf -profile local --module bestref
nextflow run test_module.nf -profile local_containers --module subtype
nextflow run test_module.nf -profile local_containers --module vadr
nextflow run test_module.nf -profile local_containers --module resistance
```

### Hopper wrappers

```bash
bash runme_test.sh <module_name>
ssh rs-fe1 'cd /fs1/jonas/src/virpipa && bash -l runme.sh resume'
ssh rs-fe1 'cd /fs1/jonas/src/virpipa && bash -l runme.sh'
```

## Reference Docs

Use these docs instead of bloating this file with procedural detail:

- `docs/pipeline_workflow.md`: current workflow structure and branch logic.
- `docs/agent_environment.md`: detailed laptop/Hopper environment facts and operational commands.
- `docs/module_test_reference.md`: fixture inventory, module test commands, and parity caveats.
- `docs/output.md`: published output contract.
- `docs/resistance_pipeline.md`: resistance annotation details.
- `docs/usage.md`: CLI parameters and launcher-facing aliases.

## Skills

The following workflow-specific guidance is better handled by skills than by this always-on file:

- VirPipa parity verification workflow.
- Hopper run failure triage.

If those skills are installed in Codex, use them when the task specifically involves parity validation or Hopper debugging.
