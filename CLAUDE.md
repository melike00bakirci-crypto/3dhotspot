# 3Dhotspot — standing operating procedure

Instructions for Claude Code sessions acting as Lead in this repo, in particular for
`hotspot3d run --gene GENE` executions. This applies to every future gene run unless the user
explicitly says otherwise for that run.

## Gene run defaults

- Never stop to ask permission mid-run. Proceed autonomously.
- Never reduce N_CAP, B, max_wall_seconds, or any other budget to force a run to fit. If a
  stage doesn't fit its budget, run what fits and report the achieved coverage honestly.
- Always report the terminal state as it is: COMPLETED, COMPLETED_NEGATIVE, UNDERPOWERED, or
  BLOCKED. Never adjust methodology, thresholds, or cohort definitions to influence the outcome.
- On stage failure: route through contract validation → pipeline-qa-debugger. If the bounded
  repair limit is reached, mark that stage BLOCKED and stop cleanly — never retry indefinitely.
- Default worker count: 16, worker_nice: 10 (per the shared-node policy already in place).
  Adjust only if the user specifies otherwise for a concurrent run.
- When running multiple genes in one session: never edit anything under `src/`, `config/`,
  `docs/`, or `.claude/` mid-batch. If a gene fails, record the diagnosis, mark it FAILED, and
  move to the next gene — never halt the batch.
- At the end of any run or batch, write one summary report: per-stage wall time, terminal
  state, all warnings, anything requiring QA intervention.
