# C++ optimization workflow

This is the working policy for further C++ optimization, recorded after issue #470 and PR #471. The owner requested that the lessons from the 6 h 36 min implementation remain part of future work. Keep numerical correctness and required integration checks; reduce late discoveries, repeated qualification, and coordination overhead.

## Start with one measurable hypothesis

1. Profile the latest qualified implementation on a representative user workload. Reuse frozen fixtures and record binary, source, provider, thread count, memory policy, and sampling window. Separate integration, output, preparation, and process lifetime. Budget the campaign from startup and comparison costs too, especially for large persisted matrix tables. Check baseline repeatability before imposing an exact whole-model numerical gate. Sampling percentages are estimates of CPU attribution, not measured end-to-end speedups.
2. Read the existing spectral-kernel-benchmarks history, closed-issue conclusions, and linked adoption decisions relevant to the proposed algorithm. Reuse the recorded investigation; do not repeat rejected experiments without identifying a changed assumption.
3. Choose one bottleneck and estimate the maximum possible whole-model benefit from its measured share. State the intended change, affected families, correctness checks, memory constraints, and go/no-go criterion before implementation. Do not pursue tiny-grid speedups unless the owner changes that preference.
4. Use a 30–45 minute investigation checkpoint as a default planning target, not a correctness deadline. At that point report evidence and decide whether to continue, narrow the experiment, or reject it. Do not turn an unsuccessful local optimization into an unrequested architecture rewrite.

## Keep the team small and responsibilities explicit

- When delegation is authorized, default to one coordinator, one implementation worker, and one independent reviewer/test worker. Add workers only for concrete independent work that shortens the critical path. Do not recursively expand the team by default.
- The coordinator owns shared interfaces, integration, verification scheduling, and the final decision. Freeze the relevant interface before parallel adapters are written; give each writer explicit file ownership or an isolated worktree.
- Use cheaper models and moderate reasoning for mechanical adapters, repository scans, bounded tests, and log triage. Use stronger reasoning for numerical algorithms, ownership/lifetime design, and subtle correctness review. Escalate based on difficulty, not by default.
- Reuse workers within an increment. Give new workers a compact task brief with paths, interfaces, exclusions, and acceptance checks rather than the entire conversation. Return a diff summary, exact tests/results, and unresolved findings.
- Only one owner schedules shared builds, MATLAB runs, and final benchmarks. Timing qualification requires an idle host; do not compete with the experiment or other workers. Do not alter a running experiment to obtain a benchmark window.
- Check free disk space before a campaign. Budget for retained representative outputs, both active pair outputs, and failure headroom; process memory is a separate budget. After successful comparison, redundant generated payloads may be removed with their hashes and comparison records preserved. If storage interrupts a campaign, retain completed evidence and resume only the incomplete workload after rechecking frozen identities.

## Detect cheap failures before expensive qualification

1. Run repository/provenance checks and the affected warning-clean compiler builds early. Check supported GCC as well as local Clang before declaring a C++ candidate ready; use actual Linux CI where local SDK constraints prevent a complete GCC build. When adding a private compiled helper, search for sibling source paths in receipt and export lists and add its fingerprint before freezing.
2. Run the affected native tests and targeted MATLAB comparisons while editing. Include consumer combinations implicated by the change: components, density references, forcing order, sampling, retries, and storage expectations where relevant. Audit existing assertions when ownership/accounting changes.
3. Finish independent correctness review and resolve findings before freezing the final runtime candidate. Do not call an intermediate candidate final while a review may still change runtime code.
4. Run the appropriate combined native suite and affected MATLAB comparisons, then freeze source, executable, fixtures, and benchmark protocol. Use short exploratory timings during development; run the full required paired campaign for the reviewed candidate.
5. Refresh source-linked receipts coherently and submit the combined candidate to required hosted integration checks. Aim for one successful final campaign and hosted submission; genuine later failures still require correction and relevant requalification. Never relax numerical tolerances or bypass required checks to meet this target.

## Reuse valid evidence and keep the record compact

- Keep one verification ledger: check, source/input identity, result, and the change that would invalidate it. Before repeating a passing gate, identify the intervening relevant change.
- Preserve exact commit provenance. Distinguish changes to compiled inputs, toolchains, dependencies, fixtures, and scientific references from evidence-only or unrelated test changes. Reuse earlier evidence only with an explicit equivalence rationale; do not redesign provenance validators inside an unrelated optimization task.
- Automate repeated preparation/receipt steps when they materially reduce manual work. Prefer one entry point over successive ad hoc commands. Do not create a broad tooling project before measuring the next bottleneck.
- Preserve failed campaigns, numerical comparisons, hashes, and compact summaries. Store large raw profiles and payloads in the benchmark artifact archive, linking their hashes from the repository. Avoid embedding repeated full logs or duplicate large payloads in the main conversation.
- Report implementation time, verification/rework time, benchmark time, and waiting separately when available. Agent count is not a throughput measure. State uncertainty rather than presenting an inferred time allocation as measured.
- For documentation-only delivery, follow the shared policy: enable auto-merge behind required checks and return without actively polling unless the user asks to wait or a failure needs diagnosis.

## Why these rules exist

Issue #470 ran approximately 13:07–19:42 local time on 2026-09-11: 3 h 20 min through initial integrated correctness, about 1 h 03 min through review and initial passing performance qualification, and about 2 h 13 min of further qualification, corrections, and merge. These are chronological intervals, not an exact active-time breakdown. Four full frozen campaigns consumed about 30 minutes of measurement, with additional preparation and receipt work. Five hosted CI attempts preceded the successful approximately 22-minute final run.

Late findings included GCC portability problems, missing density-source lifetime handling, component-energy cache identity, stale source-selection hashes, and obsolete memory assertions. The scientific bugs needed fixing; the workflow failure was discovering them after expensive qualification. Details and preserved failures are in [the #470 verification ledger](../.github/planning/issue-470-verification.md). [PR #471](https://github.com/JeffreyEarly/wave-vortex-model/pull/471) is merged and [issue #470](https://github.com/JeffreyEarly/wave-vortex-model/issues/470) is closed; historical pending statements in the ledger describe their recorded stage.
