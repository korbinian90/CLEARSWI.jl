# CLEARSWI.jl - issue list

Working notes from the architecture review of 2026-08-19..21. Uncommitted on
purpose. Cross-repository items and release ordering are in `issues-stack.md`
(X0..X7). Full evidence:
https://claude.ai/code/artifact/1dcb7a4a-4523-46fc-af23-7c5b0ecc3688

Measured on this machine unless marked *unverified*. State refers to branch
`claude/julia-repos-architecture-review-tj1zzz`.

The code is small and readable, about 800 lines. The problem is not the code, it
is the pressure on it: this is the package korec needs to extend, the package
the viewer reimplements, the package the ICE functor reimplements, and the
package with four unpushed local branches.

---

## Done on this branch

- **P2 `--qsm` cited the wrong algorithm.** It cited Kames et al. (RTS dipole
  inversion, from QSM.jl). The path actually goes through
  QuantitativeSusceptibilityMappingTGV, which is what the compiled app depends
  on, and the flag's own help text says "uses TGV QSM". Anyone following that
  citation would have described the wrong method in a paper. Now Langkammer 2015
  and Bredies 2014, plus ROMEO and, multi-echo only, MCPC-3D-S, which
  `qsm_romeo_B0` runs on the way there. With `--qsm-input` the user supplies a
  finished map, so nothing is claimed for it.
- **P3/P7** `saveconfiguration` writes arrays (sorted, for stable diffs) so the
  resolved echo times reach the file; the record moved to *after* echo-time
  resolution so it holds the values used, not the raw flags; and `writesteps`
  now also gets `citations_swi.txt`, covering only the methods that
  configuration uses. `settings_swi.txt` records the CLEARSWI,
  MriResearchTools, ROMEO and Julia versions.
- Version 1.6.2, requiring MriResearchTools `3.6`.

## Open

### C1. F6/X4 - `writesteps` names are treated as a contract but are not one
The highest-value item in this repository, because three downstream consumers
depend on it.

- In `laplacian_combine`, `filteredphase` is the 4D per-echo high-passed phase
  and `combinedphase` is the field combined *from* it.
- In `romeo_combine`, `combinedphase` is the unfiltered combination and
  `filteredphase` is the 3D filtered result.
- `laplacianslice_combine` saves the already-filtered array under the name
  `unwrappedphase` and never writes `filteredphase` at all.

`src/phase_processing.jl:88-127`. The viewer's `JULIA_STEM_ROLE` re-pointing
table exists solely to compensate for this. Writing the contract down, and
making the code obey it, retires most of that compensation, gives the ICE export
a schema to target, and unblocks the Rust parity work (X3, which is blocked on
this for the CLEAR-SWI half because the two `steps/` dumps disagree on the shape
of `phase_unwrapped`).

### C2. F3 - the four local branches
KorecSWI needs CLEARSWI 1.7 with a `Hooks` API. `master` is 1.6.1 with no
`Hooks`, no `getphasemask`, no `weighted_median`, and none of
`phase_unwrap_kernel`, `phase_mask`, `mag_sens_source`, `viewer_dir`. `origin`
has no branch carrying them. Four local branches implement the same two features
twice, and the largest measured ICE improvement (-65% ring artefact) has been
stranded unmerged since June. Gated on the public/private decision (X7.1), not
on engineering. The cost accrues while it waits.

### C3. F10 - abstract and untyped public API
`src/utility.jl:1-70`. `Data` holds `mag::AbstractArray`, an untyped header and
`TEs::AbstractVector`; `Options` has eleven fields, most untyped or `Union`-typed.
Nothing in the type system can reject a wrong option name or type.

### C4. F16 - CNR echo combination table is incomplete but accepts anything
`src/tissue.jl` has `factor` entries for 7T and a partial 3T, and accepts any
field symbol. The `T2s`, `T1`, `PD` and `factor` dictionaries are also non-const
globals (the F9 pattern), left alone here because this file sits inside the
CLEAR-SWI work under review.

### C5. F19/X2 - `eval` in the echo-time parser
`ext/ClearswiApp/argparse.jl:111,133` and `ext/ClearswiApp/caller.jl:100`.

### C6. F11 - the test data has forked
`test/data/small/{Mag,Phase}.nii` are byte-identical across MriResearchTools,
ROMEO, CompileMRI and mritools-binaries. **This repository's copies have
different checksums** under the same filenames. Same name, different data.
Reconcile when building the conformance dataset (X3).

### C7. `getTEs` disagrees with the other three implementations (F13/X1)
CLEARSWI returns `[1]` for single-echo where romeo returns the scalar `1`; tests
`isa Matrix` where romeo tests `isa AbstractMatrix`; drops romeo's `1 <` guard in
`1 < length(TEs) == neco`. Same flag, same job, different behaviour.
