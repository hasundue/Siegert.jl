# AGENTS Guide (Julia + Nix)

## Overview
- Julia + Nix project; keep changes minimal and focused; prefer small, clear PRs.

## Environment
- Julia 1.11.x; enter dev shell with `nix develop`.

## Setup
- `julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'`.
- Build rarely: `julia --project -e 'using Pkg; Pkg.build()'`.

### Pkg dependencies (normal + extras)
- Never edit `Project.toml` or `Manifest.toml` by hand. Use `Pkg` for all dependency changes (including stdlibs like `LinearAlgebra`), and let it write files. Example: `julia --project -e 'using Pkg; Pkg.add("LinearAlgebra"); Pkg.resolve()'`.
- Normal deps: `Pkg.add("PkgName")` / `Pkg.rm("PkgName")`.
- Test-only deps: `Pkg.add("PkgName"; target=:extras)` / `Pkg.rm("PkgName"; target=:extras)`. Also list them under `[targets]` as `test = ["Test", "PkgName", ...]` (keep "Test" explicit here).
- Maintain `[compat]` manually when bumping versions.

## Workflow
- Commit messages (pre-1.0): use `<scope>: body` style. Example: `dvr: return ψ, z, Λ from jacobi_dvr_basis; update tests/examples`.
- Before pushing: `treefmt -c`, then `julia --project -e 'using Pkg; Pkg.test()'`.

## Testing
- All tests: `julia --project -e 'using Pkg; Pkg.test()'`.
- Single file: `julia --project test/<file>.jl`.
- Name filter: `using Test; Test.@testset filter = t -> occursin("NAME", string(t)); include("test/runtests.jl")`.
- Prefer targeted runs while iterating. Example: `julia --project test/waves_tests.jl` instead of the full suite. Only run `Pkg.test()` before pushing.

## Code Style
- 4-space indent; ≲100 cols; trailing commas on multiline; one export per line; triple-quoted docstrings.

## Imports
- Group stdlib, external, local; one per line; prefer `using`; qualify on conflicts; alphabetize groups.

## Types & Naming
- Generic code; annotate struct fields; concrete element types in containers; avoid globals; use multiple dispatch.
- CamelCase (modules/types), snake_case (functions/vars), `is_`/`has_` predicates, UPPER_CASE constants.

## Errors
- Validate inputs; throw `ArgumentError`/`DomainError`/`MethodError`; use `@assert` for internal invariants only.

## Docs
- Every public API: signature, brief, example, edge cases/units.
- Physics implementations: always reference the original paper in code headers/docstrings with full citation (authors, title, venue, year, DOI); cite equation numbers used and call out any scaling/convention choices.
- When porting formulas, state assumptions and units; link follow-up references when extending (e.g., multi-channel, l>0).

## Nix & Formatting
- Keep `flake.nix` formatted; `treefmt` uses `nixfmt`; don’t bypass hooks.

## Research (Zotero MCP)
- Read-only. If `http://127.0.0.1:23119` is reachable (Zotero 7 Local API), use it to search/read. Assume project papers are in the `Siegert.jl` collection to constrain lookups. All edits happen in Zotero.
