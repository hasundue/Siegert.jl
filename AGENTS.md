# AGENTS Guide (Julia + Nix)

## Overview
- Julia + Nix project; keep changes minimal and focused; prefer small, clear PRs.

## Environment
- Julia 1.12.x; enter dev shell with `nix develop`.

## Setup
- `julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'`.
- Build rarely: `julia --project -e 'using Pkg; Pkg.build()'`.

## Workflow
- Before pushing: `treefmt -c`, then `julia --project -e 'using Pkg; Pkg.test()'`.

## Testing
- All tests: `julia --project -e 'using Pkg; Pkg.test()'`.
- Single file: `julia --project test/<file>.jl`.
- Name filter: `using Test; Test.@testset filter = t -> occursin("NAME", string(t)); include("test/runtests.jl")`.

## Code Style
- 2-space indent; ≲100 cols; trailing commas on multiline; one export per line; triple-quoted docstrings.

## Imports
- Group stdlib, external, local; one per line; prefer `using`; qualify on conflicts; alphabetize groups.

## Types & Naming
- Generic code; annotate struct fields; concrete element types in containers; avoid globals; use multiple dispatch.
- CamelCase (modules/types), snake_case (functions/vars), `is_`/`has_` predicates, UPPER_CASE constants.

## Errors
- Validate inputs; throw `ArgumentError`/`DomainError`/`MethodError`; use `@assert` for internal invariants only.

## Docs
- Every public API: signature, brief, example, edge cases/units.

## Nix & Formatting
- Keep `flake.nix` formatted; `treefmt` uses `nixfmt`; don’t bypass hooks.

## Research (Zotero MCP)
- Read-only. If `http://127.0.0.1:23119` is reachable (Zotero 7 Local API), use it to search/read. Assume project papers are in the `Siegert.jl` collection to constrain lookups. All edits happen in Zotero.
