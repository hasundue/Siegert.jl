# Just recipes for developer helpers
#
# Usage:
#   - Ensure you're in the Nix dev shell: `nix develop`
#   - Run `just install-lsp` once to populate the Julia LSP environment in the
#     project-local depot (`.julia`).

install-lsp:
	#!/usr/bin/env bash
	set -euo pipefail
	if [ -z "${JULIA_DEPOT_PATH:-}" ]; then
		echo "error: enter the dev shell with \`nix develop\` first" >&2
		exit 1
	fi
	julia --project=@nvim-lspconfig -e 'using Pkg; Pkg.add("LanguageServer"); Pkg.precompile()'

test *targets:
	#!/usr/bin/env bash
	set -- {{targets}}
	if [ $# -eq 0 ]; then
		julia --project -e 'using Pkg; Pkg.test()'
	else
		args=$(printf '"%s",' "$@" | sed 's/,$//')
		julia --project -e "using Pkg; Pkg.test(test_args=[$args])"
	fi
