# Just recipes for developer helpers
#
# Usage:
#   - Ensure you're in the Nix dev shell: `nix develop`
#   - Run `just install-lsp` to (re)install the shared LanguageServer in
#     `~/.julia/environments/nvim-lspconfig`. The devShell symlinks it into the
#     project depot so `julials` finds it without a per-project install.

install-lsp:
	#!/usr/bin/env bash
	set -euo pipefail
	case "$(julia --version)" in
		"julia version 1.13."*) ;;
		*)
			echo "error: enter the dev shell with \`nix develop\` first" >&2
			exit 1
			;;
	esac
	JULIA_DEPOT_PATH="${HOME}/.julia" julia --project=@nvim-lspconfig -e 'using Pkg; Pkg.add("LanguageServer"); Pkg.precompile()'

test *targets:
	#!/usr/bin/env bash
	set -- {{targets}}
	if [ $# -eq 0 ]; then
		julia --project -e 'using Pkg; Pkg.test()'
	else
		args=$(printf '"%s",' "$@" | sed 's/,$//')
		julia --project -e "using Pkg; Pkg.test(test_args=[$args])"
	fi
