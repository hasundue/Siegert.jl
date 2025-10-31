# Just recipes for developer helpers
#
# Usage:
#   - Ensure you're in the Nix dev shell: `nix develop`
#   - Run: `just install-lsp`
#
# This installs Julia LSP tooling into the `nvim-lspconfig` environment
# the server, try: `julia --project=@nvim-lspconfig ...` instead.

install-lsp:
	julia --project=@nvim-lspconfig -e 'using Pkg; Pkg.add(["LanguageServer","SymbolServer","StaticLint","JuliaFormatter"]); Pkg.precompile()'
