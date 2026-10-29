{
  description = "Experimental Julia implementation of Siegert pseudostates calculation (WIP)";
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
    git-hooks-nix = {
      url = "github:cachix/git-hooks.nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    treefmt-nix = {
      url = "github:numtide/treefmt-nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };
  outputs =
    {
      self,
      nixpkgs,
      git-hooks-nix,
      treefmt-nix,
    }:
    let
      systems = [
        "x86_64-linux"
        "x86_64-darwin"
      ];
      forEachSystem =
        f:
        nixpkgs.lib.genAttrs systems (
          system:
          f {
            inherit system;
            pkgs = import nixpkgs {
              inherit system;
              config.allowUnfree = false;
            };
          }
        );
    in
    {
      devShells = forEachSystem (
        {
          pkgs,
          system,
        }:
        let
          # Keep the formatter reproducible instead of relying on a
          # version-specific global environment (~/.julia/environments).
          juliaFormatter = pkgs.julia.withPackages [ "JuliaFormatter" ];
          treefmt = treefmt-nix.lib.mkWrapper pkgs {
            programs.nixfmt.enable = true;
            settings.formatter.julia = {
              command = "${juliaFormatter}/bin/julia";
              options = [
                "-e"
                "using JuliaFormatter; JuliaFormatter.format(collect(ARGS))"
                "--"
              ];
              includes = [ "*.jl" ];
            };
          };
          git-hooks = git-hooks-nix.lib.${system}.run {
            src = ./.;
            hooks = {
              treefmt = {
                enable = true;
                package = treefmt;
              };
            };
          };
        in
        {
          default = pkgs.mkShell {
            packages = with pkgs; [
              julia
              just
              nixfmt
              treefmt
            ];
            shellHook = ''
              # Keep project packages in the repository, but let the shared LSP
              # tooling resolve from the user depot. gitignore covers .julia/.
              julia_root="$(git rev-parse --show-toplevel 2>/dev/null || echo "$PWD")"
              export JULIA_DEPOT_PATH="$julia_root/.julia:$HOME/.julia"
              lsp_env="$julia_root/.julia/environments/nvim-lspconfig"
              if [ -e "$HOME/.julia/environments/nvim-lspconfig" ] && [ ! -e "$lsp_env" ]; then
                mkdir -p "$(dirname "$lsp_env")"
                ln -s "$HOME/.julia/environments/nvim-lspconfig" "$lsp_env"
              fi
              # Julia JLLs (libgfortran/libquadmath) are not found otherwise.
              export LD_LIBRARY_PATH="${
                pkgs.lib.makeLibraryPath [ pkgs.gfortran.cc.lib ]
              }''${LD_LIBRARY_PATH:+:}''$LD_LIBRARY_PATH"
              ${git-hooks.shellHook}
            '';
          };
        }
      );
    };
}
