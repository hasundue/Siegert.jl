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
    zotero-mcp-src = {
      url = "github:54yyyu/zotero-mcp";
      flake = false;
    };
  };
  outputs =
    {
      self,
      nixpkgs,
      git-hooks-nix,
      treefmt-nix,
      zotero-mcp-src,
    }:
    let
      systems = [
        "x86_64-linux"
        "aarch64-linux"
        "x86_64-darwin"
        "aarch64-darwin"
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
            lib = nixpkgs.lib // {
              git-hooks-nix = git-hooks-nix.lib;
              treefmt-nix = treefmt-nix.lib;
            };
          }
        );
    in
    {
      devShells = forEachSystem (
        {
          pkgs,
          lib,
          system,
        }:
        let
          # Use Julia 1.11 for compatibility with withPackages
          julia = pkgs.julia_111-bin.withPackages [
            "JuliaFormatter"
            "TestEnv"
          ];
          treefmt = lib.treefmt-nix.mkWrapper pkgs {
            programs.nixfmt.enable = true;
            settings.formatter.julia = {
              command = "${lib.getExe julia}";
              options = [
                "-e"
                "using JuliaFormatter; JuliaFormatter.format(collect(ARGS))"
                "--"
              ];
              includes = [ "*.jl" ];
            };
          };
          git-hooks = lib.git-hooks-nix.${system}.run {
            src = ./.;
            hooks = {
              treefmt = {
                enable = true;
                package = treefmt;
              };
            };
          };
          zotero-mcp = pkgs.writeShellScriptBin "zotero-mcp" ''
            export LD_LIBRARY_PATH="${
              pkgs.lib.makeLibraryPath [ pkgs.stdenv.cc.cc.lib ]
            }''${LD_LIBRARY_PATH:+:}''$LD_LIBRARY_PATH"
            exec ${pkgs.uv}/bin/uvx \
              --python ${pkgs.python312}/bin/python3 \
              --from ${zotero-mcp-src} zotero-mcp "$@"
          '';
        in
        {
          default = pkgs.mkShell {
            packages = with pkgs; [
              julia
              treefmt
              zotero-mcp
            ];
            shellHook = git-hooks.shellHook;
          };
        }
      );
    };
}
