{
  description = "R environment for transpiling slopepower (Stata -> R)";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };

        # Add any extra CRAN packages here (names match pkgs.rPackages.*).
        rPkgs = with pkgs.rPackages; [
          tidyverse
          lme4      # provides lmer()
          lmerTest  # p-values / anova for lmer models
          nlme      # lme()/gls(), correlation & variance structures
          haven     # read Stata .dta files
          testthat  # unit tests
          roxygen2  # generate NAMESPACE and man/ from #' comments
        ];

        rEnv = pkgs.rWrapper.override { packages = rPkgs; };
        rstudioEnv = pkgs.rstudioWrapper.override { packages = rPkgs; };
      in
      {
        devShells.default = pkgs.mkShell {
          packages = [ rEnv ];

          shellHook = ''
            echo "R $(R --version | head -n1 | cut -d' ' -f3) — tidyverse, lme4, lmerTest, nlme, haven, testthat, roxygen2"
          '';
        };

        # nix develop .#rstudio
        devShells.rstudio = pkgs.mkShell {
          packages = [ rstudioEnv ];
        };

        packages.default = rEnv;
      });
}
