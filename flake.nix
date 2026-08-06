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
          nlme      # lme()/gls(), correlation & variance structures
          haven     # read Stata .dta files
          testthat  # unit tests
          roxygen2  # generate NAMESPACE and man/ from #' comments
          devtools  # load_all(), document(), check(), test()
          knitr     # vignette engine
          rmarkdown # render() and the html_vignette format; shells out to pandoc
        ];

        rEnv = pkgs.rWrapper.override { packages = rPkgs; };
        rstudioEnv = pkgs.rstudioWrapper.override { packages = rPkgs; };

        # pandoc is a system binary, not an R package: rmarkdown shells out to
        # whatever it finds on PATH.
        sysPkgs = [ pkgs.pandoc ];
      in
      {
        devShells.default = pkgs.mkShell {
          packages = [ rEnv ] ++ sysPkgs;

          shellHook = ''
            echo "R $(R --version | head -n1 | cut -d' ' -f3) — tidyverse, nlme, haven, testthat, roxygen2, devtools, knitr, rmarkdown"
            echo "pandoc $(pandoc --version | head -n1 | cut -d' ' -f2)"
          '';
        };

        # nix develop .#rstudio
        devShells.rstudio = pkgs.mkShell {
          packages = [ rstudioEnv ] ++ sysPkgs;
        };

        packages.default = rEnv;
      });
}
