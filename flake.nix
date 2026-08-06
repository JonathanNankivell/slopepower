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

        # roxygen2 8.1.0 (CRAN, 2026-07-31) isn't packaged in nixpkgs yet --
        # the pinned nixpkgs-unstable still ships 8.0.0, and DESCRIPTION
        # requires 8.1.0 (see the "Regenerate NAMESPACE with roxygen2 8.1.0"
        # commit). Pin the CRAN release by hand rather than waiting on the
        # next R-package regeneration to land upstream; every other package
        # in this shell stays exactly what nixpkgs provides.
        #
        # 8.1.0 depends on rdtools, a brand-new CRAN package that also isn't
        # in nixpkgs yet, so it has to be built the same way.
        rdtools = pkgs.rPackages.buildRPackage {
          name = "r-rdtools-0.1.0";
          src = pkgs.fetchurl {
            url = "https://cran.r-project.org/src/contrib/rdtools_0.1.0.tar.gz";
            sha256 = "0jq9dk3dj4glrdn1n6nmgi1w4a20ig303d7lfv525lwaqzqaxzch";
          };
        };

        roxygen2 = pkgs.rPackages.roxygen2.overrideAttrs (old: {
          name = "r-roxygen2-8.1.0";
          version = "8.1.0";
          src = pkgs.fetchurl {
            url = "https://cran.r-project.org/src/contrib/roxygen2_8.1.0.tar.gz";
            sha256 = "14dn3ya51z9ipivsl3larmx6lgph997rwhcwmfqdpkx8bdamvqvf";
          };
          propagatedBuildInputs = old.propagatedBuildInputs ++ [ rdtools ];
          nativeBuildInputs = old.nativeBuildInputs ++ [ rdtools ];
        });

        # Add any extra CRAN packages here (names match pkgs.rPackages.*).
        rPkgs = with pkgs.rPackages; [
          tidyverse
          nlme      # lme()/gls(), correlation & variance structures
          haven     # read Stata .dta files
          testthat  # unit tests
          roxygen2  # generate NAMESPACE and man/ from #' comments -- 8.1.0, pinned above
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
