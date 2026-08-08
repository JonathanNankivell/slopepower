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
          covr      # test coverage reports
          roxygen2  # generate NAMESPACE and man/ from #' comments -- 8.1.0, pinned above
          devtools  # load_all(), document(), check(), test()
          knitr     # vignette engine
          rmarkdown # render() and the html_vignette format; shells out to pandoc
        ];

        rEnv = pkgs.rWrapper.override { packages = rPkgs; };

        # pandoc is a system binary, not an R package: rmarkdown shells out to
        # whatever it finds on PATH.
        sysPkgs = [ pkgs.pandoc ];

        # nixpkgs' own `rstudio` package builds the whole Electron IDE from
        # source (electron_41, boost, its own R, ...) — slow, and frequently
        # behind upstream. Take the prebuilt .deb from Posit directly instead,
        # same as installing it on Ubuntu, and just repoint it at our R via
        # RSTUDIO_WHICH_R.
        rstudioVersion = "2026.07.1-147";
        rstudioSrc = pkgs.fetchurl {
          url = "https://download1.rstudio.org/electron/jammy/amd64/rstudio-${rstudioVersion}-amd64.deb";
          sha256 = "3a130a7209c9c9034c00440aa4b46164bbc5b75c1cf5588c98ef22a236ac1f4b";
        };

        # The Electron shell, its GTK/X11 chrome, and the rsession/rpostback
        # backend binaries between them need this stack (found by scanning
        # every ELF the .deb ships for unresolved NEEDED entries). Also
        # threaded onto LD_LIBRARY_PATH below: autoPatchelfHook only fixes up
        # RPATH for libraries a binary statically NEEDED-links (visible to
        # ldd) -- ANGLE/Chromium dlopen()s libEGL.so.1/libGLESv2.so.2 at
        # runtime, invisible to that scan, and silently falls through to a
        # SwiftShader software path that SIGILLs on at least one real CPU.
        rstudioRuntimeLibs = with pkgs; [
          alsa-lib
          at-spi2-atk
          at-spi2-core
          atk
          cairo
          cups
          dbus
          expat
          fontconfig
          freetype
          gdk-pixbuf
          glib
          gtk3
          libgbm
          libglvnd
          libsecret
          libuuid
          libxkbcommon
          mesa
          nspr
          nss
          openssl_3
          pango
          sqlite
          stdenv.cc.cc.lib
          systemdLibs
          xorg.libX11
          xorg.libXcomposite
          xorg.libXdamage
          xorg.libXext
          xorg.libXfixes
          xorg.libXrandr
          xorg.libXtst
          xorg.libxcb
          zlib
        ];

        rstudioFromDeb = pkgs.stdenv.mkDerivation {
          pname = "rstudio";
          version = rstudioVersion;
          src = rstudioSrc;

          nativeBuildInputs = with pkgs; [
            dpkg
            autoPatchelfHook
            makeWrapper
            copyDesktopItems
          ];

          buildInputs = rstudioRuntimeLibs;

          # Two categories of unresolved NEEDED entries are expected and fine:
          #  - libR.so: rsession links it directly but never finds it in the
          #    store; RStudio resolves this itself at runtime by shelling out
          #    to RSTUDIO_WHICH_R and putting its R_HOME/lib on
          #    LD_LIBRARY_PATH before it execs rsession, exactly as it would
          #    against a system R on Debian.
          #  - the bundled GitHub Copilot native modules ship prebuilds for
          #    other platforms/arches (musl, other libjpeg/libei/pipewire
          #    ABIs) that never load on this system; Copilot is optional and
          #    not part of this project's R workflow.
          autoPatchelfIgnoreMissingDeps = [ "*" ];

          unpackPhase = ''
            runHook preUnpack
            dpkg-deb -x "$src" .
            runHook postUnpack
          '';

          installPhase = ''
            runHook preInstall

            mkdir -p $out
            cp -r usr/lib $out/lib
            cp -r usr/share $out/share
            rm -f $out/share/applications/*.desktop

            # chrome-sandbox needs to be setuid-root to do its job; the Nix
            # store can't hold that, so run the Electron sandbox disabled
            # instead (RStudio doesn't browse the open web, so this trades
            # away renderer-process isolation the same way running it
            # unpackaged straight out of the .deb without `dpkg -i` would).
            #
            # The rest works around this being a plain devShell, not NixOS,
            # so there's no /run/opengl-driver and no guarantee of a working
            # GPU/Vulkan stack -- verified by actually launching this under
            # Xvfb + gdb, not just guessed at:
            #  - --ozone-platform=x11: without it, Ozone auto-detects Wayland
            #    whenever $WAYLAND_DISPLAY is set (true on most desktops
            #    today) and reaches for a real DRM render node even inside
            #    Xvfb; when the host has none available it hard CHECK-fails
            #    (a deliberate `ud2` abort, not a hardware fault).
            #  - LIBGL_ALWAYS_SOFTWARE / LIBGL_DRIVERS_PATH /
            #    __EGL_VENDOR_LIBRARY_FILENAMES: point GL/EGL at Mesa's own
            #    llvmpipe software rasterizer. --disable-gpu alone doesn't
            #    stop Chromium's Viz compositor from needing *some* GL
            #    backend, and without a real one it falls back to its
            #    bundled SwiftShader, which reliably hit the same `ud2`
            #    abort on the CPU this was tested on.
            #  - TMPDIR=/tmp: RStudio's Electron shell binds a Unix-domain
            #    socket for its single-instance lock under $TMPDIR, and
            #    sockaddr_un.sun_path has a hard 108-byte kernel limit. Each
            #    nested `nix develop`/`nix-shell` stacks another
            #    nix-shell.XXXXXX segment onto the ambient $TMPDIR (observed
            #    208 characters deep after two nested shells), which blows
            #    past that limit and hits the exact same `ud2` abort. Pin it
            #    to something always short instead of inheriting whatever
            #    the invoking shell chain happened to set.
            #
            # RSTUDIO_WHICH_R only tells RStudio which R binary to run `R
            # RHOME` against to find R_HOME -- it does NOT get rsession the
            # declared package library. rsession dlopens libR.so directly
            # rather than re-executing rEnv's bin/R wrapper script, so the
            # R_LIBS_SITE that script would normally set (via
            # rWrapper/makeWrapper, see nixpkgs' r-modules/wrapper.nix) never
            # reaches it. Ask that same wrapper what it resolves R_LIBS_SITE
            # to and bake the answer in directly instead.
            r_libs_site="$(${rEnv}/bin/R --no-echo --vanilla -e 'cat(Sys.getenv("R_LIBS_SITE"))')"
            makeWrapper $out/lib/rstudio/rstudio $out/bin/rstudio \
              --set RSTUDIO_WHICH_R ${rEnv}/bin/R \
              --set R_LIBS_SITE "$r_libs_site" \
              --set TMPDIR /tmp \
              --set LIBGL_ALWAYS_SOFTWARE 1 \
              --set LIBGL_DRIVERS_PATH ${pkgs.mesa}/lib/dri \
              --set __EGL_VENDOR_LIBRARY_FILENAMES ${pkgs.mesa}/share/glvnd/egl_vendor.d/50_mesa.json \
              --add-flags "--no-sandbox --disable-gpu --ozone-platform=x11" \
              --prefix PATH : ${pkgs.lib.makeBinPath sysPkgs} \
              --prefix LD_LIBRARY_PATH : ${pkgs.lib.makeLibraryPath rstudioRuntimeLibs}

            runHook postInstall
          '';

          desktopItems = [
            (pkgs.makeDesktopItem {
              name = "rstudio";
              exec = "rstudio %F";
              icon = "rstudio";
              desktopName = "RStudio";
              genericName = "RStudio";
              categories = [ "Development" "IDE" ];
            })
          ];

          meta = {
            description = "RStudio Desktop IDE, from the upstream .deb";
            homepage = "https://posit.co/products/open-source/rstudio/";
            license = pkgs.lib.licenses.agpl3Only;
            platforms = [ "x86_64-linux" ];
            mainProgram = "rstudio";
          };
        };
      in
      {
        devShells.default = pkgs.mkShell {
          packages = [ rEnv ] ++ sysPkgs;

          shellHook = ''
            echo "R $(R --version | head -n1 | cut -d' ' -f3) — tidyverse, nlme, haven, testthat, covr, roxygen2, devtools, knitr, rmarkdown"
            echo "pandoc $(pandoc --version | head -n1 | cut -d' ' -f2)"
          '';
        };

        # nix develop .#rstudio
        devShells.rstudio = pkgs.mkShell {
          packages = [ rEnv rstudioFromDeb ] ++ sysPkgs;

          shellHook = ''
            export RSTUDIO_WHICH_R=${rEnv}/bin/R
            echo "rstudio $(cat ${rstudioFromDeb}/lib/rstudio/version 2>/dev/null || echo unknown) — run 'rstudio' to launch"
          '';
        };

        packages.default = rEnv;
        packages.rstudio = rstudioFromDeb;
      });
}
