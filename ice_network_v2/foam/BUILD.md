# Reproducible OpenFOAM build path (pinned: OpenFOAM ESI v1912, Ubuntu source/binary package 1912.200626-2build3)

Host used: Ubuntu 24.04.5 (noble), x86_64, gcc 13.3.0, 4 cores, 15 GB RAM, root in a container. Disk: sources 272 MB (+45 MB tarball), built user apps 17 MB, a full set of verification case directories ~1 GB (scratch, not in git).
The vendor repository script (`dl.openfoam.com/add-debian-repo.sh`) was NOT run (blocked by the sandbox permission classifier, not worked around); only Ubuntu-archive packages and the matching Ubuntu source package are used, so binaries (libraries) and headers/tools come from the SAME source version. No mixture of distributions.

```sh
# 1. binaries + libraries (stock solvers icoFoam, interFoam, scalarTransportFoam, simpleFoam, checkMesh, blockMesh ...)
apt-get update && apt-get install -y openfoam openfoam-examples           # 1912.200626-2build3 ; examples package is EMPTY of tutorials
# stock solvers need:  export WM_PROJECT_DIR=/usr/share/openfoam FOAM_ETC=/usr/share/openfoam/etc   and (to avoid an FPE abort seen in this Debian build inside PCG) unset FOAM_SIGFPE for interFoam runs
# 2. matching development sources (headers, wmake): enable deb-src for the Ubuntu archive, download the source package
cp /etc/apt/sources.list.d/ubuntu.sources /root/ubuntu.sources.bak && sed -i 's/^Types: deb$/Types: deb deb-src/' /etc/apt/sources.list.d/ubuntu.sources && apt-get update
mkdir -p /tmp/claude-0/ofsrc && cd /tmp/claude-0/ofsrc && apt-get source --download-only openfoam && dpkg-source -x openfoam_1912.200626-2build3.dsc ofdeb   # applies the Debian patches (gcc-13, hardening)
apt-get install -y flex libfl-dev libopenmpi-dev zlib1g-dev bison cmake
# 3. build the wmake tools and generate lnInclude symlinks (headers only; libraries are the installed ones in /usr/lib)
source ice_network_v2/foam/build_env.sh            # replicates debian/rules exports; FOAM_LIBBIN=/usr/lib
cd $SRC/wmake/src && make && cd $SRC/src && wmakeLnIncludeAll
# 4. build the project's applications (user binaries go to /tmp/claude-0/ofbuild/bin)
cd ice_network_v2/foam/testApp && wmake            # copy of scalarTransportFoam: build test
cd ../interThermalFoam && wmake                    # interFoam + sensible-heat equation in both phases
cd ../stefan1DFoam && wmake                        # ALE sharp-interface Stefan solver
```
Minimal execution checks (recorded): the custom-built `scalarTransportFoamBuildTest` reproduces the stock `scalarTransportFoam` output bit-for-bit on a Graetz case (`cmp` of the T field); `interThermalFoam` was validated by the cases in `freesurface/`; `stefan1DFoam` by `stefan1DFoam/stefan_of.py`.
Known build/runtime issues (this Debian build): `postProcess -func writeCellCentres` aborts (sha1 IOstream error); floating-point trapping aborts `interFoam` in `PCG::scalarSolve` early (FOAM_SIGFPE unset in `freesurface/common.py`; NaN checks in the analysis scripts); wmake emits `-Walloc-size-larger-than` warnings from `List.C` (harmless).
