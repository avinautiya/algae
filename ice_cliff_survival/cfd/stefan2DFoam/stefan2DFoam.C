/*---------------------------------------------------------------------------*\
  stefan2DFoam (ice_cliff_survival m1): 2-D heat conduction in a liquid layer bounded by a MOVING ice interface (ALE, columns of the structured mesh move in x), sensible-plus-latent bookkeeping, optional static debris shielding.
  Geometry: liquid 0 < x < s(y,t), 0 < y < Ly (one cell deep, depth dz). Hot wall x = 0 (T_w(y) fixed, supplies liquid when the density jump requires it); sides y = 0, Ly adiabatic; ice (temperate, T_m = 0, at rest) for x > s(y,t).
  Equation (ALE, conservative):  d/dt int T dV + int T (u_w - w).dS - int alpha grad T.dS = 0.  Equal-density limit u_w = 0; otherwise u_w = uniform (1 - rho_i/rho_w) <sdot> x^ with <sdot> = swept-volume rate / hot-wall area (global mass balance exact; local lateral flow NOT represented).
  Interface (n from ice into liquid): rho_i L m_n = q, q = heat flux arriving at the ice per area = the SAME discrete flux the temperature equation removes through the face (k * face snGrad, with debris shielding via a mixed BC: q = h_eff T_face).
  Kinematics: the interface is a graph x = s(y,t); a face's swept volume is V_f = m_n A_f dt; node displacement dx_j = mean of the adjacent face displacements dx_f = V_f/(dy dz) (so the TOTAL swept volume is exactly sum V_f; per-face smoothing error converges with dy).
  Heun predictor-corrector in time. Debris shielding is STATIC: h_eff(y) = k_deb/d_deb per interface row (large = bare); debris heat capacity ignored (quasi-steady layer).
  Ledgers written every step: sensible storage, hot-boundary conductive+advective supply, interface discrete conductive loss, latent rho_i L V_swept (V_swept from the geometry), GCL check against the OpenFOAM mesh flux.
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "OFstream.H"
#include "mixedFvPatchFields.H"

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary prop(IOobject("stefanProperties", runTime.constant(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE));
    const scalar rhoW = readScalar(prop.lookup("rhoW")), cpW = readScalar(prop.lookup("cpW")), kW = readScalar(prop.lookup("kW"));
    const scalar rhoI = readScalar(prop.lookup("rhoI")), Lf = readScalar(prop.lookup("L"));
    const scalar s0 = readScalar(prop.lookup("s0")), t0 = readScalar(prop.lookup("t0")), Ly = readScalar(prop.lookup("Ly")), dz = readScalar(prop.lookup("dz"));
    const label ny = readLabel(prop.lookup("ny"));
    const bool densityJump = readBool(prop.lookup("densityJump"));
    const scalarField tRow(prop.lookup("tRow")), TwRow(prop.lookup("TwRow")), lamRow(prop.lookup("lambdaRow")), bRow(prop.lookup("bRow")), hShield(prop.lookup("hShield"));
    if (tRow.size() != ny || TwRow.size() != ny || lamRow.size() != ny || bRow.size() != ny || hShield.size() != ny) FatalErrorInFunction << "tRow, TwRow, lambdaRow, bRow, hShield must have ny entries" << exit(FatalError);
    const scalar alphaW = kW/(rhoW*cpW), r = rhoI/rhoW, dy = Ly/ny;

    const label icePatch = mesh.boundaryMesh().findPatchID("ice"), hotPatch = mesh.boundaryMesh().findPatchID("hot");
    if (icePatch < 0 || hotPatch < 0) FatalErrorInFunction << "patches 'ice' and 'hot' required" << exit(FatalError);

    // reference coordinates: xi = x/s0, row index of every point (by y)
    pointField xi0(mesh.points()); labelList jRow(xi0.size());
    forAll(xi0, pi) { jRow[pi] = label(xi0[pi].y()/dy + 0.5); xi0[pi].x() /= s0; }
    scalarField sNode(ny + 1, s0);
    auto movePts = [&](const scalarField& sn) { pointField pp(xi0); forAll(pp, pi) { pp[pi].x() *= sn[jRow[pi]]; } mesh.movePoints(pp); };

    volScalarField T(IOobject("T", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    auto Iint = [&](scalar eta, scalar b) { return 0.5*std::sqrt(constant::mathematical::pi)*(::erf(eta - b) + ::erf(b)); };
    forAll(T, ci)
    {
        const label j = min(ny - 1, label(mesh.C()[ci].y()/dy)); const scalar eta = mesh.C()[ci].x()/(2*std::sqrt(alphaW*tRow[j])); const scalar Cc = TwRow[j]/Iint(lamRow[j], bRow[j]);
        T[ci] = TwRow[j] - Cc*Iint(eta, bRow[j]);
    }
    // hot wall values per row
    {
        fvPatchScalarField& Th = T.boundaryFieldRef()[hotPatch]; forAll(Th, fi) { const label j = min(ny - 1, label(mesh.Cf().boundaryField()[hotPatch][fi].y()/dy)); Th[fi] = TwRow[j]; }
    }
    mixedFvPatchScalarField& Tice = refCast<mixedFvPatchScalarField>(T.boundaryFieldRef()[icePatch]);
    labelList jFace(Tice.size()); forAll(jFace, fi) jFace[fi] = min(ny - 1, label(mesh.Cf().boundaryField()[icePatch][fi].y()/dy));
    auto setShield = [&]()      // valueFraction f = h/(h + k delta): face temperature T_f = (1-f) T_P, flux h T_f = k delta f T_P ; h >= 1e29 means bare (Dirichlet)
    {
        const scalarField& dc = Tice.patch().deltaCoeffs();
        forAll(Tice.valueFraction(), fi) { const scalar h = hShield[jFace[fi]]; Tice.valueFraction()[fi] = (h > 1e29) ? 1.0 : h/(h + kW*dc[fi]); }
        Tice.refValue() = 0.0; Tice.refGrad() = 0.0;
    };
    setShield(); T.correctBoundaryConditions();

    const scalar Ahot = gSum(mesh.magSf().boundaryField()[hotPatch]);
    auto iceFlux = [&](scalarField& qf)       // heat flux per area arriving at the ice [W/m2], per interface face
    {
        const scalarField sg = Tice.snGrad(); qf.setSize(sg.size()); forAll(sg, fi) qf[fi] = -kW*sg[fi];
    };
    auto rowDisp = [&](const scalarField& qA, const scalarField& qB, scalar dt, scalarField& dsn)      // node displacements from face-averaged fluxes; face length / dy from the current faces
    {
        scalarField dxf(ny, 0.0); const scalarField& Af = mesh.magSf().boundaryField()[icePatch]; scalarField cnt(ny, 0.0);
        forAll(qA, fi) { const label j = jFace[fi]; dxf[j] += 0.5*(qA[fi] + qB[fi])*(Af[fi]/dz)/dy*dt/(rhoI*Lf); cnt[j] += 1.0; }
        dsn.setSize(ny + 1); dsn[0] = dxf[0]; dsn[ny] = dxf[ny - 1]; for (label j = 1; j < ny; j++) dsn[j] = 0.5*(dxf[j - 1] + dxf[j]);
    };

    OFstream led(runTime.path()/"ledger.csv"); led << "t,s_mean,s_min,s_max,Vliq,Vsw_geom,Vsw_meshphi,E_sens,E_hot_cond,E_hot_adv,E_if_cond,E_latent,resid_phys,resid_disc,gcl_resid,uw" << endl;
    runTime.setTime(t0, 0); T.oldTime(); movePts(sNode);
    const scalar Vliq0 = gSum(mesh.V()); const scalar Esens0 = gSum(rhoW*cpW*T.primitiveField()*mesh.V());
    scalar Ehc = 0, Eha = 0, Eif = 0, Vmp = 0, VliqLast = Vliq0;
    scalarField qOld; iceFlux(qOld);
    Info<< "stefan2DFoam: ny = " << ny << "  cells = " << mesh.nCells() << "  densityJump = " << densityJump << endl;

    while (runTime.run())
    {
        const scalar dt = runTime.deltaTValue(); runTime++; scalarField qNew = qOld, dsn; scalar uwMean = 0;
        for (int corr = 0; corr < 3; corr++)
        {
            rowDisp(qOld, qNew, dt, dsn); scalarField sTry(sNode); forAll(sTry, j) sTry[j] += dsn[j]; movePts(sTry); setShield();
            const scalar Vliq_now = gSum(mesh.V());
            uwMean = densityJump ? (1.0 - r)*(Vliq_now - VliqLast)/dt/Ahot : 0.0;
            surfaceScalarField phi(IOobject("phiL", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh.Sf().component(vector::X)*dimensionedScalar("uw", dimVelocity, uwMean));
            surfaceScalarField phiEq(phi - mesh.phi()); phiEq.boundaryFieldRef()[icePatch] = 0.0;      // meltwater enters at T_m: NO enthalpy flux through the ice patch (with a shielded face T_f != 0 the default ALE flux would create liquid at T_f; found by the energy ledger)
            fvScalarMatrix TEqn(fvm::ddt(T) + fvm::div(phiEq, T) - fvm::laplacian(dimensionedScalar("alpha", dimViscosity, alphaW), T)); TEqn.solve(); iceFlux(qNew);
        }
        rowDisp(qOld, qNew, dt, dsn); forAll(sNode, j) sNode[j] += dsn[j];
        // accounting of this step with the converged geometry (final solve already done on the last trial geometry == this geometry up to the last corrector change; re-solve to be consistent)
        movePts(sNode); setShield(); const scalar Vliq_new = gSum(mesh.V());
        uwMean = densityJump ? (1.0 - r)*(Vliq_new - VliqLast)/dt/Ahot : 0.0;
        {
            surfaceScalarField phi(IOobject("phiL", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh.Sf().component(vector::X)*dimensionedScalar("uw", dimVelocity, uwMean));
            surfaceScalarField phiEq(phi - mesh.phi()); phiEq.boundaryFieldRef()[icePatch] = 0.0;      // meltwater enters at T_m: NO enthalpy flux through the ice patch (with a shielded face T_f != 0 the default ALE flux would create liquid at T_f; found by the energy ledger)
            fvScalarMatrix TEqn(fvm::ddt(T) + fvm::div(phiEq, T) - fvm::laplacian(dimensionedScalar("alpha", dimViscosity, alphaW), T)); TEqn.solve();
        }
        iceFlux(qOld); VliqLast = Vliq_new;
        // ledgers
        const scalar Vsw = Vliq_new - Vliq0; const scalar Esens = gSum(rhoW*cpW*T.primitiveField()*mesh.V()); const scalarField& Af = mesh.magSf().boundaryField()[icePatch];
        scalar qa = 0; forAll(qOld, fi) qa += qOld[fi]*Af[fi]; Eif += qa*dt;
        { const fvPatchScalarField& Th = T.boundaryField()[hotPatch]; const scalarField sgh = Th.snGrad(); const scalarField& Ah = mesh.magSf().boundaryField()[hotPatch]; scalar c = 0, a = 0; forAll(sgh, fi) { c += kW*sgh[fi]*Ah[fi]; a += rhoW*cpW*uwMean*Th[fi]*Ah[fi]; } Ehc += c*dt; Eha += a*dt; }
        { const scalarField& mp = mesh.phi().boundaryField()[icePatch]; Vmp += gSum(mp)*dt; }
        scalar smin = GREAT, smax = -GREAT, smean = 0; forAll(sNode, j) { smin = min(smin, sNode[j]); smax = max(smax, sNode[j]); } smean = s0 + Vsw/(Ly*dz);
        const scalar Elat = rhoI*Lf*Vsw;
        led << runTime.value() << "," << smean << "," << smin << "," << smax << "," << Vliq_new << "," << Vsw << "," << Vmp << "," << Esens << "," << Ehc << "," << Eha << "," << Eif << "," << Elat << ","
            << ((Esens - Esens0) + Elat - (Ehc + Eha)) << "," << ((Esens - Esens0) + Eif - (Ehc + Eha)) << "," << (Vsw - Vmp) << "," << uwMean << endl;
        if (runTime.writeTime()) { runTime.write(); OFstream ic(runTime.path()/runTime.timeName()/"interface.csv"); ic << "y,s" << endl; for (label j = 0; j <= ny; j++) ic << j*dy << "," << sNode[j] << endl; }
    }
    return 0;
}
