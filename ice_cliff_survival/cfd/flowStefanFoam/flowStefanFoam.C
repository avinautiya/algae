/*---------------------------------------------------------------------------*\
  flowStefanFoam (ice_cliff_survival m2): laminar single-phase incompressible flow + conservative sensible-heat transport in a RIGID-LID channel whose floor is ONE melting (temperate) ice boundary that moves.
  Geometry (2-D, one cell deep, depth dz): 0<x<Ltot, floor y=b(x,t) (initially 0), lid y=H fixed. Floor = [inert wall Lu | melting ice Li | inert wall Ld]; inert walls keep inlet/outlet openings fixed (inlet Poiseuille with fixed Q, outlet p=0, T zero-gradient).
  Equations (ALE, absolute velocity u, mesh flux w.n = meshPhi):  div u = 0 ;  d/dt u + div((u-w) u) = -grad p/rho + nu lap u ;  d/dt T + div((u-w) T) = alpha lap T  (finite volume, Euler in time, conservative: V0->V handled by ddt).
  Interface (n from ice into liquid, outward face normal Sf points into ice, ice at rest, T_m = 0):  mass:  rho_w (u_w - v).Sf = rho_i (0 - v).Sf  =>  u_w.Sf = (1 - rho_i/rho_w) meshPhi_f   PER FACE (local, exact; no global uniform inflow);  equal density: u_w = 0.
  energy: rho_i L m_n = q = k_w * (-snGrad T_f) (same discrete face flux the T-equation removes); meltwater is created at T_m: NO enthalpy flux through the ice patch (phiEq = 0 there).
  Kinematics: floor is a graph; node displacement = mean of adjacent face displacements dxf = q A_f dt /(rho_i L Dx dz) (total swept volume equals sum of face melt volumes EXACTLY); interior nodes stretch linearly between floor and lid.
  Time coupling: BACKWARD-EULER FIXED POINT per step: geometry_{k+1} = geometry_n + dt m(q_k); at each iterate the flow+heat are re-solved from the step-start state; iterate to |d geometry| < tolCoupling*scale. Coupling defect recorded.
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include <vector>
#include <algorithm>
#include "OFstream.H"
#include "mixedFvPatchFields.H"
#include "fixedValueFvPatchFields.H"

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    IOdictionary prop(IOobject("flowProperties", runTime.constant(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE));
    const scalar rhoW = readScalar(prop.lookup("rhoW")), cpW = readScalar(prop.lookup("cpW")), kW = readScalar(prop.lookup("kW")), nu = readScalar(prop.lookup("nu"));
    const scalar rhoI = readScalar(prop.lookup("rhoI")), Lf = readScalar(prop.lookup("L")), H = readScalar(prop.lookup("H")), dz = readScalar(prop.lookup("dz"));
    const scalar Lu = readScalar(prop.lookup("Lu")), Li = readScalar(prop.lookup("Li")), Ltot = readScalar(prop.lookup("Ltot")), Ubar = readScalar(prop.lookup("Ubar")), Tin = readScalar(prop.lookup("Tin"));
    const bool melt = readBool(prop.lookup("melt")), closedEnds = readBool(prop.lookup("closedEnds")), lidClosure = readBool(prop.lookup("lidClosure")), dump = readBool(prop.lookup("dump"));
    const label nx = readLabel(prop.lookup("nx")), ny = readLabel(prop.lookup("ny")), maxCoup = readLabel(prop.lookup("maxCoupling")), nPiso = readLabel(prop.lookup("nPiso"));
    const bool useDdtCorr = prop.found("ddtCorr") ? readBool(prop.lookup("ddtCorr")) : true;
    const bool jumpBC = prop.found("jumpBC") ? readBool(prop.lookup("jumpBC")) : true;
    const scalar slopeMax = prop.found("slopeMax") ? readScalar(prop.lookup("slopeMax")) : 1.0;
    const label dumpEvery = prop.found("dumpEvery") ? readLabel(prop.lookup("dumpEvery")) : 1;
    const bool convU = prop.found("convU") ? readBool(prop.lookup("convU")) : true;
    const bool momPred = prop.found("momPred") ? readBool(prop.lookup("momPred")) : true;
    const scalar tolC = readScalar(prop.lookup("tolCoupling")), r = rhoI/rhoW, alpha = kW/(rhoW*cpW), Dx = Ltot/nx, hBare = 1e30; (void)hBare;
    const label icePatch = mesh.boundaryMesh().findPatchID("ice"), inletP = mesh.boundaryMesh().findPatchID("inlet"), outletP = mesh.boundaryMesh().findPatchID("outlet"), lidP = mesh.boundaryMesh().findPatchID("lid");
    const label fuP = mesh.boundaryMesh().findPatchID("floorUp"), fdP = mesh.boundaryMesh().findPatchID("floorDown"); labelList kinP; kinP.append(icePatch); if (fuP >= 0) kinP.append(fuP); if (fdP >= 0) kinP.append(fdP);      // faces that can move: ice + the junction faces of the inert walls (shared nodes)
    const scalarField hShield(prop.found("hShield") ? scalarField(prop.lookup("hShield")) : scalarField(nx, 1e30));
    if (icePatch < 0 || inletP < 0 || outletP < 0 || lidP < 0) FatalErrorInFunction << "patches ice, inlet, outlet, lid required" << exit(FatalError);

    // reference geometry: column of each point (x/Dx), eta = y/H
    pointField p0(mesh.points()); labelList iCol(p0.size()); scalarField eta(p0.size());
    forAll(p0, pi) { iCol[pi] = label(p0[pi].x()/Dx + 0.5); eta[pi] = p0[pi].y()/H; }
    scalarField bNode(nx + 1, 0.0);
    auto movePts = [&](const scalarField& bn) { pointField pp(p0); forAll(pp, pi) pp[pi].y() = bn[iCol[pi]] + (H - bn[iCol[pi]])*eta[pi]; mesh.movePoints(pp); };
    // cell (i,j) from reference centres
    labelList cI(mesh.nCells()), cJ(mesh.nCells());      // column i = x/Dx, row j = rank of the cell centre in y within its column (reference mesh)
    { std::vector<std::vector<std::pair<double, int>>> col(nx); forAll(cI, ci) { cI[ci] = min(nx - 1, label(mesh.C()[ci].x()/Dx)); col[cI[ci]].push_back({mesh.C()[ci].y(), ci}); }
      for (auto& c : col) { std::sort(c.begin(), c.end()); for (size_t k = 0; k < c.size(); k++) cJ[c[k].second] = label(k); } }

    volScalarField p(IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    volVectorField U(IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    volScalarField T(IOobject("T", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    // initial: developed Poiseuille profile in every cell, uniform T = Tin; inlet profile
    forAll(U, ci) { const scalar e = mesh.C()[ci].y()/H; U[ci] = vector(6.0*Ubar*e*(1.0 - e), 0, 0); T[ci] = Tin; }
    if (prop.found("neumannTw"))      // Stefan-limit initial condition: similarity profile measured from the (hot, fixed-T) lid at distance x = H - y, taken at time t0 (liquid depth H = 2 lam sqrt(alpha t0))
    {
        const scalar Tw = readScalar(prop.lookup("neumannTw")), lam = readScalar(prop.lookup("neumannLam")), bb = readScalar(prop.lookup("neumannB")), t0n = readScalar(prop.lookup("neumannT0"));
        auto Iint = [&](scalar e, scalar b0) { return 0.5*std::sqrt(constant::mathematical::pi)*(::erf(e - b0) + ::erf(b0)); };
        forAll(T, ci) { const scalar x = H - mesh.C()[ci].y(); T[ci] = Tw - Tw/Iint(lam, bb)*Iint(x/(2*std::sqrt(alpha*t0n)), bb); }
    }
    if (closedEnds) U = dimensionedVector("z", dimVelocity, vector::zero);
    { fvPatchVectorField& Ui = U.boundaryFieldRef()[inletP]; forAll(Ui, fi) { const scalar e = mesh.Cf().boundaryField()[inletP][fi].y()/H; Ui[fi] = closedEnds ? vector::zero : vector(6.0*Ubar*e*(1.0 - e), 0, 0); } }
    U.correctBoundaryConditions(); T.correctBoundaryConditions();
    surfaceScalarField phi(IOobject("phi", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), fvc::flux(U));
    mixedFvPatchScalarField& Tice = refCast<mixedFvPatchScalarField>(T.boundaryFieldRef()[icePatch]);
    fixedValueFvPatchVectorField& Uice = refCast<fixedValueFvPatchVectorField>(U.boundaryFieldRef()[icePatch]);
    fixedValueFvPatchVectorField& Ulid = refCast<fixedValueFvPatchVectorField>(U.boundaryFieldRef()[lidP]);
    labelList iceCol(Tice.size()); forAll(iceCol, fi) iceCol[fi] = min(nx - 1, label(mesh.Cf().boundaryField()[icePatch][fi].x()/Dx));
    labelList iceColFirstLast(2); iceColFirstLast[0] = label(Lu/Dx + 0.5); iceColFirstLast[1] = label((Lu + Li)/Dx + 0.5);
    auto setShield = [&]() { const scalarField& dc = Tice.patch().deltaCoeffs(); forAll(Tice.valueFraction(), fi) { const scalar h = hShield[iceCol[fi]]; Tice.valueFraction()[fi] = (h > 1e29) ? 1.0 : h/(h + kW*dc[fi]); } Tice.refValue() = 0.0; Tice.refGrad() = 0.0; };
    setShield(); T.correctBoundaryConditions();
    const scalar Alid = gSum(mesh.magSf().boundaryField()[lidP]);
    auto iceFlux = [&](scalarField& qf) { const scalarField sg = Tice.snGrad(); qf.setSize(sg.size()); forAll(sg, fi) qf[fi] = -kW*sg[fi]; };
    auto faceDisp = [&](const scalarField& q, scalar dt, scalarField& bn)     // new floor nodes from face fluxes q on the CURRENT geometry (face area)
    {
        scalarField dxf(nx, 0.0); const scalarField& Af = mesh.magSf().boundaryField()[icePatch];
        if (melt) forAll(q, fi) dxf[iceCol[fi]] += q[fi]*(Af[fi]/dz)/Dx*dt/(rhoI*Lf);       // negative y displacement magnitude: floor moves DOWN by dxf
        bn = bNode; scalarField dn(nx + 1, 0.0); dn[0] = dxf[0]; dn[nx] = dxf[nx - 1]; for (label i = 1; i < nx; i++) dn[i] = 0.5*(dxf[i - 1] + dxf[i]);
        forAll(bn, i) bn[i] -= dn[i];
    };
    const scalar nuS = nu;
    OFstream led(runTime.path()/"ledger.csv");
    led << "t,dt,nCoup,coupRes,Vliq,Vsw,Vsw_face,gcl,Qin,Qout,Qfloor,Qlid,Esens,Ein,Eout,Elid,Efl,Eif,Elat,R_disc,R_phys,massLocal,qmax,umax" << endl;
    OFstream dmp(runTime.path()/"dump.csv"); if (dump) dmp << "# dump: per step blocks. S n t | B bNodes | C i j T Ux Uy | I inlet(j,Ux,Uy,T) | O outlet | F floor(i,Ux,Uy,Tf) | M floor meshPhi_dt per face" << endl;
    runTime.setTime(0.0, 0); U.oldTime(); T.oldTime(); phi.oldTime(); p.oldTime(); movePts(bNode);
    const scalar Vliq0 = gSum(mesh.V()); scalar VliqLast = Vliq0, Esens0 = gSum(rhoW*cpW*T.primitiveField()*mesh.V()); scalar Eif = 0, Ein = 0, Eout = 0, Elid = 0, Efl = 0;
    scalarField qOld; iceFlux(qOld);
    Info<< "flowStefanFoam: cells " << mesh.nCells() << " melt " << melt << " density ratio " << r << endl;
    label pRefCell = 0; scalar pRefValue = 0.0;
    while (runTime.run())
    {
        const scalar dt = runTime.deltaTValue(); volVectorField Un("Un", U); volScalarField pn("pn", p), Tn("Tn", T); surfaceScalarField phin("phin", phi);
        runTime++;
        scalarField qNew = qOld, bUsed, bProp, rVec, rPrev; scalar coupRes = 0, omega = 1.0, scale = 1e-300; label k = 0; faceDisp(qOld, dt, bUsed);      // k = 0: explicit predictor from the previous flux
        for (k = 0; k < (melt ? maxCoup : 1); k++)
        {
            movePts(bUsed); setShield();
            U = Un; p = pn; T = Tn; phi = phin;
            forAll(kinP, kp) { const label pp = kinP[kp]; const scalarField& mp = mesh.phi().boundaryField()[pp]; const vectorField& Sf = mesh.Sf().boundaryField()[pp];
              vectorField& Ub = refCast<fixedValueFvPatchVectorField>(U.boundaryFieldRef()[pp]); forAll(Ub, fi) Ub[fi] = (melt && !(r == 1.0) && jumpBC) ? (1.0 - r)*mp[fi]/magSqr(Sf[fi])*Sf[fi] : vector::zero; }
            if (lidClosure)      // TEST-ONLY global closure at the LID (not at the interface): inflow through the lid equals the net floor outflow N = sum (1-r) meshPhi_b
            { scalar N = 0; forAll(kinP, kp) { const scalarField& mp = mesh.phi().boundaryField()[kinP[kp]]; forAll(mp, fi) N += (1.0 - r)*mp[fi]; }
              const vectorField& Sfl = mesh.Sf().boundaryField()[lidP]; vectorField& Ul = Ulid; forAll(Ul, fi) Ul[fi] = -(N/Alid)*Sfl[fi]/mag(Sfl[fi]); }
            U.correctBoundaryConditions();
            surfaceScalarField phiRel("phiRel", phi - mesh.phi());
            fvVectorMatrix UEqn(fvm::ddt(U) + fvm::div(convU ? phiRel : surfaceScalarField("phiRel", 0.0*phiRel), U) - fvm::laplacian(dimensionedScalar("nu", dimViscosity, nuS), U));
            if (momPred) solve(UEqn == -fvc::grad(p));
            for (int c = 0; c < nPiso; c++)
            {
                volScalarField rAU(1.0/UEqn.A()); volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
                surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA) + (useDdtCorr ? fvc::interpolate(rAU)*fvc::ddtCorr(U, phi) : surfaceScalarField::New("zeroCorr", mesh, dimensionedScalar("z", dimVelocity*dimArea, 0.0))));
                constrainPressure(p, U, phiHbyA, rAU);
                fvScalarMatrix pEqn(fvm::laplacian(rAU, p) == fvc::div(phiHbyA)); if (closedEnds) pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve(); phi = phiHbyA - pEqn.flux(); U = HbyA - rAU*fvc::grad(p); U.correctBoundaryConditions();
            }
            surfaceScalarField phiEq("phiEq", phi - mesh.phi()); forAll(kinP, kp) phiEq.boundaryFieldRef()[kinP[kp]] = 0.0;      // new liquid enters at T_m: no enthalpy flux through any face that can move
            fvScalarMatrix TEqn(fvm::ddt(T) + fvm::div(phiEq, T) - fvm::laplacian(dimensionedScalar("alpha", dimViscosity, alpha), T)); TEqn.solve();
            iceFlux(qNew);
            if (!melt) { k++; break; }
            faceDisp(qNew, dt, bProp); rVec = bProp - bUsed; coupRes = 0; forAll(rVec, i) coupRes = max(coupRes, mag(rVec[i]));
            if (k == 0) { forAll(bProp, i) scale = max(scale, mag(bProp[i] - bNode[i])); }
            if (runTime.timeIndex() == 1) { scalar qs2 = 0; forAll(qNew, fi) qs2 += qNew[fi]*mesh.magSf().boundaryField()[icePatch][fi]; scalar vv = gSum(mesh.V()); Info<< "coup " << k << " res " << coupRes << " omega " << omega << " sumqA " << qs2 << " V " << vv << " qmid " << qNew[qNew.size()/2] << " q0 " << qNew[0] << endl; }
            if (coupRes < tolC*scale) { k++; break; }
            if (k > 0) { scalar den = 0, num = 0; forAll(rVec, i) { const scalar dri = rVec[i] - rPrev[i]; den += dri*dri; num += rPrev[i]*dri; } if (den > 0) omega = max(0.05, min(1.5, -omega*num/den)); }
            rPrev = rVec; bUsed = bUsed + omega*rVec;
        }
        const scalar defect = coupRes;
        bNode = bUsed; qOld = qNew;
        { scalar smax = 0; for (label i = 0; i < nx; i++) smax = max(smax, mag(bNode[i+1] - bNode[i])/Dx); if (smax > slopeMax) FatalErrorInFunction << "floor slope " << smax << " exceeds slopeMax " << slopeMax << ": the vertical-column mesh is no longer valid (leading-edge scarp); run aborted at t = " << runTime.value() << exit(FatalError); }
        // ---- ledgers (this step's final state) ----
        const scalar Vl = gSum(mesh.V()), Vsw = Vl - Vliq0, dV = Vl - VliqLast; const scalarField& Af = mesh.magSf().boundaryField()[icePatch];
        scalar qa = 0, qmax = 0, Vface = 0; forAll(qNew, fi) { qa += qNew[fi]*Af[fi]; qmax = max(qmax, qNew[fi]); Vface += (melt ? qNew[fi]*Af[fi]*dt/(rhoI*Lf) : 0.0); }
        Eif += qa*dt;
        scalar Qin = 0, Qout = 0, Qfl = 0, Qlid = 0, ein = 0, eout = 0, elid = 0, gclSum = 0;
        { const scalarField& pf = phi.boundaryField()[inletP]; const fvPatchScalarField& Tp = T.boundaryField()[inletP]; const scalarField sg = Tp.snGrad(); const scalarField& Ab = mesh.magSf().boundaryField()[inletP]; forAll(pf, fi) { Qin -= pf[fi]; ein -= rhoW*cpW*pf[fi]*Tp[fi]; ein += kW*sg[fi]*Ab[fi]; } }      // inlet: advective + diffusive heat in
        { const scalarField& pf = phi.boundaryField()[outletP]; const fvPatchScalarField& Tp = T.boundaryField()[outletP]; const scalarField sg = Tp.snGrad(); const scalarField& Ab = mesh.magSf().boundaryField()[outletP]; forAll(pf, fi) { Qout += pf[fi]; eout += rhoW*cpW*pf[fi]*Tp[fi]; eout -= kW*sg[fi]*Ab[fi]; } }      // outlet: advective + diffusive heat out
        forAll(kinP, kp) Qfl += gSum(phi.boundaryField()[kinP[kp]]); Qlid = gSum(phi.boundaryField()[lidP]);
        { const fvPatchScalarField& Tp = T.boundaryField()[lidP]; const scalarField sg = Tp.snGrad(); const scalarField& Ab = mesh.magSf().boundaryField()[lidP]; forAll(sg, fi) elid += kW*sg[fi]*Ab[fi]; }      // heat INTO the domain through the lid (snGrad outward = T_f - T_P)
        forAll(mesh.boundary(), pi) gclSum += gSum(mesh.phi().boundaryField()[pi]);
        scalar efl = 0; forAll(kinP, kp) { if (kinP[kp] == icePatch) continue; const fvPatchScalarField& Tp = T.boundaryField()[kinP[kp]]; const scalarField sg = Tp.snGrad(); const scalarField& Ab = mesh.magSf().boundaryField()[kinP[kp]]; forAll(sg, fi) efl += kW*sg[fi]*Ab[fi]; }      // heat INTO the domain through inert (isothermal-T_m or adiabatic) floor walls
        { const scalarField& pf = phi.boundaryField()[lidP]; const fvPatchScalarField& Tp = T.boundaryField()[lidP]; forAll(pf, fi) ein -= rhoW*cpW*pf[fi]*Tp[fi]; }      // advective heat through the lid (inflow in the lid-closure test)
        Ein += ein*dt; Eout += eout*dt; Elid += elid*dt; Efl += efl*dt;
        const scalar Esens = gSum(rhoW*cpW*T.primitiveField()*mesh.V()), Elat = rhoI*Lf*Vsw;       // lat uses ice volume lost == swept volume
        // mass: local interface residual  max |phi_b - (1-r) meshPhi_b| / (|meshPhi_b| + tiny) ; global: Qout - Qin - Qfloor-type
        scalar mloc = 0, msc = 1e-300; forAll(kinP, kp) { const scalarField& mp = mesh.phi().boundaryField()[kinP[kp]]; const scalarField& pb = phi.boundaryField()[kinP[kp]]; forAll(pb, fi) { mloc = max(mloc, mag(pb[fi] - ((melt && r != 1.0) ? (1.0 - r)*mp[fi] : 0.0))); msc = max(msc, mag(mp[fi])); } } mloc /= msc;
        led << runTime.value() << "," << dt << "," << k << "," << defect << "," << Vl << "," << Vsw << "," << Vface << "," << (dV - dt*gclSum) << "," << Qin << "," << Qout << "," << Qfl << "," << Qlid << "," << Esens << "," << Ein << "," << Eout << "," << Elid << "," << Efl << "," << Eif << "," << Elat << ","
            << ((Esens - Esens0) + Eif - (Ein - Eout + Elid + Efl)) << "," << ((Esens - Esens0) + Elat - (Ein - Eout + Elid + Efl)) << "," << mloc << "," << qmax << "," << max(mag(U)).value() << endl;
        VliqLast = Vl;
        if (dump && (runTime.timeIndex() % dumpEvery == 0 || runTime.value() > runTime.endTime().value() - 0.5*dt || runTime.timeIndex() == 1))
        {
            dmp << "S " << runTime.timeIndex() << " " << runTime.value() << endl; dmp << "B"; forAll(bNode, i) dmp << " " << bNode[i]; dmp << endl;
            forAll(U, ci) dmp << "C " << cI[ci] << " " << cJ[ci] << " " << T[ci] << " " << U[ci].x() << " " << U[ci].y() << endl;
            { const vectorField& Ub = U.boundaryField()[inletP]; const scalarField& Tb = T.boundaryField()[inletP]; forAll(Ub, fi) dmp << "I " << fi << " " << Ub[fi].x() << " " << Ub[fi].y() << " " << Tb[fi] << endl; }
            { const vectorField& Ub = U.boundaryField()[outletP]; const scalarField& Tb = T.boundaryField()[outletP]; forAll(Ub, fi) dmp << "O " << fi << " " << Ub[fi].x() << " " << Ub[fi].y() << " " << Tb[fi] << endl; }
            { const vectorField& Ub = U.boundaryField()[icePatch]; const scalarField& Tb = T.boundaryField()[icePatch]; const scalarField& mp = mesh.phi().boundaryField()[icePatch]; forAll(Ub, fi) dmp << "F " << iceCol[fi] << " " << Ub[fi].x() << " " << Ub[fi].y() << " " << Tb[fi] << " " << mp[fi]*dt << " " << qNew[fi] << endl; }
        }
        if (runTime.writeTime()) runTime.write();
    }
    return 0;
}
