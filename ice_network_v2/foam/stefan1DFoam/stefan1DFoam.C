/*---------------------------------------------------------------------------*\
  stefan1DFoam (v2): sharp-interface one-phase Stefan melting on a MOVING (ALE) 1-D mesh, built on OpenFOAM v1912 finiteVolume.
  Liquid occupies 0 < x < s(t) (hot boundary 'hot' at x = 0 held at T_w, supplying liquid when the density jump requires it); the ice (temperate, T = T_m = 0, semi-infinite, at rest) fills x > s.
  Equation (ALE, conservative):  d/dt int T dV  +  int T (u_w - w) . dS  -  int alpha grad T . dS = 0,   u_w = uniform liquid velocity, w = mesh velocity (mesh flux phiMesh),
  Stefan condition (energy):  rho_i L sdot = - k dT/dx |_{x = s}   (heat conducted from the liquid into the interface; ice has no conduction so no ice-side term),
  Mass (normal convention n points from ice into liquid, i.e. -x):  rho_w (u_w - sdot) = rho_i (0 - sdot)  =>  u_w = (1 - rho_i/rho_w) sdot  (liquid velocity in +x, toward the interface; in the packet's convention u_w.n = (rho_i/rho_w - 1) m_n with m_n = sdot).
  Interface motion: the 'ice' patch moves with sdot; interior points are stretched linearly (x = xi s); Heun predictor-corrector for s (two evaluations of the interface heat flux per step).
  Interface heat flux uses the second-order one-sided wall gradient through the first two cell values. Ledger written every step: sensible storage, hot-boundary conductive+advective supply, latent energy rho_i L (s - s0).
  Properties are the physical ones (rho_w, cp_w, k_w, rho_i, L) read from constant/stefanProperties; nothing is rescaled.
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "OFstream.H"

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary prop(IOobject("stefanProperties", runTime.constant(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE));
    const scalar rhoW = readScalar(prop.lookup("rhoW")), cpW = readScalar(prop.lookup("cpW")), kW = readScalar(prop.lookup("kW"));
    const scalar rhoI = readScalar(prop.lookup("rhoI")), Lf = readScalar(prop.lookup("L")), Tw = readScalar(prop.lookup("Tw"));
    const scalar s0 = readScalar(prop.lookup("s0"));            // initial liquid thickness (initial profile = exact solution at t0)
    const scalar t0 = readScalar(prop.lookup("t0"));
    const scalar lam = readScalar(prop.lookup("lambda")), bpar = readScalar(prop.lookup("bpar"));   // similarity constants of the exact solution (for initial profile only)
    const bool densityJump = readBool(prop.lookup("densityJump"));
    const scalar alphaW = kW/(rhoW*cpW);
    const scalar r = rhoI/rhoW;

    // exact similarity profile T(eta) = Tw - C int_0^eta exp(-(xi-b)^2) dxi ; Tw ... 0 at eta = lam
    auto Iint = [&](scalar eta) { return 0.5*std::sqrt(constant::mathematical::pi)*(::erf(eta - bpar) + ::erf(bpar)); };
    const scalar Cc = Tw/Iint(lam);

    volScalarField T(IOobject("T", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    const label icePatch = mesh.boundaryMesh().findPatchID("ice");
    const label hotPatch = mesh.boundaryMesh().findPatchID("hot");
    if (icePatch < 0 || hotPatch < 0) FatalErrorInFunction << "patches 'ice' and 'hot' required" << exit(FatalError);

    // reference (unit) coordinates: xi = x / s0 of the initial mesh; moving mesh points x = xi * s
    pointField xi0(mesh.points()); forAll(xi0, pi) { xi0[pi].x() /= s0; }     // only the x-coordinate is stretched (y,z fixed)
    scalar s = s0;
    // initial profile from the exact solution at time t0 (cell centres)
    forAll(T, ci) { const scalar eta = mesh.C()[ci].x()/(2*std::sqrt(alphaW*t0)); T[ci] = Tw - Cc*Iint(eta); }
    T.correctBoundaryConditions();

    const scalar Aface = gSum(mesh.magSf().boundaryField()[hotPatch]);                          // cross-section area [m2] (ledger per unit area)
    // interface heat flux used for the Stefan condition. MODE 'discrete' (default): the SAME flux the temperature equation removes through the interface face, k (T_P - T_b) delta  (first-order in dx, exactly energy-conserving with the latent term);
    // MODE 'quadratic': second-order one-sided gradient through the two nearest cells (more accurate flux, not exactly conservative). Selected by stefanProperties entry fluxMode (word).
    const word fluxMode(prop.lookupOrDefault<word>("fluxMode", "discrete"));
    auto ifaceFlux = [&]() -> scalar
    {
        if (fluxMode == "discrete") { return -kW*T.boundaryField()[icePatch].snGrad()()[0]; }
        const label N = mesh.nCells(); const scalar x1 = s - mesh.C()[N-1].x(), x2 = s - mesh.C()[N-2].x(); const scalar t1 = T[N-1], t2 = T[N-2];
        const scalar d0 = (t1*x2*x2 - t2*x1*x1)/(x1*x2*(x2 - x1));
        return kW*d0;
    };
    OFstream led(runTime.path()/"ledger.csv"); led << "t,s,sdot,q_if,q_hot_cond,q_hot_adv,q_if_cond_discrete,E_sens,E_latent,E_supply_cum,E_cond_out_if_cum,resid_physical,resid_discrete" << endl;
    scalar Esupply = 0, Eout = 0; const scalar Esens0 = gSum(rhoW*cpW*T.primitiveField()*mesh.V())/Aface;
    runTime.setTime(t0, 0); T.oldTime();
    Info<< "stefan1DFoam: lambda = " << lam << "  Ste = " << cpW*Tw/Lf << "  s0 = " << s0 << "  densityJump = " << densityJump << endl;
    scalar sdot = ifaceFlux()/(rhoI*Lf);
    { pointField pp(xi0); forAll(pp, pi) { pp[pi].x() *= s; } mesh.movePoints(pp); }   // consistent start

    while (runTime.run())
    {
        const scalar dt = runTime.deltaTValue(); runTime++;
        scalar sd1 = sdot, sNew = s + dt*sdot, uwMean = 0;
        for (int corr = 0; corr < 3; corr++)
        {
            { pointField pp(xi0); forAll(pp, pi) { pp[pi].x() *= sNew; } mesh.movePoints(pp); }                                                        // V0/old points are those of the start of this time step (set at the first movePoints of the step)
            uwMean = densityJump ? (1.0 - r)*0.5*(sdot + sd1) : 0.0;
            surfaceScalarField phi(IOobject("phiL", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh.Sf().component(vector::X)*dimensionedScalar("uw", dimVelocity, uwMean));
            fvScalarMatrix TEqn(fvm::ddt(T) + fvm::div(phi - mesh.phi(), T) - fvm::laplacian(dimensionedScalar("alpha", dimViscosity, alphaW), T));
            TEqn.solve();
            sd1 = ifaceFlux()/(rhoI*Lf); sNew = s + 0.5*dt*(sdot + sd1);
        }
        // final geometry with the converged interface speed, final T solve
        { pointField pp(xi0); forAll(pp, pi) { pp[pi].x() *= sNew; } mesh.movePoints(pp); } uwMean = densityJump ? (1.0 - r)*0.5*(sdot + sd1) : 0.0;
        {
            surfaceScalarField phi(IOobject("phiL", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh.Sf().component(vector::X)*dimensionedScalar("uw", dimVelocity, uwMean));
            fvScalarMatrix TEqn(fvm::ddt(T) + fvm::div(phi - mesh.phi(), T) - fvm::laplacian(dimensionedScalar("alpha", dimViscosity, alphaW), T));
            TEqn.solve();
        }
        s = sNew; sdot = ifaceFlux()/(rhoI*Lf);
        const scalar q = ifaceFlux();
        const scalar Esens = gSum(rhoW*cpW*T.primitiveField()*mesh.V())/Aface, Elat = rhoI*Lf*(s - s0);
        const scalar qhotc = kW*T.boundaryField()[hotPatch].snGrad()()[0] * (-1.0) * (-1.0);        // inflow through the hot face: k (T_b - T_P) delta  (snGrad = (T_b - T_P) delta along the outward normal)
        const scalar qhota = rhoW*cpW*uwMean*Tw;
        const scalar qifc = kW*(-T.boundaryField()[icePatch].snGrad()()[0]);                          // discrete (first-order) conductive loss through the interface face used by the matrix
        Esupply += (qhotc + qhota)*dt; Eout += qifc*dt;
        led << runTime.value() << "," << s << "," << sdot << "," << q << "," << qhotc << "," << qhota << "," << qifc << "," << Esens << "," << Elat << "," << Esupply << "," << Eout << ","
            << ((Esens - Esens0) + Elat - Esupply) << "," << ((Esens - Esens0) + Eout - Esupply) << endl;
        if (runTime.writeTime()) runTime.write();
    }
    return 0;
}
