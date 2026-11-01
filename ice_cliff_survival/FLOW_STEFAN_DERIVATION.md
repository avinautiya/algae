# Flow–heat–moving-ice coupling: derivation and discretisation (M2, `cfd/flowStefanFoam`)

**Scope.** One single-phase liquid (water), rigid lid, one melting temperate-ice boundary (the floor), prescribed inlet flow and temperature, documented outlet; initially bare ice. NOT included: mobile grains, free surface, fracture, atmosphere, glacier mechanics, buoyancy, turbulence, cold ice, debris storage. Laminar only (inputs below are assumed test values, not measurements).

## 1. Governing equations on the moving domain Ω(t)
Constant ρ_w, c_w, k_w, ν. In the lab frame with absolute velocity u and mesh (domain-boundary) velocity w:
- continuity: ∇·u = 0 on Ω(t) (domain growth is fed by phase change at the interface, not by compressing the liquid);
- momentum (ALE, conservative): d/dt∫_V u dV + ∮_{∂V} (u−w)·n u dS = −(1/ρ)∮ p n dS + ν∮ ∇u·n dS;
- sensible heat (T measured from T_m = 0): d/dt∫_V T dV + ∮ (u−w)·n T dS = α∮ ∇T·n dS.
Reynolds transport gives the ALE form for any control volume V(t) moving with w. The discrete form uses the geometric conservation law (GCL): V^{n+1} − V^n = Δt Σ_f w_f·S_f, with the face mesh flux `meshPhi_f` equal to the volume swept by face f over the step (computed from old/new points; checked each step: `gcl` column, ≲1e-19 m³).

## 2. Interface conditions (n from ice into liquid; S_f points out of the liquid into the ice)
Ice at rest (u_ice=0), interface velocity v. Mass: ρ_w(u_w−v)·n = ρ_i(0−v)·n ⇒ **u_w·n = (1−ρ_i/ρ_w)(v·n)** pointwise (equal density: u_w·n=0). With S_f: **φ_f ≡ u_w·S_f = (1−r) meshPhi_f**, r=ρ_i/ρ_w, imposed FACE BY FACE as the velocity boundary condition; tangential velocity 0 (no slip on stationary ice). The relative mass flux through the face is φ_f − meshPhi_f = −r meshPhi_f, i.e. exactly ρ_i/ρ_w of the swept volume appears as liquid (ρ_i-mass of ice melted). No uniform/global closure is used at the interface. (The 1-D unequal-density verification needs inflow at the *lid*; that lid condition is a TEST-ONLY global closure and is flagged as such.)
Energy: ρ_i L m_n = q, q = −k ∂T/∂n|_liquid side (temperate ice carries no conductive flux). Meltwater is created at T_m: the enthalpy flux of the created liquid is c_w(T_m−T_m)=0, so **no enthalpy crosses the interface patch** (also for the inert-wall junction faces, which can move because they share nodes with ice faces).
Discrete interface heat is the SAME face flux the temperature equation removes (q_f = −k snGrad T_f), and the displacement uses it: V_f = q_f A_f Δt /(ρ_i L).

## 3. Mesh motion
Floor is a graph y=b(x,t); columns stretch linearly between floor and lid (lid fixed): y = b_i + (H−b_i)η. Node displacement = mean of adjacent face displacements dxf_f = V_f/(Δx dz); then Σ node-trapezoid volumes = Σ V_f exactly (total swept volume = total melted volume; local smoothing error O(h)). Inlet and outlet openings sit on inert walls (Lu, Ld) so their cross-sections are fixed (inlet Poiseuille profile ⇒ fixed Q).
Validity limit: vertical columns require moderate floor slope; the code aborts with an explicit error when max |Δb/Δx| > slopeMax (default 1.0, i.e. 45°). A finite-height scarp forms where the melting floor meets the non-melting inert wall and steepens under refinement (slope ∝ 1/Δx); the guard is an admissibility check, not a tuned parameter.

## 4. Time coupling (backward Euler, consistent)
Per step, fixed-point on the geometry b^{n+1}: (i) predictor from q^n; (ii) restore step-start fields, move mesh to b_k, solve U–p (PISO, 3 correctors) and T; (iii) b_{k+1}^{prop} = b^n − Δ(q_k); residual r_k = b^{prop} − b_k; Aitken-relaxed update; stop at max|r| < tol·(step displacement). At convergence the heat removed from the water, the latent energy ρ_i L V_swept and the interface displacement use the same q, so `E_if − E_latent` is the coupling defect (reported, never absorbed into a ledger).

## 5. Energy accounting
Discrete (equation) residual: ΔE_sens + E_if − (E_in − E_out + E_lid + E_inert-floor) with advective + diffusive inlet/outlet terms (≈1e-13 of the terms). Physical residual replaces E_if by ρ_i L V_swept (geometry): its difference to the discrete one is exactly the coupling defect. Independent reconstruction (`flowstefan_post.py`) recomputes everything from saved fields and node heights with its own geometry.
