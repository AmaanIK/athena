//========================================================================================
// Kelvin–Helmholtz (hydro) with a passive scalar and a sinusoidal seed at the interface.
// - Two layers: (rho_lo, vx_lo) below, (rho_hi, vx_hi) above, uniform p0.
// - Seed options:
//    a) velocity: vy = A_vy * cs_ref * sin(kx) * exp[-(y-ymid)^2/(2*sigma^2)]
//    b) interface displacement: interface at y = ymid + eta * sin(kx)
//========================================================================================
#include <cmath>
#include <sstream>
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#if MAGNETIC_FIELDS_ENABLED
#include "../field/field.hpp"
#endif
#if NSCALARS > 0
#include "../scalars/scalars.hpp"
#endif

void Mesh::InitUserMeshData(ParameterInput *pin) { return; }

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // --- Layer properties (classic KH) ---
  const Real rho_lo = pin->GetOrAddReal("problem","rho_lo", 1.0);    // lower layer (y<ymid)
  const Real rho_hi = pin->GetOrAddReal("problem","rho_hi", 2.0);    // upper layer (y>=ymid)
  const Real p0     = pin->GetOrAddReal("problem","p0",     2.5);    // uniform pressure
  const Real vx_lo  = pin->GetOrAddReal("problem","vx_lo", -0.5);
  const Real vx_hi  = pin->GetOrAddReal("problem","vx_hi", +0.5);

  // --- Passive scalar colors (mass fractions) ---
  const Real color_lo = pin->GetOrAddReal("problem","color_lo", 0.0);
  const Real color_hi = pin->GetOrAddReal("problem","color_hi", 1.0);

  // --- Sinusoidal seed controls ---
  const std::string seed_kind = pin->GetOrAddString("problem","seed_kind","velocity");
  const int    m_mode   = pin->GetOrAddInteger("problem","mode_m", 1);      // integer mode count across box
  const Real   phase    = pin->GetOrAddReal("problem","phase", 0.0);        // radians
  const Real   sigma_y  = pin->GetOrAddReal("problem","sigma_y", 0.02);     // envelope width (fraction of Ly) for velocity seed
  const Real   A_vy     = pin->GetOrAddReal("problem","A_vy", 1.0e-2);      // dimensionless vy amplitude (fraction of cs_ref)
  const Real   eta_y    = pin->GetOrAddReal("problem","eta_y", 0.0);        // interface displacement amplitude (absolute, same units as y)

  const Real gamma = peos->GetGamma();
  const Real gm1   = gamma - 1.0;
  if (NON_BAROTROPIC_EOS && !(gamma > 1.0)) {
    std::stringstream msg; msg << "KH requires gamma > 1 for adiabatic EOS.";
    ATHENA_ERROR(msg);
  }

  // Domain metrics
  const Real xmin = pmy_mesh->mesh_size.x1min, xmax = pmy_mesh->mesh_size.x1max;
  const Real ymin = pmy_mesh->mesh_size.x2min, ymax = pmy_mesh->mesh_size.x2max;
  const Real Lx = xmax - xmin, Ly = ymax - ymin;
  const Real ymid = 0.5*(ymin + ymax);
  const Real kx = 2.0*PI*static_cast<Real>(m_mode) / Lx;      // sinusoid wavenumber

  // Reference sound speed for velocity amplitude
  const Real rho_ref = 0.5*(rho_lo + rho_hi);
  const Real cs_ref  = std::sqrt(gamma * p0 / rho_ref);
  const Real sig_abs = std::max(1e-6, sigma_y * Ly);          // absolute width for envelope

  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      const Real y = pcoord->x2v(j);

      // If we use interface displacement, compute the local interface height:
      Real y_int = ymid;
      if (seed_kind == "interface") {
        // Sine-shaped interface: y_int(x) = ymid + eta_y * sin(kx*x + phase)
        // (eta_y may be 0 if user left default; then it reduces to a flat interface)
        y_int = ymid; // will set per-cell inside i-loop because it depends on x
      }

      for (int i = is; i <= ie; ++i) {
        const Real x = pcoord->x1v(i);

        // Pick layer by either (a) flat interface at ymid or (b) displaced interface
        bool lower;
        if (seed_kind == "interface" && std::abs(eta_y) > 0.0) {
          const Real y_interface = ymid + eta_y * std::sin(kx*(x - xmin) + phase);
          lower = (y < y_interface);
        } else {
          lower = (y < ymid);
        }

        const Real rho = lower ? rho_lo : rho_hi;
        const Real vx  = lower ? vx_lo  : vx_hi;
        Real vy = 0.0;

        // Velocity seed (vy) localized near the interface with Gaussian envelope
        if (seed_kind == "velocity") {
          const Real envelope = std::exp(-0.5 * SQR( (y - ymid) / sig_abs ));
          vy = A_vy * cs_ref * std::sin(kx * (x - xmin) + phase) * envelope;
        }

        // --- write conserved vars ---
        phydro->u(IDN, k, j, i) = rho;
        phydro->u(IM1, k, j, i) = rho * vx;   // shear along x
        phydro->u(IM2, k, j, i) = rho * vy;   // seed in y
        phydro->u(IM3, k, j, i) = 0.0;

        if (NON_BAROTROPIC_EOS) {
          Real Etot = p0/gm1 + 0.5*rho*(vx*vx + vy*vy);
          phydro->u(IEN, k, j, i) = Etot;
        }

#if NSCALARS > 0
        const int n = 0;
        const Real clr = lower ? color_lo : color_hi;
        pscalars->s(n, k, j, i) = rho * clr;
        pscalars->r(n, k, j, i) = clr;
#endif
      }
    }
  }
}
