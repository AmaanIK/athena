//========================================================================================
// Athena++ minimal problem generator: uniform rho, uniform pressure, zero velocity.
//========================================================================================

#include <cmath>
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

void Mesh::InitUserMeshData(ParameterInput *pin) { return; }

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  const Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  const Real p0   = pin->GetOrAddReal("problem", "p0",   1.0);
  const Real gamma = peos->GetGamma();
  const Real gm1   = gamma - 1.0;
  const Real xmin = pmy_mesh->mesh_size.x1min;
  const Real xmax = pmy_mesh->mesh_size.x1max;
  const Real Lx   = xmax - xmin;
  const Real ymin = pmy_mesh->mesh_size.x2min;
  const Real ymax = pmy_mesh->mesh_size.x2max;
  const Real Ly   = ymax - ymin;
  // Parameters for Gaussian density perturbation
  const Real x2c  = pin->GetOrAddReal("problem", "x2c", (ymin + ymax) * 0.5); // Center of perturbation
  const Real x1c  = pin->GetOrAddReal("problem", "x1c", (xmin + xmax) * 0.5); // Center of perturbation
  const Real A    = pin->GetOrAddReal("problem", "A",   1.0e-3);    // Amplitude of perturbation
  const Real sigma= pin->GetOrAddReal("problem", "sigma", 0.05*Lx); // Gaussian w/ width sigma
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        const Real dx = pcoord->x1v(i) - x1c;
        const Real dy = pcoord->x2v(j) - x2c;
        const Real r2 = dx*dx + dy*dy;
        const Real rho = rho0 ;
        phydro->u(IDN, k, j, i) = rho;
        phydro->u(IM1, k, j, i) = 0.0;
        phydro->u(IM2, k, j, i) = 0.0;
        phydro->u(IM3, k, j, i) = 0.0;

        if (NON_BAROTROPIC_EOS) {
          const Real p = p0 * std::pow(rho / rho0, gamma) * (1.0 + A * std::exp(-0.5 * (r2/(sigma*sigma)))); //Isentropic EOS
          phydro->u(IEN, k, j, i) = p/gm1; 
        }
      }
    }
  }
}
