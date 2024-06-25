--- Import mesh

meshgen1 = mesh.MeshGenerator.Create
({
  inputs =
  {
    mesh.FromFileMeshGenerator.Create
    ({
      filename = "/home/staff/p/pablogarcia44/opensn/pincell/Pincell_nogap_0.msh"
    })
  }
})
mesh.MeshGenerator.Execute(meshgen1)
mesh.ExportToPVTU("pincell_mesh_only")


-- Create materials

mat_names = {"moderator","clad","fuel"}


materials = {}

Nmat = #mat_names 

for imat = 1, Nmat do
    materials[imat] = mat.AddMaterial(mat_names[imat])
    mat.SetProperty(materials[imat], TRANSPORT_XSECTIONS, OPENMC_XSLIB, "mgxs_10_pcm.h5", 294, mat_names[imat])
end

-- Setup Physics

-- Angular quadrature
pquad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,4, 32)
aquad.OptimizeForPolarSymmetry(pquad, 4.0*math.pi)

-- Solver
num_groups = 172

--############################################### Setup Physics
lbs_block = {
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, num_groups - 1 },
      angular_quadrature_handle = pquad,
      inner_linear_method = "krylov_richardson",
      l_max_its = 20,
      l_abs_tol = 1e-6,
      angle_aggregation_type = "polar",
    },
}

lbs_options = {
    boundary_conditions = {
        { name = "xmin", type = "reflecting" },
        { name = "xmax", type = "reflecting" },
        { name = "ymin", type = "reflecting" },
        { name = "ymax", type = "reflecting" },
        { name = "zmin", type = "reflecting" },
        { name = "zmax", type = "reflecting" },

      },
    scattering_order = 0,
  verbose_inner_iterations = true,
  verbose_outer_iterations = true,
  power_field_function_on = true,
  power_default_kappa = 1.0,
  power_normalization = 1.0,
  save_angular_flux = true,
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys1, lbs_options)

-- k_solver = lbs.PowerIterationKEigen.Create({
--   lbs_solver_handle = phys1,
--   k_tol = 1.0e-8,
-- })
k_solver = lbs.PowerIterationKEigenSMM.Create({
  lbs_solver_handle = phys1,
  accel_pi_verbose = true,
  k_tol = 1.0e-8,
  accel_pi_k_tol = 1.0e-8,
  accel_pi_max_its = 30,
  diff_sdm = "pwlc",
})
solver.Initialize(k_solver)
solver.Execute(k_solver)


fflist,count = lbs.GetScalarFieldFunctionList(phys1)
fieldfunc.ExportToVTKMulti(fflist,"block_mesh")
