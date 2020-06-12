# Terzaghis problem in a two layers system. The analytical solution of this problem
# can be found on the github repository

[Mesh]
  [the_mesh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 20
    xmin = -1
    xmax = 1
    ymin = -1
    ymax = 1
    zmin = -10
    zmax = 10
  [../]
    [./aquifer_01]
      type = SubdomainBoundingBoxGenerator
      block_id = 1
      bottom_left = '0 0 0'
      top_right = '0 0 10'
      input = the_mesh
    [../]
    [./aquifer_02]
      type = SubdomainBoundingBoxGenerator
      block_id = 2
      bottom_left = '0 0 -10'
      top_right = '0 0 0'
      input = aquifer_01
    [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  PorousFlowDictator = dictator
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
  [./porepressure]
    scaling = 1E11
  [../]
[]

[BCs]
  [./confinex]
    type = PresetBC
    variable = disp_x
    value = 0
    boundary = 'left right'
  [../]
  [./confiney]
    type = PresetBC
    variable = disp_y
    value = 0
    boundary = 'bottom top'
  [../]
  [./basefixed]
    type = PresetBC
    variable = disp_z
    value = 0
    boundary = back
  [../]
  [./topdrained]
    type = DirichletBC
    variable = porepressure
    value = 0
    boundary = front
  [../]
  [./topload]
    type = NeumannBC
    variable = disp_z
    value = -1000
    boundary = front
  [../]
[]

[Modules]
  [./FluidProperties]
    [./the_simple_fluid]
      type = SimpleFluidProperties
      thermal_expansion = 0.0
      bulk_modulus = 2.2E9
      viscosity = 1E-3
      density0 = 1000.0
    [../]
  [../]
[]

[PorousFlowBasicTHM]
  coupling_type = HydroMechanical
  displacements = 'disp_x disp_y disp_z'
  multiply_by_density = false
  porepressure = porepressure
  biot_coefficient = 0.9
  gravity = '0 0 0'
  fp = the_simple_fluid
[]


[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    bulk_modulus = 1.0E7
    shear_modulus = 6.9E7
  [../]
  [./strain]
    type = ComputeSmallStrain
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [./porosity]
    type = PorousFlowPorosityConst
    porosity = 0.2
  [../]
  [./biot_modulus]
    type = PorousFlowConstantBiotModulus
    biot_coefficient = 0.9
    fluid_bulk_modulus = 2.2E9
    solid_bulk_compliance = 1.0E-7
  [../]
  [./permeability_01]
    type = PorousFlowPermeabilityConst
    permeability = '1.0E-10 0 0   0 1.0E-10 0   0 0 1.0E-10'
    block = 1
  [../]
  [./permeability_02]
    type = PorousFlowPermeabilityConst
    permeability = '1.0E-12 0 0   0 1.0E-12 0   0 0 1.0E-12'
    block = 2
  [../]
[]

[Postprocessors]
  [./p0]
    type = PointValue
    outputs = csv
    point = '0 0 0'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p1]
    type = PointValue
    outputs = csv
    point = '0 0 1'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p2]
    type = PointValue
    outputs = csv
    point = '0 0 2'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p3]
    type = PointValue
    outputs = csv
    point = '0 0 3'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p4]
    type = PointValue
    outputs = csv
    point = '0 0 4'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p5]
    type = PointValue
    outputs = csv
    point = '0 0 5'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p6]
    type = PointValue
    outputs = csv
    point = '0 0 6'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p7]
    type = PointValue
    outputs = csv
    point = '0 0 7'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p8]
    type = PointValue
    outputs = csv
    point = '0 0 8'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p9]
    type = PointValue
    outputs = csv
    point = '0 0 9'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p99]
    type = PointValue
    outputs = csv
    point = '0 0 10'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./-p1]
    type = PointValue
    outputs = csv
    point = '0 0 -10'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./-p2]
    type = PointValue
    outputs = csv
    point = '0 0 -9'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./-p3]
    type = PointValue
    outputs = csv
    point = '0 0 -8'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./-p4]
    type = PointValue
    outputs = csv
    point = '0 0 -7'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./-p5]
    type = PointValue
    outputs = csv
    point = '0 0 -6'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./-p6]
    type = PointValue
    outputs = csv
    point = '0 0 -5'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./-p7]
    type = PointValue
    outputs = csv
    point = '0 0 -4'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./-p8]
    type = PointValue
    outputs = csv
    point = '0 0 -3'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./-p9]
    type = PointValue
    outputs = csv
    point = '0 0 -2'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./-p99]
    type = PointValue
    outputs = csv
    point = '0 0 -1'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./dt]
    type = FunctionValuePostprocessor
    outputs = console
    function = if(t<0,0.0001,10)
  [../]
[]

[Preconditioning]
  [./mumps_is_best_for_parallel_jobs]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  start_time = -0.0001
  end_time = 1000
  nl_rel_tol = 1E-3
  [./TimeStepper]
    type = PostprocessorDT
    postprocessor = dt
    dt = 0.0001
  [../]
[]

[Outputs]
  execute_on = 'timestep_end'
  file_base = gold/Terzaghis_problem_2_layers
  [./csv]
    type = CSV
  [../]
[]
