# Terzaghis problem in a multy layers system. The analytical solution of this problem
# can be found on the github repository

[Mesh]
  [the_mesh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 30
    xmin = -1
    xmax = 1
    ymin = -1
    ymax = 1
    zmin = -20
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
    [./aquifer_03]
      type = SubdomainBoundingBoxGenerator
      block_id = 3
      bottom_left = '0 0 -20'
      top_right = '0 0 -10'
      input = aquifer_02
    [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  PorousFlowDictator = dictator
  biot_coefficient = 0.9
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
    type = FunctionDirichletBC
    variable = porepressure
    function =  0
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
  gravity = '0 0 0'
  fp = the_simple_fluid
[]


[Materials]
  [./elasticity_tensor_01]
    type = ComputeIsotropicElasticityTensor
    bulk_modulus = 1.0E7
    shear_modulus = 6.0E7
    block = 1
  [../]
  [./elasticity_tensor_02]
    type = ComputeIsotropicElasticityTensor
    bulk_modulus = 2.5E7
    shear_modulus = 7.0E7
    block = 2
  [../]
  [./elasticity_tensor_03]
    type = ComputeIsotropicElasticityTensor
    bulk_modulus = 5.0E7
    shear_modulus = 8.0E7
    block = 3
  [../]
  [./strain]
    type = ComputeSmallStrain
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [./porosity_01]
    type = PorousFlowPorosityConst
    porosity = 0.3
    block = 1
  [../]
  [./porosity_02]
    type = PorousFlowPorosityConst
    porosity = 0.25
    block = 2
  [../]
  [./porosity_03]
    type = PorousFlowPorosityConst
    porosity = 0.2
    block = 3
  [../]
  [./biot_modulus_01]
    type = PorousFlowConstantBiotModulus
    biot_coefficient = 0.9
    fluid_bulk_modulus = 2.2E9
    solid_bulk_compliance = 1.0E-7
    block = 1
  [../]
  [./biot_modulus_02]
    type = PorousFlowConstantBiotModulus
    biot_coefficient = 0.9
    fluid_bulk_modulus = 2.2E9
    solid_bulk_compliance = 4.0E-8
    block = 2
  [../]
  [./biot_modulus_03]
    type = PorousFlowConstantBiotModulus
    biot_coefficient = 0.9
    fluid_bulk_modulus = 2.2E9
    solid_bulk_compliance = 2.0E-8
    block = 3
  [../]
  [./permeability_01]
    type = PorousFlowPermeabilityConst
    permeability = '1.0E-10 0 0   0 1.0E-10 0   0 0 1.0E-10'
    block = 1
  [../]
  [./permeability_02]
    type = PorousFlowPermeabilityConst
    permeability = '1.0E-11 0 0   0 1.0E-11 0   0 0 1.0E-11'
    block = 2
  [../]
  [./permeability_03]
    type = PorousFlowPermeabilityConst
    permeability = '1.0E-13 0 0   0 1.0E-13 0   0 0 1.0E-13'
    block = 3
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
    point = '0 0 -20'
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
    point = '0 0 -7'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./-p4]
    type = PointValue
    outputs = csv
    point = '0 0 -6'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./-p5]
    type = PointValue
    outputs = csv
    point = '0 0 -5'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./-p6]
    type = PointValue
    outputs = csv
    point = '0 0 -4'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./-p7]
    type = PointValue
    outputs = csv
    point = '0 0 -3'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./-p8]
    type = PointValue
    outputs = csv
    point = '0 0 -2'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./-p9]
    type = PointValue
    outputs = csv
    point = '0 0 -1'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./-p10]
    type = PointValue
    outputs = csv
    point = '0 0 -19'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./-p11]
  type = PointValue
  outputs = csv
  point = '0 0 -18'
  variable = porepressure
  use_displaced_mesh = false
[../]
[./-p12]
  type = PointValue
  outputs = csv
  point = '0 0 -17'
  variable = porepressure
  use_displaced_mesh = false
[../]
[./-p13]
  type = PointValue
  outputs = csv
  point = '0 0 -16'
  variable = porepressure
  use_displaced_mesh = false
[../]
[./-p14]
  type = PointValue
  outputs = csv
  point = '0 0 -15'
  variable = porepressure
  use_displaced_mesh = false
[../]
[./-p15]
  type = PointValue
  outputs = csv
  point = '0 0 -14'
  variable = porepressure
  use_displaced_mesh = false
[../]
[./-p16]
  type = PointValue
  outputs = csv
  point = '0 0 -13'
  variable = porepressure
  use_displaced_mesh = false
[../]
[./-p17]
  type = PointValue
  outputs = csv
  point = '0 0 -12'
  variable = porepressure
  use_displaced_mesh = false
[../]
[./-p18]
  type = PointValue
  outputs = csv
  point = '0 0 -11'
  variable = porepressure
  use_displaced_mesh = false
[../]
[./-p19]
  type = PointValue
  outputs = csv
  point = '0 0 -10'
  variable = porepressure
  use_displaced_mesh = false
[../]
[./-p20]
  type = PointValue
  outputs = csv
  point = '0 0 -8'
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
  [./andy]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  start_time = -0.0001
  end_time = 1000
  nl_rel_tol = 1E-6
  [./TimeStepper]
    type = PostprocessorDT
    postprocessor = dt
    dt = 0.0001
  [../]
[]

[Outputs]
  execute_on = 'timestep_end'
  file_base = gold/Terzaghis_problem_multy_layers
  [./csv]
    type = CSV
  [../]
[]
