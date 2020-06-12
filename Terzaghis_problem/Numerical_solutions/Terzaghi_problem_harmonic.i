# Classical Terzaghi's problem. The analytical solution of this problem
# can be found on the github repository

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 20
  xmin = -1
  xmax = 1
  ymin = -1
  ymax = 1
  zmin = 0
  zmax = 100
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

[Functions]
  [./harmonic_stress]
    type = ParsedFunction
    value = '-1000 * sin(pi * t / 1500)^2'
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
    type = FunctionNeumannBC
    variable = disp_z
    function = harmonic_stress
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
    type = PorousFlowPorosityConst # only the initial value of this is used
    porosity = 0.2
  [../]
  [./biot_modulus]
    type = PorousFlowConstantBiotModulus
    biot_coefficient = 0.9
    fluid_bulk_modulus = 2.2E9
    solid_bulk_compliance = 1.0E-7
  [../]
  [./permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1.0E-11 0 0   0 1.0E-11 0   0 0 1.0E-11'
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
    point = '0 0 10'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p2]
    type = PointValue
    outputs = csv
    point = '0 0 20'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p3]
    type = PointValue
    outputs = csv
    point = '0 0 30'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p4]
    type = PointValue
    outputs = csv
    point = '0 0 40'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p5]
    type = PointValue
    outputs = csv
    point = '0 0 50'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p6]
    type = PointValue
    outputs = csv
    point = '0 0 60'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p7]
    type = PointValue
    outputs = csv
    point = '0 0 70'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p8]
    type = PointValue
    outputs = csv
    point = '0 0 80'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p9]
    type = PointValue
    outputs = csv
    point = '0 0 90'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./p99]
    type = PointValue
    outputs = csv
    point = '0 0 100'
    variable = porepressure
    use_displaced_mesh = false
  [../]
  [./dt]
    type = FunctionValuePostprocessor
    outputs = console
    function = if(t<0,0.001,50)
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
  start_time = -0.001
  end_time = 1000
  nl_rel_tol = 1E-3
  [./TimeStepper]
    type = PostprocessorDT
    postprocessor = dt
    dt = 0.001
  [../]
[]

[Outputs]
  execute_on = 'timestep_end'
  file_base = gold/Terzaghis_problem_harmonic
  [./csv]
    type = CSV
  [../]
[]
