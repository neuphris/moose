[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [square]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 3
    ny = 3
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]

[AuxVariables]
  [nl_strain_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_zz]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Postprocessors]
  [react_z]
    type = MaterialTensorIntegral
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
  []
  [strain_zz]
    type = ElementalVariableValue
    elementid = 0
    variable = nl_strain_zz
  []
  [stress_xx_el_0]
    type = ElementalVariableValue
    elementid = 0
    variable = stress_xx
  []
  [stress_zz_el_0]
    type = ElementalVariableValue
    elementid = 0
    variable = stress_zz
  []
[]

[AuxKernels]

  [strain_zz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = nl_strain_zz
    index_i = 2
    index_j = 2
  []
  [stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  []
  [stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
  []
[]

[Functions]
  [pull]
    type = PiecewiseLinear
    x = '0     1   100'
    y = '0  0.2  0.00'
  []
[]

[BCs]
  [bottomx]
    type = DirichletBC
    boundary = left
    variable = disp_x
    value = 0.0
  []
  [bottomy]
    type = DirichletBC
    boundary = left
    variable = disp_y
    value = 0.0
  []
  [disp]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'right'
    function = pull
  []
[]

[Materials]
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    poissons_ratio = 0.3
    youngs_modulus = 1e6
  []
  [stress]
    type = ComputeFiniteStrainElasticStress
    #type = ComputeStrainIncrementBasedStress
    compute = false
  []
  [strain]
    type = ComputePlaneFiniteStrain
    compute = false
  []
  [wps]
    type = PlaneStressWrapper
    stress_calculator = stress
    strain_calculator = strain
  []
[]

[Kernels]
  [disp_x]
    type = StressDivergenceTensors
    variable = disp_x
    component = 0
  []
  [disp_y]
    type = StressDivergenceTensors
    variable = disp_y
    component = 1
  []
[]

[Executioner]
  type = Transient

  solve_type = PJFNK
  line_search = none

  # controls for linear iterations
  l_max_its = 100
  l_tol = 1e-06

  # controls for nonlinear iterations
  nl_max_its = 15
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-10

  # time control
  start_time = 0.0
  dt = 1
  dtmin = 1
  end_time = 1
[]

[Outputs]
  exodus = true
[]
