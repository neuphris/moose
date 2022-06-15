[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [square]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 1
    ny = 1
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0
  []
  [allnodes]
    type = BoundingBoxNodeSetGenerator
    input = square
    bottom_left = '0.0 0.0 0.0'
    top_right = '1.0 1.0 0.0'
    new_boundary = 101
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
  [strain_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_yy]
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
  [stress_xx_el_0]
    type = ElementalVariableValue
    elementid = 0
    variable = stress_xx
  []
  [strain_xx_el_0]
    type = ElementalVariableValue
    elementid = 0
    variable = strain_xx
  []
  [strain_yy_el_0]
    type = ElementalVariableValue
    elementid = 0
    variable = strain_yy
  []
  [stress_zz_el_0]
    type = ElementalVariableValue
    elementid = 0
    variable = stress_zz
  []
  [strain_zz]
    type = ElementalVariableValue
    elementid = 0
    variable = nl_strain_zz
  []
[]

[AuxKernels]
  [nl_strain_zz3]
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
  [strain_xx]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_xx
    index_i = 0
    index_j = 0
  []
  [strain_yy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_yy
    index_i = 1
    index_j = 1
  []
  [stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
  []
[]

[BCs]
  [leftx]
    type = DirichletBC
    boundary = left
    variable = disp_x
    value = 0.0
  []
  [lefty]
    type = DirichletBC
    boundary = left
    variable = disp_y
    value = 0.0
  []
  [disp]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'right'
    function = displacement
  []
[]

[Functions]
  [displacement]
    type = PiecewiseLinear
    x = '0.0 1.0'
    y = '0.0 0.5'
  []
[]

[Materials]
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    poissons_ratio = 0.3
    youngs_modulus = 4
  []
  [stress]
    type = ComputeLinearElasticStress
    compute = false
  []
  [strain]
    type = ComputePlaneSmallStrain
    compute = false
  []
  [wps]
    type = PlaneStressWrapper
    stress_calculator = stress
    strain_calculator = strain
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  # controls for nonlinear iterations
  nl_max_its = 15
  nl_rel_tol = 1e-06
  # time control
  start_time = 0.0
  dt = 1
  dtmin = 1
  end_time = 1.0
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

[Outputs]
  exodus = true
[]
