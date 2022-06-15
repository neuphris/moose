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
    xmax = 120.0
    ymin = 0.0
    ymax = 160.0
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
  [stress_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_xy]
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
  [dispx]
    type = PointValue
    point = '120 0 0'
    variable = disp_x
  []
  [dispy]
    type = PointValue
    point = '120 0 0'
    variable = disp_y
  []
  [stress_xx_el_0]
    type = ElementalVariableValue
    elementid = 0
    variable = stress_xx
  []
  [stress_yy_el_0]
    type = ElementalVariableValue
    elementid = 0
    variable = stress_yy
  []
  [stress_zz_el_0]
    type = ElementalVariableValue
    elementid = 0
    variable = stress_zz
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
  [strain_zz_el_0]
    type = ElementalVariableValue
    elementid = 0
    variable = nl_strain_zz
  []
[]

[AuxKernels]
  [strain_zz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = nl_strain_zz
    index_i = 2
    index_j = 2
    selected_qp = 3
  []
  [stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
    selected_qp = 3
  []
  [stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
    selected_qp = 3
  []
  [stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
    selected_qp = 3
  []
  [strain_xx]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_xx
    index_i = 0
    index_j = 0
    selected_qp = 3
  []
  [strain_yy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_yy
    index_i = 1
    index_j = 1
    selected_qp = 3
  []
  [stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
    selected_qp = 3
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

[]

[NodalKernels]
  [force_right]
    type = ConstantRate
    variable = disp_x
    boundary = right
    rate = 800
  []
[]

[Materials]
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    poissons_ratio = 0.25
    youngs_modulus = 1.08e6
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
  nl_max_its = 30
  nl_rel_tol = 1e-8
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
