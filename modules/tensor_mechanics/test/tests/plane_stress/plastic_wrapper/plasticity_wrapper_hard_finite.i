[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [square]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 1
    ny = 1
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]

[AuxVariables]
  [strain_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_xy]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_zz]
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
  [eff_plastic_strain]
    order = CONSTANT
    family = MONOMIAL
  []
  [vonmises]
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
  [stress_yy_el_0]
    type = ElementalVariableValue
    elementid = 0
    variable = stress_yy
  []
  [stress_xy_el_0]
    type = ElementalVariableValue
    elementid = 0
    variable = stress_xy
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
  [vonmises]
    type = ElementalVariableValue
    variable = vonmises
    elementid = 0
  []
  [eff_plastic_strain]
    type = ElementalVariableValue
    variable = eff_plastic_strain
    elementid = 0
  []

[]

[AuxKernels]
  [stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  []
  [stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  []
  [stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
  []
  [stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
  []
  [strain_xx]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_xx
    index_i = 0
    index_j = 0
  []
  [strain_xy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_xy
    index_i = 0
    index_j = 1
  []
  [strain_yy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_yy
    index_i = 1
    index_j = 1
  []
  [strain_zz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_zz
    index_i = 2
    index_j = 2
  []
  [vonmises]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    scalar_type = VonMisesStress
    execute_on = timestep_end
  []
  [eff_plastic_strain]
    type = MaterialRealAux
    property = effective_plastic_strain
    variable = eff_plastic_strain
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
    youngs_modulus = 40
  []
  [isotropic_plasticity]
    type = IsotropicPlasticityStressUpdate
    yield_stress = 4
    hardening_constant = 4
  []
  [radial_return_stress]
    type = ComputeMultipleInelasticStress
    #tangent_operator = elastic
    tangent_operator = nonlinear
    inelastic_models = 'isotropic_plasticity'
    perform_finite_strain_rotations = false
    compute = false
  []
  [strain]
    type = ComputePlaneFiniteStrain
    compute = false
  []
  [wps]
    type = PlaneStressWrapper
    stress_calculator = radial_return_stress
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
  l_tol = 1e-08

  # controls for nonlinear iterations
  nl_max_its = 20
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8

  # time control
  start_time = 0.0
  dt = 0.01
  dtmin = 0.01
  end_time = 1
[]

[Outputs]
  exodus = true
  csv = true
[]
