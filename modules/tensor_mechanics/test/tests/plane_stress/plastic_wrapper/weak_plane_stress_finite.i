[GlobalParams]
  displacements = 'disp_x disp_y'
  out_of_plane_strain = strain_zz
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
  [strain_zz]
  []
[]

[AuxVariables]
  [nl_strain_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_zz0]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Postprocessors]
  [strain_zz_el_0]
    type = ElementalVariableValue
    elementid = 0
    variable = strain_zz
  []
  [stresszz_el_0]
    type = ElementalVariableValue
    elementid = 0
    variable = stress_zz0
  []
  # [stresszz_el_1]
  #   type = ElementalVariableValue
  #   elementid = 1
  #   variable = stress_zz0
  # []
  # [stresszz_el_2]
  #   type = ElementalVariableValue
  #   elementid = 2
  #   variable = stress_zz0
  # []
  # [stresszz_el_3]
  #   type = ElementalVariableValue
  #   elementid = 3
  #   variable = stress_zz0
  # []
[]

[AuxKernels]
  [strain_zz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = nl_strain_zz
    index_i = 2
    index_j = 2
  []
  [stress_zz0]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz0
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
  []
  [strain]
    type = ComputePlaneFiniteStrain
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
  [solid_z]
    type = WeakPlaneStress
    variable = strain_zz
  []
[]

[Executioner]
  type = Transient

  solve_type = PJFNK
  line_search = none

  # petsc_options_iname = '-ksp_type'
  # petsc_options_value = 'lu'

  # controls for linear iterations
  l_max_its = 100
  l_tol = 1e-06

  # controls for nonlinear iterations
  nl_max_its = 50
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12

  # time control
  start_time = 0.0
  dt = 1.0
  dtmin = 1.0
  end_time = 1.0
[]

[Outputs]
  exodus = true
[]
