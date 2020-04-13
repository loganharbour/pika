[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 10
    ny = 5
    xmax = 5
    ymax = 5
  []
[]

[Variables/phase]
  [./InitialCondition]
    type = FunctionIC
    function = '(x > 2.25) * 0.5 + (x > 2.5) * 0.5'
  []
[]

[AuxVariables]
  [refractive_index]
    [InitialCondition]
      type = FunctionIC
      function = '1.0 + (x > 2.5) * 0.3'
    []
  []
  [deposition]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[RayBCs/kill]
  type = KillRayBC
  boundary = 'top right bottom left'
[]

[RayKernels]
  [refraction]
    type = RefractionRayKernel
    phase = phase
    refractive_index = refractive_index
  []
  [deposition]
    type = OpticDepositionRayKernel
    variable = deposition
    depends_on = refraction
  []
[]

[UserObjects/lots]
  type = SimpleOpticRayStudy
  start_points = '0 2.5 0'
  directions = '1 0.5 0'
  data = '1 2 3'
  energy_groups = 3
  execute_on = 'initial'
[]

[Executioner]
  type = Steady
[]

[Problem]
  solve = false
[]

[Outputs]
  exodus = true
  csv = true
[]
