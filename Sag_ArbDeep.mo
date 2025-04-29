model Sag_ArbDeep
  inner Modelica.Mechanics.MultiBody.World world annotation(Placement(visible = true, transformation(origin = {-130, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  // Parameters
  parameter Integer n = 10 "Number of springs";
  parameter Integer m = n - 1 "Number of bodies";
  parameter Modelica.Units.SI.Position[3] f1_r = {-2, 0, 0} "Left fixed point";
  parameter Modelica.Units.SI.Position[3] f2_r = {2, 0, 0} "Right fixed point";
  parameter Modelica.Units.SI.TranslationalSpringConstant c = 100 "Spring constant";
  parameter Modelica.Units.SI.Mass m_t = 20;
  final parameter Modelica.Units.SI.Length totalDistance = Modelica.Math.Vectors.length(f1_r - f2_r);
  final parameter Modelica.Units.SI.Length s_unstretched = ((totalDistance / n) * 1.1);
  final parameter Real s_list[n] = {s_unstretched for i in 1:n};
  final parameter Modelica.Units.SI.Mass m_s = m_t / m;
  final parameter Modelica.Units.SI.Acceleration g = world.g "Gravity magnitude";
  final parameter Modelica.Units.SI.Position body_r_0[m, 3](each fixed = false);
  //
  //
  //
  //final parameter Modelica.Units.SI.Position pos[m, 3] = Helpers.computeStaticPositions(m = m, L = totalDistance, c = c, m_s = m_s, g = g, s_unstretched = s_unstretched);
  //
  //
  final parameter Modelica.Units.SI.Position body_r_A[m, 3] = solveSpringMass(P1 = f1_r, P2 = f2_r, n = n, s0 = s_list, c = c, total_mass = m_t, g = g, max_sag = 10.0);
  //
  Modelica.Mechanics.MultiBody.Parts.Fixed fixed1(r = f1_r, shapeType = "sphere", length = 0.1, width = 0.1, height = 0.1, r_shape = {fixed1.r[1] + (fixed1.length / 2), fixed1.r[2], fixed1.r[3]}, lengthDirection = {-1, 0, 0}, color = {0, 180, 0}) annotation(Placement(visible = true, transformation(origin = {-70, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  //
  Modelica.Mechanics.MultiBody.Parts.Fixed fixed2(r = f2_r, shapeType = "sphere", length = 0.1, width = 0.1, height = 0.1, r_shape = {fixed2.r[1] - (fixed2.length / 2), fixed2.r[2], fixed2.r[3]}, lengthDirection = {1, 0, 0}, color = {0, 180, 0}) annotation(Placement(visible = true, transformation(origin = {70, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  //
  Modelica.Mechanics.MultiBody.Forces.Spring springs[n](each c = c, each s_unstretched = s_unstretched) annotation(Placement(visible = true, transformation(origin = {-10, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  //
  Modelica.Mechanics.MultiBody.Parts.Body bodies[m](each m = m_s, each r_0.fixed = false, r_0.start = body_r_0) annotation(Placement(visible = true, transformation(origin = {0, -30}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  //
  //

  // WSM searches in Resources/Library
  //
  function solveSpringMass
    input Real P1[3];
    input Real P2[3];
    input Integer n;
    input Real s0[n];
    input Real c;
    input Real total_mass;
    input Real g;
    input Real max_sag;
    output Real positions[n - 1, 3];
  //algorithm
    //for i in 1:n - 1 loop
    //  positions[i, 1] := s0[i];
    //  positions[i, 2] := s0[i];
    //  positions[i, 3] := s0[i];
    //end for;
  external "C"
    solve_spring_mass_c(P1, P2, n, s0, c, total_mass, g, max_sag, positions)
      annotation(
       Library = "spring_solver",
       LibraryDirectory = "modelica://Ropes/Resources/Library"
     );
  end solveSpringMass;
  //
  //
  //
initial equation
  //body_r_0 = pos_global;
  //body_r_0 = pos;
  for i in 1:m loop
    der(bodies[i].r_0) = zeros(3);
    // Zero initial velocity
    der(bodies[i].v_0) = zeros(3);
    // Zero initial acceleration
  end for;
  //
  //
equation
  // Connect fixed1 to first spring
  connect(fixed1.frame_b, springs[1].frame_a);
  // Connect springs and bodies
  for i in 1:m loop
    connect(springs[i].frame_b, bodies[i].frame_a);
    connect(bodies[i].frame_a, springs[i + 1].frame_a);
  end for;
  // Connect last spring to fixed2
  connect(springs[n].frame_b, fixed2.frame_b);
  //
  annotation(uses(Modelica(version = "4.0.0")), experiment(StopTime = 10.0), Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {0, 114, 195}, fillColor = {255, 255, 255}, extent = {{-100, -100}, {100, 100}}, radius = 25), Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name")}));
end Sag_ArbDeep;
