model Sag_Chat
  inner Modelica.Mechanics.MultiBody.World world annotation(
    Placement(visible = true, transformation(origin = {-130, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  // Parameters
  parameter Integer n = 10 "Number of springs";
  parameter Integer m = n - 1 "Number of bodies";
  parameter Modelica.Units.SI.Position[3] f1_r = {-2, 0, 0} "Left fixed point";
  parameter Modelica.Units.SI.Position[3] f2_r = {2, 0, 0} "Right fixed point";
  parameter Modelica.Units.SI.TranslationalSpringConstant c = 100 "Spring constant";
  parameter Modelica.Units.SI.Mass m_t = 20;
  
  // 🔄 New orientation-aware code
  final parameter Modelica.Units.SI.Position d_vec[3] = f2_r - f1_r;
  final parameter Modelica.Units.SI.Length L = Modelica.Math.Vectors.length(d_vec);
  final parameter Real e_x[3] = d_vec / L;
  final parameter Real e_z[3] = {0, 0, 1};
  final parameter Real e_y[3] = Modelica.Math.Vectors.normalize(Modelica.Math.Vectors.cross(e_z, e_x));
  final parameter Real e_z2[3] = Modelica.Math.Vectors.cross(e_x, e_y);

  final parameter Modelica.Units.SI.Length totalDistance = L;
  final parameter Modelica.Units.SI.Length s_unstretched = ((totalDistance/n)*1.1);
  final parameter Modelica.Units.SI.Mass m_s = m_t/m;
  final parameter Modelica.Units.SI.Acceleration g = 9.81 "Gravity magnitude";
  final parameter Modelica.Units.SI.Position body_r_0[m, 3](each fixed = false);

  // 💡 Update: compute local (x,y) and rotate into 3D
  final parameter Modelica.Units.SI.Position pos[m, 3] = 
    {f1_r + computeStaticPositions(m, L, c, m_s, g, s_unstretched)[i, 1]*e_x 
                  + computeStaticPositions(m, L, c, m_s, g, s_unstretched)[i, 2]*e_y 
                  for i in 1:m};

  Modelica.Mechanics.MultiBody.Parts.Fixed fixed1(r = f1_r, shapeType = "sphere", length = 0.1, width = 0.1, height = 0.1, r_shape = {fixed1.r[1] + (fixed1.length/2), fixed1.r[2], fixed1.r[3]}, lengthDirection = {-1, 0, 0}, color = {0, 180, 0}) annotation(
    Placement(visible = true, transformation(origin = {-70, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  Modelica.Mechanics.MultiBody.Parts.Fixed fixed2(r = f2_r, shapeType = "sphere", length = 0.1, width = 0.1, height = 0.1, r_shape = {fixed2.r[1] - (fixed2.length/2), fixed2.r[2], fixed2.r[3]}, lengthDirection = {1, 0, 0}, color = {0, 180, 0}) annotation(
    Placement(visible = true, transformation(origin = {70, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));

  Modelica.Mechanics.MultiBody.Forces.Spring springs[n](each c = c, each s_unstretched = s_unstretched) annotation(
    Placement(visible = true, transformation(origin = {-10, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  Modelica.Mechanics.MultiBody.Parts.Body bodies[m](each m = m_s,
                                                    each r_0.fixed = false,
                                                    r_0.start = body_r_0) annotation(
    Placement(visible = true, transformation(origin = {0, -30}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));

  function computeStaticPositions
    input Integer m "Number of bodies";
    input Modelica.Units.SI.Length L "Total horizontal span";
    input Modelica.Units.SI.TranslationalSpringConstant c "Spring stiffness";
    input Modelica.Units.SI.Mass m_s "Mass per body";
    input Modelica.Units.SI.Acceleration g "Gravity";
    input Modelica.Units.SI.Length s_unstretched "Unstretched spring length";
    input Integer maxIter = 100 "Maximum iterations";
    input Real tol = 1e-6 "Convergence tolerance";
    output Modelica.Units.SI.Position pos[m, 2] "Local [x, y] positions";
  protected
    Real x[m], y[m];
    Real F[2*m];
    Real F_perturbed[2*m];
    Real J[2*m, 2*m];
    Real delta[2*m];
    Real dx = L/(m + 1);
    Boolean converged;
    Real h = 1e-6;
    Real xL, yL, dxL, dyL, lenL, TL, FxL, FyL;
    Real xR, yR, dxR, dyR, lenR, TR, FxR, FyR;
    Integer i, k;
    Real x_temp[m];
    Real y_temp[m];
  algorithm
    for i in 1:m loop
      x[i] := -L/2 + dx*i;
      y[i] := -0.1;
    end for;
    converged := false;
    for iter in 1:maxIter loop
      for i in 1:m loop
        xL := if i == 1 then -L/2 else x[i - 1];
        yL := if i == 1 then 0 else y[i - 1];
        dxL := x[i] - xL;
        dyL := y[i] - yL;
        lenL := sqrt(dxL^2 + dyL^2);
        TL := c*(lenL - s_unstretched);
        FxL := TL*dxL/lenL;
        FyL := TL*dyL/lenL;

        xR := if i == m then L/2 else x[i + 1];
        yR := if i == m then 0 else y[i + 1];
        dxR := xR - x[i];
        dyR := yR - y[i];
        lenR := sqrt(dxR^2 + dyR^2);
        TR := c*(lenR - s_unstretched);
        FxR := TR*dxR/lenR;
        FyR := TR*dyR/lenR;

        F[2*i - 1] := FxR - FxL;
        F[2*i]     := FyR - FyL - m_s*g;
      end for;

      if Modelica.Math.Vectors.norm(F) < tol then
        converged := true;
        break;
      end if;

      for k in 1:2*m loop
        x_temp := x;
        y_temp := y;
        if mod(k, 2) == 1 then
          x_temp[div(k + 1, 2)] := x_temp[div(k + 1, 2)] + h;
        else
          y_temp[div(k, 2)] := y_temp[div(k, 2)] + h;
        end if;

        for i in 1:m loop
          xL := if i == 1 then -L/2 else x_temp[i - 1];
          yL := if i == 1 then 0 else y_temp[i - 1];
          dxL := x_temp[i] - xL;
          dyL := y_temp[i] - yL;
          lenL := sqrt(dxL^2 + dyL^2);
          TL := c*(lenL - s_unstretched);
          FxL := TL*dxL/lenL;
          FyL := TL*dyL/lenL;

          xR := if i == m then L/2 else x_temp[i + 1];
          yR := if i == m then 0 else y_temp[i + 1];
          dxR := xR - x_temp[i];
          dyR := yR - y_temp[i];
          lenR := sqrt(dxR^2 + dyR^2);
          TR := c*(lenR - s_unstretched);
          FxR := TR*dxR/lenR;
          FyR := TR*dyR/lenR;

          F_perturbed[2*i - 1] := FxR - FxL;
          F_perturbed[2*i] := FyR - FyL - m_s*g;
        end for;

        for r in 1:2*m loop
          J[r, k] := (F_perturbed[r] - F[r])/h;
        end for;
      end for;

      delta := Modelica.Math.Matrices.solve(J, -F);
      for i in 1:m loop
        x[i] := x[i] + delta[2*i - 1];
        y[i] := y[i] + delta[2*i];
      end for;
    end for;

    for i in 1:m loop
      pos[i, 1] := x[i];
      pos[i, 2] := y[i];
    end for;
  end computeStaticPositions;

  initial equation
    body_r_0 = pos;
    for i in 1:m loop
      der(bodies[i].r_0) = zeros(3);
      der(bodies[i].v_0) = zeros(3);
    end for;

  equation
    connect(fixed1.frame_b, springs[1].frame_a);
    for i in 1:m loop
      connect(springs[i].frame_b, bodies[i].frame_a);
      connect(bodies[i].frame_a, springs[i + 1].frame_a);
    end for;
    connect(springs[n].frame_b, fixed2.frame_b);

  annotation(
    uses(Modelica(version = "4.0.0")),
    experiment(StopTime = 10.0),
    Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})),
    Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}),
      graphics = {Rectangle(visible = true, lineColor = {0, 114, 195}, fillColor = {255, 255, 255}, extent = {{-100, -100}, {100, 100}}, radius = 25),
                  Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name")}));
end Sag_Chat;
