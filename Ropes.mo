package Ropes
  package Interfaces
    connector CableFrame
      SI.Position r_0[3];
      Real e_pre[3] "The normalized direction
of the last cable element in world
coordinates";
      Real e_next[3] "The normalized direction
of the next cable element in world
coordinates";
      flow SI.Force f[3] "Cut-force resolved in WORLD frame.";
      annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, origin = {60, 16.667}, textColor = {64, 64, 64}, extent = {{-160, 33.333}, {40, 73.333}}, textString = "%name")}));
    end CableFrame;

    connector CF_a
      extends CableFrame;
      annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Ellipse(visible = true, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-40, -20}, {40, 20}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Ellipse(visible = true, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-40, -20}, {40, 20}})}));
    end CF_a;

    connector CF_b
      extends CableFrame;
      annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Ellipse(visible = true, lineColor = {0, 0, 255}, fillColor = {0, 0, 255}, lineThickness = 5, extent = {{-40, -20}, {40, 20}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Ellipse(visible = true, lineColor = {0, 0, 255}, fillColor = {0, 0, 255}, lineThickness = 5, extent = {{-40, -20}, {40, 20}})}));
    end CF_b;
    annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {0.248, 0.044}, lineColor = {56, 56, 56}, fillColor = {128, 202, 255}, fillPattern = FillPattern.Solid, points = {{99.752, 100}, {99.752, 59.956}, {99.752, -50}, {100, -100}, {49.752, -100}, {-19.752, -100.044}, {-100.248, -100}, {-100.248, -50}, {-90.248, 29.956}, {-90.248, 79.956}, {-40.248, 79.956}, {-20.138, 79.813}, {-0.248, 79.956}, {19.752, 99.956}, {39.752, 99.956}, {59.752, 99.956}}, smooth = Smooth.Bezier), Polygon(visible = true, origin = {0, -13.079}, lineColor = {192, 192, 192}, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.HorizontalCylinder, points = {{100, -86.921}, {50, -86.921}, {-50, -86.921}, {-100, -86.921}, {-100, -36.921}, {-100, 53.079}, {-100, 103.079}, {-50, 103.079}, {0, 103.079}, {20, 83.079}, {50, 83.079}, {100, 83.079}, {100, 33.079}, {100, -36.921}}, smooth = Smooth.Bezier), Polygon(visible = true, origin = {0, -10.704}, lineColor = {113, 113, 113}, fillColor = {255, 255, 255}, points = {{100, -89.296}, {50, -89.296}, {-50, -89.296}, {-100, -89.296}, {-100, -39.296}, {-100, 50.704}, {-100, 100.704}, {-50, 100.704}, {0, 100.704}, {20, 80.704}, {50, 80.704}, {100, 80.704}, {100, 30.704}, {100, -39.296}}, smooth = Smooth.Bezier), Rectangle(visible = true, origin = {-30, 20}, lineColor = {56, 56, 56}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid, extent = {{-22.5, -22.5}, {22.5, 22.5}}), Ellipse(visible = true, origin = {35, 20}, lineColor = {56, 56, 56}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid, extent = {{-22.5, -22.5}, {22.5, 22.5}}), Polygon(visible = true, origin = {35.018, -44.722}, lineColor = {56, 56, 56}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid, points = {{-15.018, -5.278}, {-22.518, 2.222}, {-0.018, 22.222}, {22.482, 2.222}, {14.982, -5.278}, {4.982, 4.722}, {4.982, -22.778}, {-4.937, -22.778}, {-4.937, 4.722}}), Polygon(visible = true, origin = {-40, -54.167}, lineColor = {56, 56, 56}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid, points = {{-12.5, 31.667}, {32.5, 9.167}, {-12.5, -15.833}})}));
  end Interfaces;

  package Helpers
    function computeStaticPositionsArbitrary
      input Integer m "Number of bodies";
      input Modelica.Units.SI.Position[3] f1_r "Left fixed point";
      input Modelica.Units.SI.Position[3] f2_r "Right fixed point";
      input Modelica.Units.SI.TranslationalSpringConstant c "Spring stiffness";
      input Modelica.Units.SI.Mass m_s "Mass per body";
      input Modelica.Units.SI.Acceleration g_v[3] "Gravity";
      input Modelica.Units.SI.Length s_unstretched "Unstretched spring length";
      input Integer maxIter = 100 "Maximum iterations";
      input Real tol = 1e-6 "Convergence tolerance";
      output Modelica.Units.SI.Position pos[m, 3] "Equilibrium positions";
    protected
      Real f1_to_f2[3] = f2_r - f1_r;
      Real L = Modelica.Math.Vectors.norm(f1_to_f2);
      Real r[m, 3] "Position vectors of the bodies";
      Real r_perturbed[m, 3];
      Real F[3 * m] "Residual vector";
      Real F_perturbed[3 * m] "Perturbed residual vector";
      Real J[3 * m, 3 * m] "Jacobian matrix";
      Real delta[3 * m] "Newton-Raphson update";
      Real h = 1e-6 "Perturbation step size";
      Boolean converged = false;
      Integer i, k, body, coord;
      Real alpha, perturbation, g_norm;
      Real[3] g_dir;
      Real xL, yL, zL, dxL, dyL, dzL, lenL, TL, FxL, FyL, FzL;
      Real xR, yR, zR, dxR, dyR, dzR, lenR, TR, FxR, FyR, FzR;
    algorithm
      // Initial guess: positions along the line between f1_r and f2_r with a small perturbation in the direction of gravity
      g_norm := Modelica.Math.Vectors.norm(g_v);
      for i in 1:m loop
        alpha := i / (m + 1);
        for j in 1:3 loop
          r[i, j] := f1_r[j] + alpha * f1_to_f2[j];
        end for;
        if g_norm > 0 then
          g_dir := g_v / g_norm;
          perturbation := 0.1 * s_unstretched;
          for j in 1:3 loop
            r[i, j] := r[i, j] + perturbation * g_dir[j];
          end for;
        end if;
      end for;
      // Newton-Raphson iteration
      for iter in 1:maxIter loop
        // Compute residuals
        for i in 1:m loop
          // Left spring
          if i == 1 then
            xL := f1_r[1];
            yL := f1_r[2];
            zL := f1_r[3];
          else
            xL := r[i - 1, 1];
            yL := r[i - 1, 2];
            zL := r[i - 1, 3];
          end if;
          dxL := r[i, 1] - xL;
          dyL := r[i, 2] - yL;
          dzL := r[i, 3] - zL;
          lenL := sqrt(dxL ^ 2 + dyL ^ 2 + dzL ^ 2);
          TL := c * (lenL - s_unstretched);
          FxL := TL * dxL / lenL;
          FyL := TL * dyL / lenL;
          FzL := TL * dzL / lenL;
          // Right spring
          if i == m then
            xR := f2_r[1];
            yR := f2_r[2];
            zR := f2_r[3];
          else
            xR := r[i + 1, 1];
            yR := r[i + 1, 2];
            zR := r[i + 1, 3];
          end if;
          dxR := xR - r[i, 1];
          dyR := yR - r[i, 2];
          dzR := zR - r[i, 3];
          lenR := sqrt(dxR ^ 2 + dyR ^ 2 + dzR ^ 2);
          TR := c * (lenR - s_unstretched);
          FxR := TR * dxR / lenR;
          FyR := TR * dyR / lenR;
          FzR := TR * dzR / lenR;
          // Residuals (F = spring forces + gravity)
          F[3 * i - 2] := FxR - FxL + m_s * g_v[1];
          F[3 * i - 1] := FyR - FyL + m_s * g_v[2];
          F[3 * i] := FzR - FzL + m_s * g_v[3];
        end for;
        // Check convergence
        if Modelica.Math.Vectors.norm(F) < tol then
          converged := true;
          break;
        end if;
        // Compute numerical Jacobian
        for k in 1:3 * m loop
          body := div(k - 1, 3) + 1;
          coord := mod(k - 1, 3) + 1;
          // Perturb coordinate
          r_perturbed := r;
          r_perturbed[body, coord] := r_perturbed[body, coord] + h;
          // Compute perturbed residuals
          for i in 1:m loop
            // Left spring
            if i == 1 then
              xL := f1_r[1];
              yL := f1_r[2];
              zL := f1_r[3];
            else
              xL := r_perturbed[i - 1, 1];
              yL := r_perturbed[i - 1, 2];
              zL := r_perturbed[i - 1, 3];
            end if;
            dxL := r_perturbed[i, 1] - xL;
            dyL := r_perturbed[i, 2] - yL;
            dzL := r_perturbed[i, 3] - zL;
            lenL := sqrt(dxL ^ 2 + dyL ^ 2 + dzL ^ 2);
            TL := c * (lenL - s_unstretched);
            FxL := TL * dxL / lenL;
            FyL := TL * dyL / lenL;
            FzL := TL * dzL / lenL;
            // Right spring
            if i == m then
              xR := f2_r[1];
              yR := f2_r[2];
              zR := f2_r[3];
            else
              xR := r_perturbed[i + 1, 1];
              yR := r_perturbed[i + 1, 2];
              zR := r_perturbed[i + 1, 3];
            end if;
            dxR := xR - r_perturbed[i, 1];
            dyR := yR - r_perturbed[i, 2];
            dzR := zR - r_perturbed[i, 3];
            lenR := sqrt(dxR ^ 2 + dyR ^ 2 + dzR ^ 2);
            TR := c * (lenR - s_unstretched);
            FxR := TR * dxR / lenR;
            FyR := TR * dyR / lenR;
            FzR := TR * dzR / lenR;
            // Perturbed residuals
            F_perturbed[3 * i - 2] := FxR - FxL + m_s * g_v[1];
            F_perturbed[3 * i - 1] := FyR - FyL + m_s * g_v[2];
            F_perturbed[3 * i] := FzR - FzL + m_s * g_v[3];
          end for;
          // Fill Jacobian column
          for j in 1:3 * m loop
            J[j, k] := (F_perturbed[j] - F[j]) / h;
          end for;
        end for;
        // Solve for update
        delta := Modelica.Math.Matrices.solve(J, -F);
        // Update positions
        for i in 1:m loop
          r[i, 1] := r[i, 1] + delta[3 * i - 2];
          r[i, 2] := r[i, 2] + delta[3 * i - 1];
          r[i, 3] := r[i, 3] + delta[3 * i];
        end for;
      end for;
      // Assign output
      for i in 1:m loop
        pos[i, :] := r[i, :];
      end for;
    end computeStaticPositionsArbitrary;

    function computeStaticPositions
      input Integer m "Number of bodies";
      input Modelica.Units.SI.Length L "Total horizontal span";
      input Modelica.Units.SI.TranslationalSpringConstant c "Spring stiffness";
      input Modelica.Units.SI.Mass m_s "Mass per body";
      input Modelica.Units.SI.Acceleration g "Gravity";
      input Modelica.Units.SI.Length s_unstretched "Unstretched spring length";
      input Integer maxIter = 100 "Maximum iterations";
      input Real tol = 1e-6 "Convergence tolerance";
      output Modelica.Units.SI.Position pos[m, 3] "Equilibrium [x, y] positions";
    protected
      Real x[m], y[m];
      Real F[2 * m];
      Real F_perturbed[2 * m];
      Real J[2 * m, 2 * m];
      Real delta[2 * m];
      Real dx = L / (m + 1);
      Boolean converged;
      Real h = 1e-6;
      Real L_left, L_right;
      Real T_left, T_right;
      Real Fx_left, Fx_right;
      Real Fy_left, Fy_right;
      Real xL, yL, dxL, dyL, lenL, TL, FxL, FyL;
      Real xR, yR, dxR, dyR, lenR, TR, FxR, FyR;
      Integer i, k;
      Real x_temp[m];
      Real y_temp[m];
    algorithm
      // Initial guess
      for i in 1:m loop
        x[i] := -L / 2 + dx * i;
        y[i] := -0.1;
        // some initial sag
      end for;
      converged := false;
      for iter in 1:maxIter loop
        // Compute residuals
        for i in 1:m loop
          // Left spring
          xL := if i == 1 then -L / 2 else x[i - 1];
          yL := if i == 1 then 0 else y[i - 1];
          dxL := x[i] - xL;
          dyL := y[i] - yL;
          lenL := sqrt(dxL ^ 2 + dyL ^ 2);
          TL := c * (lenL - s_unstretched);
          FxL := TL * dxL / lenL;
          FyL := TL * dyL / lenL;
          // Right spring
          xR := if i == m then L / 2 else x[i + 1];
          yR := if i == m then 0 else y[i + 1];
          dxR := xR - x[i];
          dyR := yR - y[i];
          lenR := sqrt(dxR ^ 2 + dyR ^ 2);
          TR := c * (lenR - s_unstretched);
          FxR := TR * dxR / lenR;
          FyR := TR * dyR / lenR;
          // Residuals
          F[2 * i - 1] := FxR - FxL;
          F[2 * i] := FyR - FyL - m_s * g;
        end for;
        // Check convergence
        if Modelica.Math.Vectors.norm(F) < tol then
          converged := true;
          break;
        end if;
        // Numerical Jacobian
        for k in 1:2 * m loop
          // Copy original x/y
          x_temp := x;
          y_temp := y;
          if mod(k, 2) == 1 then
            x_temp[div(k + 1, 2)] := x_temp[div(k + 1, 2)] + h;
          else
            y_temp[div(k, 2)] := y_temp[div(k, 2)] + h;
          end if;
          // Recompute F_perturbed
          for i in 1:m loop
            xL := if i == 1 then -L / 2 else x_temp[i - 1];
            yL := if i == 1 then 0 else y_temp[i - 1];
            dxL := x_temp[i] - xL;
            dyL := y_temp[i] - yL;
            lenL := sqrt(dxL ^ 2 + dyL ^ 2);
            TL := c * (lenL - s_unstretched);
            FxL := TL * dxL / lenL;
            FyL := TL * dyL / lenL;
            xR := if i == m then L / 2 else x_temp[i + 1];
            yR := if i == m then 0 else y_temp[i + 1];
            dxR := xR - x_temp[i];
            dyR := yR - y_temp[i];
            lenR := sqrt(dxR ^ 2 + dyR ^ 2);
            TR := c * (lenR - s_unstretched);
            FxR := TR * dxR / lenR;
            FyR := TR * dyR / lenR;
            F_perturbed[2 * i - 1] := FxR - FxL;
            F_perturbed[2 * i] := FyR - FyL - m_s * g;
          end for;
          // Fill Jacobian column
          for r in 1:2 * m loop
            J[r, k] := (F_perturbed[r] - F[r]) / h;
          end for;
        end for;
        // Solve for update
        delta := Modelica.Math.Matrices.solve(J, -F);
        // Apply update
        for i in 1:m loop
          x[i] := x[i] + delta[2 * i - 1];
          y[i] := y[i] + delta[2 * i];
        end for;
      end for;
      // Final output
      for i in 1:m loop
        pos[i, 1] := x[i];
        pos[i, 2] := y[i];
        pos[i, 3] := 0;
      end for;
    end computeStaticPositions;
    annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {0.248, 0.044}, lineColor = {56, 56, 56}, fillColor = {128, 202, 255}, fillPattern = FillPattern.Solid, points = {{99.752, 100}, {99.752, 59.956}, {99.752, -50}, {100, -100}, {49.752, -100}, {-19.752, -100.044}, {-100.248, -100}, {-100.248, -50}, {-90.248, 29.956}, {-90.248, 79.956}, {-40.248, 79.956}, {-20.138, 79.813}, {-0.248, 79.956}, {19.752, 99.956}, {39.752, 99.956}, {59.752, 99.956}}, smooth = Smooth.Bezier), Polygon(visible = true, origin = {0, -13.079}, lineColor = {192, 192, 192}, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.HorizontalCylinder, points = {{100, -86.921}, {50, -86.921}, {-50, -86.921}, {-100, -86.921}, {-100, -36.921}, {-100, 53.079}, {-100, 103.079}, {-50, 103.079}, {0, 103.079}, {20, 83.079}, {50, 83.079}, {100, 83.079}, {100, 33.079}, {100, -36.921}}, smooth = Smooth.Bezier), Polygon(visible = true, origin = {0, -10.704}, lineColor = {113, 113, 113}, fillColor = {255, 255, 255}, points = {{100, -89.296}, {50, -89.296}, {-50, -89.296}, {-100, -89.296}, {-100, -39.296}, {-100, 50.704}, {-100, 100.704}, {-50, 100.704}, {0, 100.704}, {20, 80.704}, {50, 80.704}, {100, 80.704}, {100, 30.704}, {100, -39.296}}, smooth = Smooth.Bezier), Ellipse(visible = true, origin = {0, -15}, lineColor = {56, 56, 56}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid, extent = {{-60, -60}, {60, 60}}), Ellipse(visible = true, origin = {0.257, -15.258}, lineColor = {56, 56, 56}, fillColor = {244, 244, 244}, fillPattern = FillPattern.Solid, extent = {{-45.257, -45.257}, {45.257, 45.257}}), Rectangle(visible = true, origin = {-0.219, -15.813}, rotation = 45, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid, extent = {{-46.497, -8.082}, {46.497, 8.082}}), Line(visible = true, origin = {6.058, -20.995}, points = {{31.492, 31.535}, {-31.492, -31.535}}), Line(visible = true, origin = {-5.572, -9.534}, points = {{31.492, 31.535}, {-31.492, -31.535}})}));
  end Helpers;

  model Sag_Arb
    inner Modelica.Mechanics.MultiBody.World world annotation(Placement(visible = true, transformation(origin = {-130, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    // Parameters
    parameter Integer n = 10 "Number of springs";
    parameter Integer m = n - 1 "Number of bodies";
    parameter Modelica.Units.SI.Position[3] f1_r = {15, 0, -7} "Left fixed point";
    parameter Modelica.Units.SI.Position[3] f2_r = {-15, 0, 7} "Right fixed point";
    parameter Modelica.Units.SI.TranslationalSpringConstant c = 100 "Spring constant";
    parameter Modelica.Units.SI.Mass m_t = 20;
    final parameter Modelica.Units.SI.Length totalDistance = Modelica.Math.Vectors.length(f1_r - f2_r);
    final parameter Modelica.Units.SI.Length s_unstretched = ((totalDistance / n) * 1.1);
    final parameter Modelica.Units.SI.Mass m_s = m_t / m;
    final parameter Modelica.Units.SI.Acceleration g = world.g "Gravity magnitude";
    final parameter Real g_world[3] = world.g * world.n;
    final parameter Modelica.Units.SI.Position body_r_0[m, 3](each fixed = false);
    Modelica.Units.SI.Position body_r_A[m, 3];
    //
    //
    // Setup local frame
    final parameter Modelica.Units.SI.Position d_vec[3] = f2_r - f1_r;
    final parameter Modelica.Units.SI.Length L = Modelica.Math.Vectors.length(d_vec);
    final parameter Modelica.Units.SI.Position center[3] = f1_r + 0.5 * d_vec;
    final parameter Real e_x[3] = d_vec / L;
    final parameter Real e_z_temp[3] = Modelica.Math.Vectors.normalize(cross(world.n, e_x));
    final parameter Real e_y[3] = Modelica.Math.Vectors.normalize(cross(e_z_temp, e_x));
    final parameter Real e_z[3] = cross(e_x, e_y);
    final parameter Real R_local[3, 3] = [e_x, e_y, e_z];
    // columns: local x, y, z
    //
    final parameter Modelica.Units.SI.Position pos[m, 3] = Helpers.computeStaticPositions(m = m, L = totalDistance, c = c, m_s = m_s, g = g, s_unstretched = s_unstretched);
    //
    //
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
  initial equation
    body_r_0 = body_r_A;
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
  end Sag_Arb;

  model Sag_1
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
    final parameter Modelica.Units.SI.Mass m_s = m_t / m;
    final parameter Modelica.Units.SI.Acceleration g = world.g "Gravity magnitude";
    final parameter Modelica.Units.SI.Position body_r_0[m, 3](each fixed = false);
    final parameter Modelica.Units.SI.Position pos[m, 3] = Helpers.computeStaticPositions(m = m, L = totalDistance, c = c, m_s = m_s, g = g, s_unstretched = s_unstretched);
    //
    //
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
    //
  initial equation
    body_r_0 = pos;
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
  end Sag_1;

  model Sag_ArbDeep
    inner Modelica.Mechanics.MultiBody.World world annotation(Placement(visible = true, transformation(origin = {-130, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    // Parameters
    parameter Integer n = 10 "Number of springs";
    parameter Integer m = n - 1 "Number of bodies";
    parameter Modelica.Units.SI.Position[3] f1_r = {15, 0, -7} "Left fixed point";
    parameter Modelica.Units.SI.Position[3] f2_r = {-15, 0, 7} "Right fixed point";
    parameter Modelica.Units.SI.TranslationalSpringConstant c = 100 "Spring constant";
    parameter Modelica.Units.SI.Mass m_t = 20;
    final parameter Modelica.Units.SI.Length totalDistance = Modelica.Math.Vectors.length(f1_r - f2_r);
    final parameter Modelica.Units.SI.Length s_unstretched = ((totalDistance / n) * 1.1);
    final parameter Modelica.Units.SI.Mass m_s = m_t / m;
    final parameter Modelica.Units.SI.Acceleration g = world.g "Gravity magnitude";
    final parameter Modelica.Units.SI.Position body_r_0[m, 3](each fixed = false);
    //
    //
    //
    final parameter Modelica.Units.SI.Position pos[m, 3] = Helpers.computeStaticPositions(m = m, L = totalDistance, c = c, m_s = m_s, g = g, s_unstretched = s_unstretched);
    //
    //
    final parameter Modelica.Units.SI.Position body_r_A[m, 3] = C_StaticPositions(n = n, f1_r = f1_r, f2_r = f2_r, totalMass = m_t, springConstant = c, s_unstretched = s_unstretched, gravity = g);
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

    function C_StaticPositions
      input Integer n;
      input Real f1[3], f2[3];
      input Real totalMass;
      input Real springConstant;
      input Real s_unstretched;
      input Real gravity;
      output Real pos[n - 1, 3];
    external "C"
      CStaticPositions_c(n, f1, f2, totalMass, springConstant, s_unstretched, gravity, pos)
        annotation(
          Library = "MySolverLib"
        );
    end C_StaticPositions;

    //
    //
  initial equation
    //body_r_0 = pos_global;
    body_r_0 = pos;
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
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {0.248, 0.044}, lineColor = {56, 56, 56}, fillColor = {128, 202, 255}, fillPattern = FillPattern.Solid, points = {{99.752, 100}, {99.752, 59.956}, {99.752, -50}, {100, -100}, {49.752, -100}, {-19.752, -100.044}, {-100.248, -100}, {-100.248, -50}, {-90.248, 29.956}, {-90.248, 79.956}, {-40.248, 79.956}, {-20.138, 79.813}, {-0.248, 79.956}, {19.752, 99.956}, {39.752, 99.956}, {59.752, 99.956}}, smooth = Smooth.Bezier), Polygon(visible = true, origin = {0, -13.079}, lineColor = {192, 192, 192}, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.HorizontalCylinder, points = {{100, -86.921}, {50, -86.921}, {-50, -86.921}, {-100, -86.921}, {-100, -36.921}, {-100, 53.079}, {-100, 103.079}, {-50, 103.079}, {0, 103.079}, {20, 83.079}, {50, 83.079}, {100, 83.079}, {100, 33.079}, {100, -36.921}}, smooth = Smooth.Bezier), Polygon(visible = true, origin = {0, -10.704}, lineColor = {113, 113, 113}, fillColor = {255, 255, 255}, points = {{100, -89.296}, {50, -89.296}, {-50, -89.296}, {-100, -89.296}, {-100, -39.296}, {-100, 50.704}, {-100, 100.704}, {-50, 100.704}, {0, 100.704}, {20, 80.704}, {50, 80.704}, {100, 80.704}, {100, 30.704}, {100, -39.296}}, smooth = Smooth.Bezier)}));
end Ropes;
