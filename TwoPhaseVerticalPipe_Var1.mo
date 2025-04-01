model TwoPhaseVerticalPipe_Var1 "Vertical two-phase pipe with replaceable phase transfer models"
  extends Modelica.Fluid.Interfaces.PartialTwoPort(redeclare package Medium = Modelica.Media.Water.StandardWater);
  import Modelica.Constants.g_n;
  import Modelica.Units.SI;
  replaceable model PhaseTransfer = PhaseTransferModels.Schrage constrainedby PhaseTransferModels.PartialPhaseTransfer "Phase change mechanism" annotation(choicesAllMatching = true);
  //
  //Geometry parameters
  parameter SI.Length L = 15 "Total pipe height" annotation(Dialog(tab = "General", group = "Geometry"));
  parameter SI.Diameter D = 0.1 "Pipe diameter" annotation(Dialog(tab = "General", group = "Geometry"));
  parameter SI.Area A_evap = 1 "Area of Phase Change" annotation(Dialog(tab = "General", group = "Geometry"));
  final parameter SI.Area A = Modelica.Constants.pi * (D / 2) ^ 2 "Cross Area";
  final SI.Height level (stateSelect= StateSelect.prefer)"Height of Water in Pipe";
  SI.Mass m_l "Liquid mass (state variable)";
  SI.Mass m_v "Vapor mass (state variable)";
  SI.Mass m_t "Total mass";
  final SI.Velocity v(start = 0, fixed = true) "Velocity";
  //SI.Acceleration a "Acceleration";
  // Heat Port
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort annotation(Placement(visible = true, transformation(origin = {0, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  final Medium.Temperature T_heatPort(start = 293.15) = heatPort.T "heatPort Temperature";
  SI.HeatFlowRate q_F = heatPort.Q_flow "heatPort QFlow";
  //
  SI.Pressure p_delta = port_a.p - port_b.p;
  // Energies
  final SI.Energy U( stateSelect= StateSelect.prefer);
  final SI.Temperature T_U;
  //
  // Instantiation of replacable model
  PhaseTransfer phaseTransfer(p_vapor = port_b.p, q_F = heatPort.Q_flow, m_l = m_l, T_l = T_U, A_evap = A_evap, A_rosh = A_rosh) "Instantiated phase transfer model";
  //
  //
  // Find the CHF Critical Heat Flow
  SI.Area A_rosh = Modelica.Constants.pi * D  * level;
  SI.HeatFlowRate CHF = A_rosh * 0.131 * (phaseTransfer.h_v - phaseTransfer.h_l) * phaseTransfer.rho_v * sqrt(g_n * (phaseTransfer.rho_l - phaseTransfer.rho_v) / phaseTransfer.sigma);
  //
initial equation
  level = 5;
  der(v) = 0;
  // Start from rest
equation
  // Mass definitions
  m_l = A * level * phaseTransfer.rho_l;
  m_v = A * (L - level) * phaseTransfer.rho_v;
  m_t = m_l + m_v;
  // Velocity, Acceleration
  v = der(level);
  //a = der(v);
  // Force balance
  //m_t * der(v) + v * der(m_t) = p_delta*A - m_t*g_n;
  //der(m_t * v) = (p_delta * A) - (m_t * g_n);
  der(v) = (p_delta*A - m_t*g_n - v*der(m_t)) / m_t;
  // connect heatPort.T to medium
  T_heatPort = T_U;
  //phaseTransfer.T_med;
  // Energy Balance
  der(U) = port_a.m_flow * actualStream(port_a.h_outflow) + port_b.m_flow * actualStream(port_b.h_outflow) + (q_F - phaseTransfer.q_phase);
  //U = m_v * phaseTransfer.h_v + m_l * phaseTransfer.h_l;
  T_U = U / ((m_l * phaseTransfer.cp_l) + (m_v * phaseTransfer.cp_v)) + 273.15;
  // Mass balance using phase transfer output
  der(m_l) = port_a.m_flow - phaseTransfer.m_dot_phase;
  der(m_v) = port_b.m_flow + phaseTransfer.m_dot_phase;

  //Enthalpy equations
  port_a.h_outflow = phaseTransfer.h_l;
  port_b.h_outflow = phaseTransfer.h_v;
  // Physical constraint
  assert(q_F < CHF, "Heat Flow above critical Heat Flow: q_F " +String(q_F) + "CHF: " + String(CHF));
  assert(level >= 0 and level <= L, "Water height out of bounds: " + String(level));
  annotation(Dialog(group = "Geometry"), uses(Modelica(version = "4.0.0")), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, origin = {-0.177, 0}, fillColor = {255, 255, 255}, extent = {{-80.177, -20}, {80.177, 20}}), Line(visible = true, origin = {-90, 0}, points = {{-10, 0}, {10, 0}}), Line(visible = true, origin = {90, 0}, points = {{-10, 0}, {10, 0}}), Rectangle(visible = true, origin = {-33.32, 0}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-46.68, -20}, {46.68, 20}})}), Documentation(info = "<html><div>
<h2>Corrected Vertical Two-Phase Pipe Model</h2>
</div>
<div>
<p>A vertical pipe that sucks a medium from the lower port up to a level determined by the gravity and the density of the medium.
 A heatport ist provided to introduce a heat-flow into the liquid phase.
 Phase change can be simulated by different phase change models)</p>
</div></html>"), Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {0, 114, 195}, fillColor = {255, 255, 255}, extent = {{-100, -100}, {100, 100}}, radius = 25), Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name")}));
end TwoPhaseVerticalPipe_Var1;