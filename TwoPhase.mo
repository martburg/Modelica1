package TwoPhase
  package Examples
    extends Modelica.Icons.ExamplesPackage;

    model Base_TwoPhaseExample
      extends Modelica.Icons.Example;
      import Modelica.Fluid.Sources.Boundary_pT;
      import Modelica.Media.Water.StandardWater;
      inner Modelica.Fluid.System system;
      Boundary_pT bottom(redeclare package Medium = StandardWater, p = 1.0e5, T = 358.15, nPorts = 1) annotation(Placement(visible = true, transformation(origin = {25, -70}, extent = {{-10, -10}, {10, 10}}, rotation = -270)));
      Boundary_pT top(redeclare package Medium = StandardWater, p = 0.05e5, T = 293.15, nPorts = 1, use_p_in = true, use_T_in = false) annotation(Placement(visible = true, transformation(origin = {25, 50}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Fluid.Valves.ValveLinear valveLinear1(redeclare package Medium = StandardWater, port_a_T.start = 293.15, port_b_T.start = 293.15, dp_nominal = 1000, m_flow_nominal = 100) annotation(Placement(visible = true, transformation(origin = {25, -30}, extent = {{-10, -10}, {10, 10}}, rotation = -270)));
      Modelica.Blocks.Sources.Constant const(k = 1) annotation(Placement(visible = true, transformation(origin = {-25, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp ramp1(duration = 200, offset = 85, height = 0, startTime = 100) annotation(Placement(visible = true, transformation(origin = {-65, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp ramp2(offset = 0, startTime = 400, height = 100000, duration = 500) annotation(Placement(visible = true, transformation(origin = {-65, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1 annotation(Placement(visible = true, transformation(origin = {-25, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Components.TwoPhaseVerticalPipe P_1(redeclare replaceable model PhaseTransfer = Components.PhaseTransferModels.Schrage) annotation(Placement(visible = true, transformation(origin = {25, 10}, extent = {{-10, -10}, {10, 10}}, rotation = -270)));
      Utilities.T_to_p_sat t_to_p_sat1 annotation(Placement(visible = true, transformation(origin = {-30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(valveLinear1.port_a, bottom.ports[1]) annotation(Line(visible = true, origin = {25, -50}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(const.y, valveLinear1.opening) annotation(Line(visible = true, origin = {1.5, -30}, points = {{-15.5, 0}, {15.5, 0}}, color = {0, 0, 127}));
      connect(ramp2.y, prescribedHeatFlow1.Q_flow) annotation(Line(visible = true, origin = {-44.5, 10}, points = {{-9.5, 0}, {9.5, 0}}, color = {0, 0, 127}));
      connect(top.ports[1], P_1.port_b) annotation(Line(visible = true, origin = {25, 30}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(P_1.port_a, valveLinear1.port_b) annotation(Line(visible = true, origin = {25, -10}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(P_1.heatPort, prescribedHeatFlow1.port) annotation(Line(visible = true, origin = {3.5, 10}, points = {{18.5, 0}, {-18.5, 0}}, color = {191, 0, 0}));
      connect(ramp1.y, t_to_p_sat1.T_target) annotation(Line(visible = true, origin = {-46.833, 70}, points = {{-7.167, 0}, {7.167, 0}}, color = {0, 0, 127}));
      connect(t_to_p_sat1.P_sat, top.p_in) annotation(Line(visible = true, origin = {14.889, 67.333}, points = {{-36.222, 2.667}, {18.111, 2.667}, {18.111, -5.333}}, color = {0, 0, 127}));
      annotation(uses(Modelica(version = "4.0.0")), Documentation, __Wolfram(ControlPanels(Panel(identifier = "charming-hardy", title = "TestPanel", elements = {InputField(variable = system.T_ambient)}))), experiment(StopTime = 3000, Interval = 0.1, Tolerance = 1e-6), Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {0, 114, 195}, fillColor = {255, 255, 255}, extent = {{-100, -100}, {100, 100}}, radius = 25), Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name")}));
    end Base_TwoPhaseExample;

    model Evap_Cond
      extends Modelica.Icons.Example;
      import Modelica.Fluid.Sources.Boundary_pT;
      import Modelica.Media.Water.StandardWater;
      inner Modelica.Fluid.System system;
      Boundary_pT bottom(redeclare package Medium = StandardWater, p = 1.0e5, T = 293.15, nPorts = 1, use_p_in = true) annotation(Placement(visible = true, transformation(origin = {-30, -20}, extent = {{10, -10}, {-10, 10}}, rotation = 270)));
      Modelica.Fluid.Valves.ValveLinear valveLinear1(redeclare package Medium = StandardWater, port_a_T.start = 293.15, port_b_T.start = 293.15, dp_nominal = 1000, m_flow_nominal = 5) annotation(Placement(visible = true, transformation(origin = {-30, 20}, extent = {{-10, -10}, {10, 10}}, rotation = -270)));
      Modelica.Blocks.Sources.Constant const(k = 1) annotation(Placement(visible = true, transformation(origin = {-70, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1(T_ref = 353.15) annotation(Placement(visible = true, transformation(origin = {-70, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Components.TwoPhaseVerticalPipe P_1(redeclare replaceable model PhaseTransfer = Components.PhaseTransferModels.Schrage, A_evap = 10, L = 15, D = 1) annotation(Placement(visible = true, transformation(origin = {-30, 60}, extent = {{-10, -10}, {10, 10}}, rotation = -270)));
      Boundary_pT bottom1(redeclare package Medium = StandardWater, p = 1.0e5, T = 353.15, nPorts = 1, use_p_in = true) annotation(Placement(visible = true, transformation(origin = {30, -20}, extent = {{-10, -10}, {10, 10}}, rotation = -270)));
      Modelica.Fluid.Valves.ValveLinear valveLinear2(redeclare package Medium = StandardWater, port_a_T.start = 293.15, port_b_T.start = 293.15, dp_nominal = 1000, m_flow_nominal = 5) annotation(Placement(visible = true, transformation(origin = {30, 20}, extent = {{10, -10}, {-10, 10}}, rotation = 270)));
      Modelica.Blocks.Sources.Constant const1(k = 1) annotation(Placement(visible = true, transformation(origin = {70, 20}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow2 annotation(Placement(visible = true, transformation(origin = {70, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Components.TwoPhaseVerticalPipe P_2(redeclare replaceable model PhaseTransfer = Components.PhaseTransferModels.Schrage, A_evap = 10, L = 15, D = 1) annotation(Placement(visible = true, transformation(origin = {30, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 270)));
      Modelica.Blocks.Math.UnitConversions.From_bar from_bar1 annotation(Placement(visible = true, transformation(origin = {-70, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.UnitConversions.From_bar from_bar2 annotation(Placement(visible = true, transformation(origin = {70, -60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Utilities.Ramp_L R_Lower_p(startTime = 10, offset = 0.1, duration = 100, height = 0.9) annotation(Placement(visible = true, transformation(origin = {110, -60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Utilities.Ramp_L L_Lower_p(startTime = 10, offset = 0.1, duration = 100, height = 0.9) annotation(Placement(visible = true, transformation(origin = {-110, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Utilities.Ramp_L L_Heat(height = 0, duration = 500, startTime = 500) annotation(Placement(visible = true, transformation(origin = {-110, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Utilities.Ramp_L R_Heat(height = 0, duration = 500, startTime = 1000) annotation(Placement(visible = true, transformation(origin = {110, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    equation
      connect(valveLinear1.port_a, bottom.ports[1]) annotation(Line(visible = true, origin = {-30, 0}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(const.y, valveLinear1.opening) annotation(Line(visible = true, origin = {-48.5, 20}, points = {{-10.5, 0}, {10.5, 0}}, color = {0, 0, 127}));
      connect(P_1.port_a, valveLinear1.port_b) annotation(Line(visible = true, origin = {-30, 40}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(P_1.heatPort, prescribedHeatFlow1.port) annotation(Line(visible = true, origin = {-46.5, 60}, points = {{13.5, 0}, {-13.5, 0}}, color = {191, 0, 0}));
      connect(valveLinear2.port_a, bottom1.ports[1]) annotation(Line(visible = true, origin = {30, 0}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(const1.y, valveLinear2.opening) annotation(Line(visible = true, origin = {48.5, 20}, points = {{10.5, 0}, {-10.5, 0}}, color = {0, 0, 127}));
      connect(P_2.port_a, valveLinear2.port_b) annotation(Line(visible = true, origin = {30, 40}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(P_2.heatPort, prescribedHeatFlow2.port) annotation(Line(visible = true, origin = {46.5, 60}, points = {{-13.5, 0}, {13.5, 0}}, color = {191, 0, 0}));
      connect(from_bar2.y, bottom1.p_in) annotation(Line(visible = true, origin = {34.333, -50.667}, points = {{24.667, -9.333}, {-12.333, -9.333}, {-12.333, 18.667}}, color = {0, 0, 127}));
      connect(from_bar1.y, bottom.p_in) annotation(Line(visible = true, origin = {-34.333, -50.667}, points = {{-24.667, -9.333}, {12.333, -9.333}, {12.333, 18.667}}, color = {0, 0, 127}));
      connect(R_Lower_p.y, from_bar2.u) annotation(Line(visible = true, origin = {90.5, -60}, points = {{8.5, 0}, {-8.5, 0}}, color = {0, 0, 127}));
      connect(L_Lower_p.y, from_bar1.u) annotation(Line(visible = true, origin = {-90.5, -60}, points = {{-8.5, 0}, {8.5, 0}}, color = {0, 0, 127}));
      connect(L_Heat.y, prescribedHeatFlow1.Q_flow) annotation(Line(visible = true, origin = {-89.5, 60}, points = {{-9.5, 0}, {9.5, 0}}, color = {0, 0, 127}));
      connect(prescribedHeatFlow2.Q_flow, R_Heat.y) annotation(Line(visible = true, origin = {89.5, 60}, points = {{-9.5, 0}, {9.5, 0}}, color = {0, 0, 127}));
      connect(P_1.port_b, P_2.port_b) annotation(Line(visible = true, origin = {0, 75}, points = {{-30, -5}, {-30, 5}, {30, 5}, {30, -5}}, color = {0, 127, 255}));
      annotation(uses(Modelica(version = "4.0.0")), Documentation(figures = {Figure(title = "Level_p_T", identifier = "affectionate-huygens", plots = {Plot(curves = {Curve(y = P_1.level), Curve(y = P_2.level)}), Plot(curves = {Curve(y = P_2.port_b.p), Curve(y = P_1.port_b.p)}), Plot(curves = {Curve(y = P_2.T_U), Curve(y = P_1.T_U)})})}), __Wolfram(ControlPanels(Panel(identifier = "charming-hardy", title = "TestPanel", elements = {InputField(variable = system.T_ambient)}), Panel(identifier = "loving-galois", title = "Untitled"))), experiment(StopTime = 3000, Interval = 0.1, Tolerance = 1e-6), Diagram(coordinateSystem(extent = {{-210, -148.5}, {210, 148.5}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {0, 114, 195}, fillColor = {255, 255, 255}, extent = {{-100, -100}, {100, 100}}, radius = 25), Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name")}));
    end Evap_Cond;

    model Peters_Desalination
      extends Modelica.Icons.Example;
      import Modelica.Fluid.Sources.Boundary_pT;
      import Modelica.Media.Water.StandardWater;
      inner Modelica.Fluid.System system;
      Boundary_pT bottom(redeclare package Medium = StandardWater, p = 1.0e5, T = 353.15, nPorts = 1, use_p_in = true) annotation(Placement(visible = true, transformation(origin = {-30, -20}, extent = {{10, -10}, {-10, 10}}, rotation = 270)));
      Modelica.Fluid.Valves.ValveLinear valveLinear1(redeclare package Medium = StandardWater, port_a_T.start = 293.15, port_b_T.start = 293.15, dp_nominal = 1000, m_flow_nominal = 5) annotation(Placement(visible = true, transformation(origin = {-30, 20}, extent = {{-10, -10}, {10, 10}}, rotation = -270)));
      Modelica.Blocks.Sources.Constant const(k = 1) annotation(Placement(visible = true, transformation(origin = {-70, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1(T_ref = 353.15) annotation(Placement(visible = true, transformation(origin = {-70, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Components.TwoPhaseVerticalPipe P_1(redeclare replaceable model PhaseTransfer = Components.PhaseTransferModels.Schrage, A_evap = 10, L = 15, D = 1) annotation(Placement(visible = true, transformation(origin = {-30, 60}, extent = {{-10, -10}, {10, 10}}, rotation = -270)));
      Boundary_pT bottom1(redeclare package Medium = StandardWater, p = 1.0e5, T = 353.15, nPorts = 1, use_p_in = true) annotation(Placement(visible = true, transformation(origin = {30, -20}, extent = {{-10, -10}, {10, 10}}, rotation = -270)));
      Modelica.Fluid.Valves.ValveLinear valveLinear2(redeclare package Medium = StandardWater, port_a_T.start = 293.15, port_b_T.start = 293.15, dp_nominal = 1000, m_flow_nominal = 5) annotation(Placement(visible = true, transformation(origin = {30, 20}, extent = {{10, -10}, {-10, 10}}, rotation = 270)));
      Modelica.Blocks.Sources.Constant const1(k = 1) annotation(Placement(visible = true, transformation(origin = {70, 20}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow2 annotation(Placement(visible = true, transformation(origin = {70, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Components.TwoPhaseVerticalPipe P_2(redeclare replaceable model PhaseTransfer = Components.PhaseTransferModels.Schrage, A_evap = 10, L = 15, D = 10) annotation(Placement(visible = true, transformation(origin = {30, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 270)));
      Modelica.Blocks.Math.UnitConversions.From_bar from_bar1 annotation(Placement(visible = true, transformation(origin = {-70, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.UnitConversions.From_bar from_bar2 annotation(Placement(visible = true, transformation(origin = {70, -60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Utilities.Ramp_L R_Lower_p(startTime = 10, offset = 0.1, duration = 100, height = 0.9) annotation(Placement(visible = true, transformation(origin = {110, -60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Utilities.Ramp_L L_Lower_p(startTime = 10, offset = 0.1, duration = 100, height = 0.9) annotation(Placement(visible = true, transformation(origin = {-110, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Utilities.Ramp_L L_Heat(height = 100e3, duration = 500, startTime = 500) annotation(Placement(visible = true, transformation(origin = {-110, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Utilities.Ramp_L R_Heat(height = -100e3, duration = 500, startTime = 1000) annotation(Placement(visible = true, transformation(origin = {110, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    equation
      connect(valveLinear1.port_a, bottom.ports[1]) annotation(Line(visible = true, origin = {-30, 0}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(const.y, valveLinear1.opening) annotation(Line(visible = true, origin = {-48.5, 20}, points = {{-10.5, 0}, {10.5, 0}}, color = {0, 0, 127}));
      connect(P_1.port_a, valveLinear1.port_b) annotation(Line(visible = true, origin = {-30, 40}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(P_1.heatPort, prescribedHeatFlow1.port) annotation(Line(visible = true, origin = {-46.5, 60}, points = {{13.5, 0}, {-13.5, 0}}, color = {191, 0, 0}));
      connect(valveLinear2.port_a, bottom1.ports[1]) annotation(Line(visible = true, origin = {30, 0}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(const1.y, valveLinear2.opening) annotation(Line(visible = true, origin = {48.5, 20}, points = {{10.5, 0}, {-10.5, 0}}, color = {0, 0, 127}));
      connect(P_2.port_a, valveLinear2.port_b) annotation(Line(visible = true, origin = {30, 40}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(P_2.heatPort, prescribedHeatFlow2.port) annotation(Line(visible = true, origin = {46.5, 60}, points = {{-13.5, 0}, {13.5, 0}}, color = {191, 0, 0}));
      connect(from_bar2.y, bottom1.p_in) annotation(Line(visible = true, origin = {34.333, -50.667}, points = {{24.667, -9.333}, {-12.333, -9.333}, {-12.333, 18.667}}, color = {0, 0, 127}));
      connect(from_bar1.y, bottom.p_in) annotation(Line(visible = true, origin = {-34.333, -50.667}, points = {{-24.667, -9.333}, {12.333, -9.333}, {12.333, 18.667}}, color = {0, 0, 127}));
      connect(R_Lower_p.y, from_bar2.u) annotation(Line(visible = true, origin = {90.5, -60}, points = {{8.5, 0}, {-8.5, 0}}, color = {0, 0, 127}));
      connect(L_Lower_p.y, from_bar1.u) annotation(Line(visible = true, origin = {-90.5, -60}, points = {{-8.5, 0}, {8.5, 0}}, color = {0, 0, 127}));
      connect(L_Heat.y, prescribedHeatFlow1.Q_flow) annotation(Line(visible = true, origin = {-89.5, 60}, points = {{-9.5, 0}, {9.5, 0}}, color = {0, 0, 127}));
      connect(prescribedHeatFlow2.Q_flow, R_Heat.y) annotation(Line(visible = true, origin = {89.5, 60}, points = {{-9.5, 0}, {9.5, 0}}, color = {0, 0, 127}));
      connect(P_1.port_b, P_2.port_b) annotation(Line(visible = true, origin = {0, 75}, points = {{-30, -5}, {-30, 5}, {30, 5}, {30, -5}}, color = {0, 127, 255}));
      annotation(uses(Modelica(version = "4.0.0")), Documentation(figures = {Figure(title = "Level_p_T", identifier = "affectionate-huygens", plots = {Plot(curves = {Curve(y = P_1.level), Curve(y = P_2.level)}), Plot(curves = {Curve(y = P_2.port_b.p), Curve(y = P_1.port_b.p)}), Plot(curves = {Curve(y = P_2.T_U), Curve(y = P_1.T_U)})}), Figure(title = "Level_p_T1_mflow", identifier = "adoring-fibonacci", plots = {Plot(curves = {Curve(y = P_1.level), Curve(y = P_2.level)}), Plot(curves = {Curve(y = P_2.port_b.p), Curve(y = P_1.port_b.p)}), Plot(curves = {Curve(y = P_2.T_U), Curve(y = P_1.T_U)}), Plot(curves = {Curve(y = P_1.port_b.m_flow)})})}), __Wolfram(ControlPanels(Panel(identifier = "charming-hardy", title = "TestPanel", elements = {InputField(variable = system.T_ambient)}), Panel(identifier = "loving-galois", title = "Untitled"))), experiment(StopTime = 3000, Interval = 0.1, Tolerance = 1e-6), Diagram(coordinateSystem(extent = {{-210, -148.5}, {210, 148.5}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {0, 114, 195}, fillColor = {255, 255, 255}, extent = {{-100, -100}, {100, 100}}, radius = 25), Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name")}));
    end Peters_Desalination;

    model Evap_Pipe_Cond
      extends Modelica.Icons.Example;
      import Modelica.Fluid.Sources.Boundary_pT;
      import Modelica.Media.Water.StandardWater;
      inner Modelica.Fluid.System system;
      Boundary_pT bottom(redeclare package Medium = StandardWater, p = 1.0e5, T = 353.15, nPorts = 1, use_p_in = true) annotation(Placement(visible = true, transformation(origin = {-30, -20}, extent = {{10, -10}, {-10, 10}}, rotation = 270)));
      Modelica.Fluid.Valves.ValveLinear valveLinear1(redeclare package Medium = StandardWater, port_a_T.start = 293.15, port_b_T.start = 293.15, dp_nominal = 1000, m_flow_nominal = 25) annotation(Placement(visible = true, transformation(origin = {-30, 20}, extent = {{-10, -10}, {10, 10}}, rotation = -270)));
      Modelica.Blocks.Sources.Constant const(k = 1) annotation(Placement(visible = true, transformation(origin = {-70, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1(T_ref = 353.15) annotation(Placement(visible = true, transformation(origin = {-70, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Components.TwoPhaseVerticalPipe P_1(redeclare replaceable model PhaseTransfer = Components.PhaseTransferModels.Schrage, A_evap = 100, L = 15, D = 1) annotation(Placement(visible = true, transformation(origin = {-30, 60}, extent = {{-10, -10}, {10, 10}}, rotation = -270)));
      Boundary_pT bottom1(redeclare package Medium = StandardWater, p = 1.0e5, T = 353.15, nPorts = 1, use_p_in = true) annotation(Placement(visible = true, transformation(origin = {30, -20}, extent = {{-10, -10}, {10, 10}}, rotation = -270)));
      Modelica.Fluid.Valves.ValveLinear valveLinear2(redeclare package Medium = StandardWater, port_a_T.start = 293.15, port_b_T.start = 293.15, dp_nominal = 1000, m_flow_nominal = 25) annotation(Placement(visible = true, transformation(origin = {30, 20}, extent = {{10, -10}, {-10, 10}}, rotation = 270)));
      Modelica.Blocks.Sources.Constant const1(k = 1) annotation(Placement(visible = true, transformation(origin = {70, 20}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow2 annotation(Placement(visible = true, transformation(origin = {70, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Components.TwoPhaseVerticalPipe P_2(redeclare replaceable model PhaseTransfer = Components.PhaseTransferModels.Schrage, A_evap = 100, L = 15, D = 1) annotation(Placement(visible = true, transformation(origin = {30, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 270)));
      Modelica.Blocks.Math.UnitConversions.From_bar from_bar1 annotation(Placement(visible = true, transformation(origin = {-70, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.UnitConversions.From_bar from_bar2 annotation(Placement(visible = true, transformation(origin = {70, -60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Utilities.Ramp_L R_Lower_p(startTime = 10, offset = 0.1, duration = 100, height = 0.9) annotation(Placement(visible = true, transformation(origin = {110, -60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Utilities.Ramp_L L_Lower_p(startTime = 10, offset = 0.1, duration = 100, height = 0.9) annotation(Placement(visible = true, transformation(origin = {-110, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Utilities.Ramp_L L_Heat(height = 1e7, duration = 1000, startTime = 500) annotation(Placement(visible = true, transformation(origin = {-110, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Utilities.Ramp_L R_Heat(height = -1.1e7, duration = 1000, startTime = 500) annotation(Placement(visible = true, transformation(origin = {110, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Fluid.Pipes.StaticPipe pipe(redeclare replaceable package Medium = StandardWater, length = 1, diameter = 0.1) annotation(Placement(visible = true, transformation(origin = {0, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(valveLinear1.port_a, bottom.ports[1]) annotation(Line(visible = true, origin = {-30, 0}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(const.y, valveLinear1.opening) annotation(Line(visible = true, origin = {-48.5, 20}, points = {{-10.5, 0}, {10.5, 0}}, color = {0, 0, 127}));
      connect(P_1.port_a, valveLinear1.port_b) annotation(Line(visible = true, origin = {-30, 40}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(P_1.heatPort, prescribedHeatFlow1.port) annotation(Line(visible = true, origin = {-46.5, 60}, points = {{13.5, 0}, {-13.5, 0}}, color = {191, 0, 0}));
      connect(valveLinear2.port_a, bottom1.ports[1]) annotation(Line(visible = true, origin = {30, 0}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(const1.y, valveLinear2.opening) annotation(Line(visible = true, origin = {48.5, 20}, points = {{10.5, 0}, {-10.5, 0}}, color = {0, 0, 127}));
      connect(P_2.port_a, valveLinear2.port_b) annotation(Line(visible = true, origin = {30, 40}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(P_2.heatPort, prescribedHeatFlow2.port) annotation(Line(visible = true, origin = {46.5, 60}, points = {{-13.5, 0}, {13.5, 0}}, color = {191, 0, 0}));
      connect(from_bar2.y, bottom1.p_in) annotation(Line(visible = true, origin = {34.333, -50.667}, points = {{24.667, -9.333}, {-12.333, -9.333}, {-12.333, 18.667}}, color = {0, 0, 127}));
      connect(from_bar1.y, bottom.p_in) annotation(Line(visible = true, origin = {-34.333, -50.667}, points = {{-24.667, -9.333}, {12.333, -9.333}, {12.333, 18.667}}, color = {0, 0, 127}));
      connect(R_Lower_p.y, from_bar2.u) annotation(Line(visible = true, origin = {90.5, -60}, points = {{8.5, 0}, {-8.5, 0}}, color = {0, 0, 127}));
      connect(L_Lower_p.y, from_bar1.u) annotation(Line(visible = true, origin = {-90.5, -60}, points = {{-8.5, 0}, {8.5, 0}}, color = {0, 0, 127}));
      connect(L_Heat.y, prescribedHeatFlow1.Q_flow) annotation(Line(visible = true, origin = {-89.5, 60}, points = {{-9.5, 0}, {9.5, 0}}, color = {0, 0, 127}));
      connect(prescribedHeatFlow2.Q_flow, R_Heat.y) annotation(Line(visible = true, origin = {89.5, 60}, points = {{-9.5, 0}, {9.5, 0}}, color = {0, 0, 127}));
      connect(pipe.port_a, P_1.port_b) annotation(Line(visible = true, origin = {-23.333, 76.667}, points = {{13.333, 3.333}, {-6.667, 3.333}, {-6.667, -6.667}}, color = {0, 127, 255}));
      connect(pipe.port_b, P_2.port_b) annotation(Line(visible = true, origin = {23.333, 76.667}, points = {{-13.333, 3.333}, {6.667, 3.333}, {6.667, -6.667}}, color = {0, 127, 255}));
      annotation(uses(Modelica(version = "4.0.0")), Documentation(figures = {Figure(title = "Level_p_T", identifier = "affectionate-huygens", plots = {Plot(curves = {Curve(y = P_1.level), Curve(y = P_2.level)}), Plot(curves = {Curve(y = P_2.port_b.p), Curve(y = P_1.port_b.p)}), Plot(curves = {Curve(y = P_2.T_U), Curve(y = P_1.T_U)})})}, info = "<html><div>
<div>
<p>parameter P_1.L is set to &lt;b&gt;String(P1_L)&lt;/b&gt;</p>
</div>
</div></html>"), __Wolfram(ControlPanels(Panel(identifier = "charming-hardy", title = "TestPanel", elements = {InputField(variable = system.T_ambient)}), Panel(identifier = "loving-galois", title = "Untitled"))), experiment(StopTime = 3000, Interval = 0.1, Tolerance = 1e-6), Diagram(coordinateSystem(extent = {{-210, -148.5}, {210, 148.5}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {0, 114, 195}, fillColor = {255, 255, 255}, extent = {{-100, -100}, {100, 100}}, radius = 25), Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name")}));
    end Evap_Pipe_Cond;

    model Evap_Pipe_Cond_E7
      extends Modelica.Icons.Example;
      import Modelica.Fluid.Sources.Boundary_pT;
      import Modelica.Media.Water.StandardWater;
      inner Modelica.Fluid.System system;
      Boundary_pT bottom(redeclare package Medium = StandardWater, p = 1.0e5, T = 353.15, nPorts = 1, use_p_in = true) annotation(Placement(visible = true, transformation(origin = {-30, -20}, extent = {{10, -10}, {-10, 10}}, rotation = 270)));
      Modelica.Fluid.Valves.ValveLinear valveLinear1(redeclare package Medium = StandardWater, port_a_T.start = 293.15, port_b_T.start = 293.15, dp_nominal = 1000, m_flow_nominal = 25) annotation(Placement(visible = true, transformation(origin = {-30, 20}, extent = {{-10, -10}, {10, 10}}, rotation = -270)));
      Modelica.Blocks.Sources.Constant const(k = 1) annotation(Placement(visible = true, transformation(origin = {-70, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1(T_ref = 353.15) annotation(Placement(visible = true, transformation(origin = {-70, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Components.TwoPhaseVerticalPipe P_1(redeclare replaceable model PhaseTransfer = Components.PhaseTransferModels.Schrage, A_evap = 100, L = 15, D = 1) annotation(Placement(visible = true, transformation(origin = {-30, 60}, extent = {{-10, -10}, {10, 10}}, rotation = -270)));
      Boundary_pT bottom1(redeclare package Medium = StandardWater, p = 1.0e5, T = 353.15, nPorts = 1, use_p_in = true) annotation(Placement(visible = true, transformation(origin = {30, -20}, extent = {{-10, -10}, {10, 10}}, rotation = -270)));
      Modelica.Fluid.Valves.ValveLinear valveLinear2(redeclare package Medium = StandardWater, port_a_T.start = 293.15, port_b_T.start = 293.15, dp_nominal = 1000, m_flow_nominal = 25) annotation(Placement(visible = true, transformation(origin = {30, 20}, extent = {{10, -10}, {-10, 10}}, rotation = 270)));
      Modelica.Blocks.Sources.Constant const1(k = 1) annotation(Placement(visible = true, transformation(origin = {70, 20}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow2 annotation(Placement(visible = true, transformation(origin = {70, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Components.TwoPhaseVerticalPipe P_2(redeclare replaceable model PhaseTransfer = Components.PhaseTransferModels.Schrage, A_evap = 100, L = 15, D = 1) annotation(Placement(visible = true, transformation(origin = {30, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 270)));
      Modelica.Blocks.Math.UnitConversions.From_bar from_bar1 annotation(Placement(visible = true, transformation(origin = {-70, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.UnitConversions.From_bar from_bar2 annotation(Placement(visible = true, transformation(origin = {70, -60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Utilities.Ramp_L R_Lower_p(startTime = 10, offset = 0.1, duration = 100, height = 0.9) annotation(Placement(visible = true, transformation(origin = {110, -60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Utilities.Ramp_L L_Lower_p(startTime = 10, offset = 0.1, duration = 100, height = 0.9) annotation(Placement(visible = true, transformation(origin = {-110, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Utilities.Ramp_L L_Heat(height = 1e7, duration = 1000, startTime = 500) annotation(Placement(visible = true, transformation(origin = {-110, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Utilities.Ramp_L R_Heat(height = -1.1e7, duration = 1000, startTime = 500) annotation(Placement(visible = true, transformation(origin = {110, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Fluid.Pipes.StaticPipe pipe(redeclare replaceable package Medium = StandardWater, length = 1, diameter = 0.1) annotation(Placement(visible = true, transformation(origin = {0, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(valveLinear1.port_a, bottom.ports[1]) annotation(Line(visible = true, origin = {-30, 0}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(const.y, valveLinear1.opening) annotation(Line(visible = true, origin = {-48.5, 20}, points = {{-10.5, 0}, {10.5, 0}}, color = {0, 0, 127}));
      connect(P_1.port_a, valveLinear1.port_b) annotation(Line(visible = true, origin = {-30, 40}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(P_1.heatPort, prescribedHeatFlow1.port) annotation(Line(visible = true, origin = {-46.5, 60}, points = {{13.5, 0}, {-13.5, 0}}, color = {191, 0, 0}));
      connect(valveLinear2.port_a, bottom1.ports[1]) annotation(Line(visible = true, origin = {30, 0}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(const1.y, valveLinear2.opening) annotation(Line(visible = true, origin = {48.5, 20}, points = {{10.5, 0}, {-10.5, 0}}, color = {0, 0, 127}));
      connect(P_2.port_a, valveLinear2.port_b) annotation(Line(visible = true, origin = {30, 40}, points = {{0, 10}, {0, -10}}, color = {0, 127, 255}));
      connect(P_2.heatPort, prescribedHeatFlow2.port) annotation(Line(visible = true, origin = {46.5, 60}, points = {{-13.5, 0}, {13.5, 0}}, color = {191, 0, 0}));
      connect(from_bar2.y, bottom1.p_in) annotation(Line(visible = true, origin = {34.333, -50.667}, points = {{24.667, -9.333}, {-12.333, -9.333}, {-12.333, 18.667}}, color = {0, 0, 127}));
      connect(from_bar1.y, bottom.p_in) annotation(Line(visible = true, origin = {-34.333, -50.667}, points = {{-24.667, -9.333}, {12.333, -9.333}, {12.333, 18.667}}, color = {0, 0, 127}));
      connect(R_Lower_p.y, from_bar2.u) annotation(Line(visible = true, origin = {90.5, -60}, points = {{8.5, 0}, {-8.5, 0}}, color = {0, 0, 127}));
      connect(L_Lower_p.y, from_bar1.u) annotation(Line(visible = true, origin = {-90.5, -60}, points = {{-8.5, 0}, {8.5, 0}}, color = {0, 0, 127}));
      connect(L_Heat.y, prescribedHeatFlow1.Q_flow) annotation(Line(visible = true, origin = {-89.5, 60}, points = {{-9.5, 0}, {9.5, 0}}, color = {0, 0, 127}));
      connect(prescribedHeatFlow2.Q_flow, R_Heat.y) annotation(Line(visible = true, origin = {89.5, 60}, points = {{-9.5, 0}, {9.5, 0}}, color = {0, 0, 127}));
      connect(pipe.port_a, P_1.port_b) annotation(Line(visible = true, origin = {-23.333, 76.667}, points = {{13.333, 3.333}, {-6.667, 3.333}, {-6.667, -6.667}}, color = {0, 127, 255}));
      connect(pipe.port_b, P_2.port_b) annotation(Line(visible = true, origin = {23.333, 76.667}, points = {{-13.333, 3.333}, {6.667, 3.333}, {6.667, -6.667}}, color = {0, 127, 255}));
      annotation(uses(Modelica(version = "4.0.0")), Documentation(figures = {Figure(title = "Level_p_T", identifier = "affectionate-huygens", plots = {Plot(curves = {Curve(y = P_1.level), Curve(y = P_2.level)}), Plot(curves = {Curve(y = P_2.port_b.p), Curve(y = P_1.port_b.p)}), Plot(curves = {Curve(y = P_2.T_U), Curve(y = P_1.T_U)})})}), __Wolfram(ControlPanels(Panel(identifier = "charming-hardy", title = "TestPanel", elements = {InputField(variable = system.T_ambient)}), Panel(identifier = "loving-galois", title = "Untitled"))), experiment(StopTime = 3000, Interval = 0.1, Tolerance = 1e-6), Diagram(coordinateSystem(extent = {{-210, -148.5}, {210, 148.5}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {0, 114, 195}, fillColor = {255, 255, 255}, extent = {{-100, -100}, {100, 100}}, radius = 25), Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name")}));
    end Evap_Pipe_Cond_E7;
  end Examples;

  package Components
    extends Modelica.Icons.VariantsPackage;

    model TwoPhaseVerticalPipe "Vertical two-phase pipe with replaceable phase transfer models"
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
      final SI.Height level "Height of Water in Pipe";
      SI.Mass m_l "Liquid mass (state variable)";
      SI.Mass m_v "Vapor mass (state variable)";
      SI.Mass m_t "Total mass";
      final SI.Velocity v(start = 0, fixed = true) "Velocity";
      SI.Acceleration a "Acceleration";
      // Heat Port
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort annotation(Placement(visible = true, transformation(origin = {0, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      final Medium.Temperature T_heatPort(start = 293.15) = heatPort.T "heatPort Temperature";
      SI.HeatFlowRate q_F = heatPort.Q_flow "heatPort QFlow";
      //
      SI.Pressure p_delta = port_a.p - port_b.p;
      // Energies
      final SI.Energy U;
      SI.Energy U_in;
      SI.Energy U_out;
      SI.EnergyFlowRate U_in_flow;
      SI.EnergyFlowRate U_out_flow;
      final SI.Temperature T_U;
      //
      // Instantiation of replacable model
      PhaseTransfer phaseTransfer(p_vapor = port_b.p, q_F = heatPort.Q_flow, m_l = m_l, T_l = T_U, A_evap = A_evap) "Instantiated phase transfer model";
      //
      //
    initial equation
      port_a.p = port_b.p + phaseTransfer.rho_l * g_n * level + phaseTransfer.rho_v * g_n * (L - level);
      U_out = 21.02 "Inner Energy at 20° and 1 bar";
      U_in = 21.02 "Inner Energy at 20° and 1 bar";
      T_U = 5 + 273.5;
    equation
      // Mass definitions
      m_l = A * level * phaseTransfer.rho_l;
      m_v = A * (L - level) * phaseTransfer.rho_v;
      m_t = m_l + m_v;
      // Velocity, Acceleration
      v = der(level);
      a = der(v);
      // Force balance
      m_t * a = (p_delta * A) - (m_t * g_n);
      // connect heatPort.T to medium
      T_heatPort = phaseTransfer.T_med;
      // Energy Balance
      der(U_in) = U_in_flow;
      der(U_out) = U_out_flow;
      U_in_flow = port_a.m_flow * actualStream(port_a.h_outflow);
      U_out_flow = port_b.m_flow * actualStream(port_b.h_outflow);
      der(U) = U_in_flow + U_out_flow + (q_F - phaseTransfer.q_phase);
      T_U = U / ((m_l * phaseTransfer.cp_l) + (m_v * phaseTransfer.cp_v)) + 273.15;
      // Mass balance using phase transfer output
      der(m_l) = port_a.m_flow - phaseTransfer.m_dot_phase;
      der(m_v) = port_b.m_flow + phaseTransfer.m_dot_phase;
      //Enthalpy equations
      port_a.h_outflow = phaseTransfer.h_l;
      port_b.h_outflow = phaseTransfer.h_v;
      // Physical constraint
      //assert(level >= 0 and level <= L, "Water height out of bounds: " + String(level));
      annotation(Dialog(group = "Geometry"), uses(Modelica(version = "4.0.0")), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, origin = {-0.177, 0}, fillColor = {255, 255, 255}, extent = {{-80.177, -20}, {80.177, 20}}), Line(visible = true, origin = {-90, 0}, points = {{-10, 0}, {10, 0}}), Line(visible = true, origin = {90, 0}, points = {{-10, 0}, {10, 0}}), Rectangle(visible = true, origin = {-33.32, 0}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-46.68, -20}, {46.68, 20}})}), Documentation(info = "<html><div>
<h2>Corrected Vertical Two-Phase Pipe Model</h2>
</div>
<div>
<p>A vertical pipe that sucks a medium from the lower port up to a level determined by the gravity and the density of the medium.
 A heatport ist provided to introduce a heat-flow into the liquid phase.
 Phase change can be simulated by different phase change models)</p>
</div></html>"), Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {0, 114, 195}, fillColor = {255, 255, 255}, extent = {{-100, -100}, {100, 100}}, radius = 25), Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name")}));
    end TwoPhaseVerticalPipe;

    package PhaseTransferModels
      extends Modelica.Icons.Package;

      partial model PartialPhaseTransfer
        replaceable package Medium = Modelica.Media.Water.StandardWater constrainedby Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching = true);
        import Modelica.Constants.g_n;
        import Modelica.Units.SI;
        // Inputs from main model
        input SI.Pressure p_vapor "Vapor pressure (port_b.p)";
        input SI.HeatFlowRate q_F "Total heat input";
        input SI.Mass m_l "Liquid Mass";
        input SI.Temperature T_l "Temperature";
        input SI.Area A_evap = 1 "Evaporation area";
        // Calulating from state leads to the system not beeing index reducible
        parameter SI.SpecificHeatCapacity cp_l = 4186 "Liquid specific heat";
        parameter SI.SpecificHeatCapacity cp_v = 1900 "Vapor specific heat";
        parameter SI.MolarMass MM = 0.01801528 "Molar Mass of Water";
        // Calculating from state leads to shattering so as an approx value in the
        parameter SI.DynamicViscosity mu_l = 1.0e-3 "Dynamic Viscosity";
        // Saturation Properties of Vapor at p_vapor
        Medium.SaturationProperties sat_v = Medium.setSat_p(p_vapor) "Saturation Properties at Saturation pressure ";
        SI.SpecificEnthalpy h_v = Medium.dewEnthalpy(sat_v) "Vapor  Specific Enthalpy";
        SI.Density rho_v = Medium.dewDensity(sat_v) "Vapor Density";
        SI.Temperature Tsat_v = sat_v.Tsat "Vapor Temperature";
        SI.Pressure p_sat_v = sat_v.psat "Vapor Pressure";
        // Saturation Properties of Liquid at T_U
        Medium.SaturationProperties sat_l = Medium.setSat_T(T_l) "Saturation Properties at Saturation Temperature";
        SI.SpecificEnthalpy h_l = Medium.bubbleEnthalpy(sat_l) "Liquid Specific Enthalpy";
        SI.Density rho_l = Medium.bubbleDensity(sat_l) "Liquid Density ";
        SI.Pressure p_sat_l = sat_l.psat "Liquid Pressure";
        SI.SurfaceTension sigma = Medium.surfaceTension(sat_l) "Surface Tension";
        Medium.ThermodynamicState state_liquid = Medium.setState_pT(p_sat_l, T_l);
        //SI.DynamicViscosity mu_l = Medium.dynamicViscosity(state_liquid) " Liquid Viscosity";
        SI.ThermalConductivity k_l = Medium.thermalConductivity(state_liquid) "Liquid Thermal conductivity";
        //SI.SpecificHeatCapacity cp_l = delay(Medium.specificHeatCapacityCp(state_liquid),0.5);
        // Outputs to main model
        output SI.MassFlowRate m_dot_phase "Phase change mass flow";
        output SI.HeatFlowRate q_phase "Heat used for phase change";
        output SI.Temperature T_med = Medium.saturationTemperature(p_sat_v);
        //equation
        //m_dot_phase =0;
        //q_phase = 0;
      end PartialPhaseTransfer;

      model Evaporation "Evaporation"
        extends PartialPhaseTransfer;
        parameter SI.LinearTemperatureCoefficient k = 1000 "Heat transfer coeff";
      equation
        q_phase = k * A_evap * (T_l - Tsat_v);
        m_dot_phase = q_phase / (h_v - h_l);
      end Evaporation;

      model Energy "Energy-driven phase change"
        extends PartialPhaseTransfer;
      equation
        if T_l < Tsat_v then
          // No phase change below saturation
          q_phase = 0;
          m_dot_phase = 0;
        else
          // All heat contributes to phase change
          q_phase = q_F;
          m_dot_phase = q_phase / (h_v - h_l);
        end if;
      end Energy;

      model NoChange "No Phase Change"
        extends PartialPhaseTransfer;
      equation
        q_phase = 0;
        m_dot_phase = 0;
      end NoChange;

      model Roshenov "Roshenov's boiling model"
        extends PartialPhaseTransfer;
        import Modelica.Constants.g_n;
        parameter Real C_boiling = 0.013 "Rohsenow constant";
        parameter Real n_boiling = 1 "Exponent";
        parameter Real C_sf = 0.006 "Surface-fluid Factor";
        parameter Real g0 = 1 "force conversation factor kgm/Ns²";
        // q_f_roshenov shatters if we calulate Pr so as an approx we set it to fixed value
        //Real                                                                                Pr = Medium.prandtlNumber(state_liquid) "Prandtl Number";
        Real Pr = 1.76;
      equation
        q_phase = A_evap * ((cp_l * (T_l - Tsat_v)) / ((h_v - h_l) * Pr ^ n_boiling * C_sf)) ^ 3 * mu_l * (h_v - h_l) * sqrt((g_n * (rho_l - rho_v)) / (g0 * sigma));
        //q_phase = C_boiling * sqrt(mu_l * g_n / sigma) * (h_v - h_l) * Pr ^ n_boiling * (T_l - Tsat_v) ^ 3;
        m_dot_phase = q_phase / (h_v - h_l);
        annotation(Documentation(info = "<html><div>
<div>
<div>
<p>information</p>
</div>
</div>
</div>
<div>
<p><a href=\"https://www.nuclear-power.com/nuclear-engineering/heat-transfer/boiling-and-condensation/nucleate-boiling-correlations-rohsenow-correlation/\">Nucleat Boiling</a></p>
<p>&nbsp;</p>
<p><img src=\"https://nuclear-power.com/wp-content/uploads/2018/01/Rohsenow-correlation-nucleate-boiling.png\" width=\"403\" height=\"87\" /></p>
<p>&nbsp;</p>
<p>where</p>
<ul>
<li>q &ndash; nucleate pool boiling heat flux [W/m<sup>2</sup>]</li>
<li>c<sub>1</sub>&nbsp;&mdash; specific heat of liquid J/kg K</li>
<li>&Delta;T &mdash; excess temperature &deg;C or K</li>
<li>h<sub>fg</sub>&nbsp;&nbsp;&ndash; enthalpy of vaporization, J/kg</li>
<li>Pr &mdash; Prandtl number of liquid</li>
<li>n &mdash; experimental constant equal to 1 for water and 1.7 for other fluids</li>
<li>C<sub>sf</sub>&nbsp;&mdash; surface fluid factor, for example, water and nickel have a C<sub>sf</sub>&nbsp;of 0.006</li>
<li>&mu;<sub>1</sub>&nbsp;&mdash; dynamic viscosity of the liquid kg/m.s</li>
<li>g &ndash; gravitational acceleration m/s<sup>2</sup></li>
<li>g<sub>0</sub>&nbsp;&mdash; force conversion factor kgm/Ns<sup>2</sup></li>
<li>&rho;<sub>1</sub>&nbsp;&mdash; density of the liquid kg/m<sup>3</sup></li>
<li>&rho;<sub>v</sub>&nbsp;&mdash; density of vapor kg/m<sup>3</sup></li>
<li>&sigma; &mdash; surface tension-liquid-vapor interface N/m</li>
</ul>
<p>As can be seen, &Delta;T &prop; (q)<sup>⅓</sup>. This very important proportionality shows increasing ability of interface to transfer heat.</p>
<p>&nbsp;</p>
<p>model NucleateBoiling<br />&nbsp; parameter Real c1 \"Specific heat of liquid [J/kg.K]\";<br />&nbsp; parameter Real deltaT \"Excess temperature [K]\";<br />&nbsp; parameter Real h_fg \"Enthalpy of vaporization [J/kg]\";<br />&nbsp; parameter Real Pr \"Prandtl number of liquid\";<br />&nbsp; parameter Real n \"Experimental constant\";<br />&nbsp; parameter Real C_sf \"Surface-fluid factor\";<br />&nbsp; parameter Real mu1 \"Dynamic viscosity of liquid [kg/m.s]\";<br />&nbsp; parameter Real g \"Gravitational acceleration [m/s^2]\";<br />&nbsp; parameter Real g0 \"Force conversion factor [kgm/Ns^2]\";<br />&nbsp; parameter Real rho1 \"Density of liquid [kg/m^3]\";<br />&nbsp; parameter Real rho_v \"Density of vapor [kg/m^3]\";<br />&nbsp; parameter Real sigma \"Surface tension [N/m]\";<br />&nbsp; Real q \"Nucleate boiling heat flux [W/m^2]\";<br />equation<br />&nbsp; q = ((c1 * deltaT) / (h_fg * Pr^n * C_sf))^3 * mu1 * h_fg * sqrt((g * (rho1 - rho_v)) / (g0 * sigma));<br />end NucleateBoiling;</p>
</div></html>"));
      end Roshenov;

      model Schrage "Schrage's equation"
        extends PartialPhaseTransfer;
        parameter Real C_evap = 0.8 "Evaporation coefficient";
      equation
        m_dot_phase = (2 * C_evap / (2 - C_evap)) * ((p_sat_l - p_vapor) * A_evap) * MM / sqrt(2 * Modelica.Constants.pi * Modelica.Constants.R * T_l * MM);
        q_phase = m_dot_phase * (h_v - h_l);
        // Heat used for evaporation
        annotation(Documentation(info = "<html><h1 class=\"wi-article-title article-title-main\">Revisiting the Schrage Equation for Kinetically Limited Evaporation and Condensation&nbsp;</h1>
<div class=\"wi-authors\">
<div class=\"al-authors-list\">
<div class=\"al-author-name\"><a class=\"linked-name js-linked-name stats-author-info-trigger\">Geoffrey Vaartstra</a><span class=\"al-author-delim\">,</span></div>
&nbsp;
<div class=\"al-author-name\"><a class=\"linked-name js-linked-name stats-author-info-trigger\">Zhengmao Lu</a><span class=\"al-author-delim\">,</span></div>
&nbsp;
<div class=\"al-author-name\"><a class=\"linked-name js-linked-name stats-author-info-trigger\">John H. Lienhard</a><span class=\"al-author-delim\">,</span></div>
&nbsp;
<div class=\"al-author-name\"><a class=\"linked-name js-linked-name stats-author-info-trigger\">Evelyn N. Wang</a></div>
</div>
</div>
<p>&nbsp;</p>
<p>Information</p>
<div class=\"pub-history-row clearfix\">
<div class=\"ww-citation-primary\"><em>J. Heat Transfer</em>. Aug 2022, 144(8): 080802 (21 pages)</div>
</div>
<div class=\"ww-citation-legacy-id\"><span class=\"ww-citation-legacy-label\">Paper No:&nbsp;</span><span class=\"ww-citation-legacy-article-id\">HT-21-1860</span>&nbsp;<a class=\"ww-doi-link\" href=\"https://doi.org/10.1115/1.4054382\" target=\"_blank\" rel=\"noopener\">https://doi.org/10.1115/1.4054382</a></div>
<div class=\"pub-history-row clearfix\">
<div class=\"ww-publish-date-wrap\"><span class=\"publish-date-label\"><strong>Published Online:</strong>&nbsp;May 24, 2022</span></div>
<div class=\"ww-publish-date-wrap\">&nbsp;</div>
</div>
<div class=\"wi-authors\">
<div class=\"al-authors-list\">
<div class=\"al-author-name\">&nbsp;</div>
</div>
</div></html>"));
      end Schrage;
    end PhaseTransferModels;
  end Components;

  package Utilities
    extends Modelica.Icons.UtilitiesPackage;

    model T_to_p_sat
      //extends Modelica.Icons.Function;
      package Medium = Modelica.Media.Water.StandardWater;
      // Calculate Psat of Water with input Temperature
      Modelica.Blocks.Interfaces.RealInput T_target annotation(Placement(visible = true, transformation(origin = {-150, -0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-96.667, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Medium.SaturationProperties sat_target = Medium.setSat_T(T_target + 273.15) "Saturation Properties at T_target";
      Modelica.Blocks.Interfaces.RealOutput P_sat = sat_target.psat annotation(Placement(visible = true, transformation(origin = {150, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {86.667, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      annotation(uses(Modelica(version = "4.0.0")), Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, textColor = {0, 0, 255}, extent = {{-150, 105}, {150, 145}}, textString = "%name"), Ellipse(visible = true, lineColor = {108, 88, 49}, fillColor = {139, 174, 255}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {100, 100}}), Text(visible = true, textColor = {108, 88, 49}, extent = {{-70, -70}, {70, 70}}, textString = "T-p")}));
    end T_to_p_sat;

    block Ramp_L "Generate ramp signal"
      parameter Real height = 1 "Height of ramps" annotation(Dialog(groupImage = "modelica://Modelica/Resources/Images/Blocks/Sources/Ramp.png"));
      parameter Modelica.Units.SI.Time duration(min = 0.0, start = 2) "Duration of ramp (= 0.0 gives a Step)";
      extends Modelica.Blocks.Interfaces.SignalSource;
    equation
      y = offset + (if time < startTime then 0 else if time < startTime + duration then (time - startTime) * height / duration else height);
      annotation(Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, initialScale = 0.1, grid = {10, 10}), graphics = {Line(visible = true, points = {{-80, 68}, {-80, -80}}, color = {192, 192, 192}), Polygon(visible = true, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid, points = {{-80, 90}, {-88, 68}, {-72, 68}, {-80, 90}}), Line(visible = true, points = {{-90, -70}, {82, -70}}, color = {192, 192, 192}), Polygon(visible = true, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid, points = {{90, -70}, {68, -62}, {68, -78}, {90, -70}}), Line(visible = true, points = {{-80, -70}, {-40, -70}, {31, 38}}), Text(visible = true, origin = {50, -20}, extent = {{-150, -150}, {150, -110}}, textString = "duration = %duration", fontName = "ISOCPEUR"), Line(visible = true, points = {{31, 38}, {86, 38}}), Text(visible = true, origin = {50, -80}, extent = {{-150, -150}, {150, -110}}, textString = "startTime = %startTime", fontName = "ISOCPEUR"), Text(visible = true, origin = {50, -52.165}, extent = {{-150, -150}, {150, -110}}, textString = "offset  = %offset", fontName = "ISOCPEUR"), Text(visible = true, origin = {50, 10}, extent = {{-150, -150}, {150, -110}}, textString = "height = %height", fontName = "ISOCPEUR")}), Documentation(info = "<html>
<p>
The Real output y is a ramp signal:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Blocks/Sources/Ramp.png\"
     alt=\"Ramp.png\">
</p>

<p>
If parameter duration is set to 0.0, the limiting case of a Step signal is achieved.
</p>
</html>"));
    end Ramp_L;
    annotation(Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {0, 114, 195}, fillColor = {255, 255, 255}, extent = {{-100, -100}, {100, 100}}, radius = 25), Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name")}));
  end Utilities;
  annotation(uses(Modelica(version = "4.0.0")), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, origin = {0, -30}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-79.157, -30}, {79.157, 30}}), Ellipse(visible = true, origin = {-58.981, -9.957}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-4.778, -3.849}, {4.778, 3.849}}), Ellipse(visible = true, origin = {-43.186, -16.24}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-5.177, -3.76}, {5.177, 3.76}}), Ellipse(visible = true, origin = {-20.621, -11.55}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-4.115, -3.584}, {4.115, 3.584}}), Ellipse(visible = true, origin = {1.413, -13.363}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-5.442, -3.363}, {5.442, 3.363}}), Ellipse(visible = true, origin = {27.694, -18.983}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-5.708, -4.38}, {5.708, 4.38}}), Ellipse(visible = true, origin = {60.347, -16.24}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-5.177, -3.76}, {5.177, 3.76}}), Ellipse(visible = true, origin = {-61.503, -39.026}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-5.708, -3.982}, {5.708, 3.982}}), Ellipse(visible = true, origin = {-33.629, -30.664}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-4.911, -3.053}, {4.911, 3.053}}), Ellipse(visible = true, origin = {-3.1, -31.195}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-6.239, -4.115}, {6.239, 4.115}}), Ellipse(visible = true, origin = {28.889, -38.893}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-6.637, -4.115}, {6.637, 4.115}}), Ellipse(visible = true, origin = {54.951, -31.593}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-5.794, -3.717}, {5.794, 3.717}}), Ellipse(visible = true, origin = {-20.488, -43.893}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-5.575, -3.893}, {5.575, 3.893}}), Ellipse(visible = true, origin = {-64.823, 6.104}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-3.98, -2.389}, {3.98, 2.389}}), Ellipse(visible = true, origin = {-51.017, 12.565}, rotation = -5.174, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-4.247, -2.565}, {4.247, 2.565}}), Ellipse(visible = true, origin = {-37.479, 1.989}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-6.636, -4.115}, {6.636, 4.115}}), Ellipse(visible = true, origin = {-19.692, 13.494}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-4.778, -3.494}, {4.778, 3.494}}), Ellipse(visible = true, origin = {-37.037, 22.476}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-3.806, -2.476}, {3.806, 2.476}}), Ellipse(visible = true, origin = {-0.976, 5.795}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-4.38, -4.205}, {4.38, 4.205}}), Ellipse(visible = true, origin = {-13.542, 28.538}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-2.699, -1.462}, {2.699, 1.462}}), Ellipse(visible = true, origin = {7.519, 17.034}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-5.442, -2.966}, {5.442, 2.966}}), Ellipse(visible = true, origin = {27.163, 6.724}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-6.239, -3.276}, {6.239, 3.276}}), Ellipse(visible = true, origin = {37.517, 23.272}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-6.239, -3.272}, {6.239, 3.272}}), Ellipse(visible = true, origin = {52.252, 9.289}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-6.904, -3.717}, {6.904, 3.717}}), Rectangle(visible = true, origin = {0, 25}, fillColor = {255, 255, 255}, extent = {{-79.157, -25}, {79.157, 25}}), Ellipse(visible = true, origin = {-62.874, 38.493}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-7.126, -3.982}, {7.126, 3.982}}), Ellipse(visible = true, origin = {-32.254, 36.194}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-4.911, -3.806}, {4.911, 3.806}}), Ellipse(visible = true, origin = {7.831, 36.061}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-7.831, -3.939}, {7.831, 3.939}}), Ellipse(visible = true, origin = {36.768, 39.024}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-6.768, -3.717}, {6.768, 3.717}}), Ellipse(visible = true, origin = {63.272, 33.45}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-6.728, -3.45}, {6.728, 3.45}}), Ellipse(visible = true, origin = {44.333, 0.133}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-7.964, -4.646}, {7.964, 4.646}}), Ellipse(visible = true, origin = {13.097, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-3.097, -2.389}, {3.097, 2.389}})}));
end TwoPhase;
