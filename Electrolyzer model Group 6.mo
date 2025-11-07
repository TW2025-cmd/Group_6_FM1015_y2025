model ElectrolyzerStack
// Inputs
  parameter Real i_pp = 2.0;// A/cm2
  // Constants
  // From table 5
  constant Real Perm_O2 = 2.26e-14;// [mol/s m Pa], Membrane Oxygen permeability
  constant Real Perm_H2 = 5.31e-14;// [mol/s m Pa], Membrane Hydrogen permeability
  constant Real M_H2O = 18.0;// [g/mol], Molar mass H2O
  constant Real M_O2 = 32.0;// [g/mol], Molar mass O2
  constant Real M_H2 = 2.0;// [g/mol], Molar mass H2
  constant Real F = 9.648e4;// [C/mol],Faraday's constant
  constant Real v_e = 2.0;// hidden electron coeff.
  constant Real v_H2O = -1.0;// stoic coeff. H2O
  constant Real v_H2 = 1.0;// stoic coeff. H2
  constant Real v_O2 = 0.5;// stoic coeff. O2
  // Controller constants
  constant Real K_p_a = 1e3;// Controller gain (Anode)
  constant Real K_p_c = 20.0;// Controller gain (Cathode)
  constant Real u_v_a = 0.0;// Control signal anode
  constant Real u_v_c = 0.0;// Control signal cathode
  // Other constants
  constant Real Va = 30e-3;// [m3], converted from [L]
  constant Real Vc = 5e-3;// [m3], converted from [L]
  constant Real delta_m = 2.0e-4;// [m], converted from micro meter, membrane thickness
  constant Real T = 60 + 273.15;// [K], Temperature
  constant Real Nc = 34;// Number of electrolycers in stack
  constant Real A_m = 90;// [cm^2], superficial membrane area
  constant Real p_scale = 1e5;// [Pa]
  constant Real p_a_e = 1e5;// [Pa]
  constant Real p_c_e = 20e5;// [Pa]
  constant Real rho_scale = 1.0;//[kg/m3]
  constant Real md_v_a_c = 2.1604e-3;// [kg/s]
  constant Real md_v_c_c = 0.0434e-3;// [kg/s]
  // Initial Parameters
  parameter Real n_H2O_a0 = 555.56;// mol From problem 2
  parameter Real n_H2_a0 = 0.0;// mol
  parameter Real n_O2_a0 = 0.7224;// mol From problem 3
  parameter Real n_H2O_c0 = 92.59;// 92.59 mol From problem 2
  parameter Real n_H2_c0 = 3.586;// mol From problem 3
  parameter Real n_O2_c0 = 0.0;// mol
  // Differential Variables
  // Anode differential Variables
  Real n_H2O_a(start = n_H2O_a0);
  Real n_H2_a(start = n_H2_a0);
  Real n_O2_a(start = n_O2_a0);
  // Cathode differential Variables
  Real n_H2O_c(start = n_H2O_c0);
  Real n_H2_c(start = n_H2_c0);
  Real n_O2_c(start = n_O2_c0);
  // VT Flash outputs for anode
  Real V_molar_L_a, V_molar_g_a;// Gas and Liquid molar volume
  Real x_H2O_a, x_H2_a, x_O2_a;//
  Real y_H2O_a, y_H2_a, y_O2_a;//
  Real p_tot_a;// Total pressure
  // VT Flash outputs for cathode
  Real V_molar_L_c, V_molar_g_c;// Gas and Liquid molar volume
  Real x_H2O_c, x_H2_c, x_O2_c;//
  Real y_H2O_c, y_H2_c, y_O2_c;//
  Real p_tot_c;// Total pressure
  // Algebraic Variables
  Real V_l_a_ref, V_l_c_ref;
  Real nd_v_a_e, nd_v_c_e;
  Real nd_H2O_i, nd_H2O_i_ref, nd_H2O_ai, nd_H2O_ae;
  Real nd_O2_ai, nd_O2_ae;
  Real nd_H2_ai, nd_H2_ae;
  Real nd_H2O_ci, nd_H2O_ce;
  Real nd_O2_ci, nd_O2_ce;
  Real nd_H2_ci, nd_H2_ce;
  Real md_p_cr, md_p_cr_ref;
  Real nd_H2O_cr, nd_O2_cr, nd_H2_cr;
  Real nd_d, nd_d_c;
  Real nd_H2O_g, nd_O2_g, nd_H2_g;
  Real nd_H2O_d, nd_O2_d, nd_H2_d;
  Real nd_O2_delta, nd_H2_delta;
  Real md_v_a_e, md_v_c_e;
  Real p_O2_c, p_O2_a, p_H2_c, p_H2_a;
  Real del_p_a, del_p_c;
  Real Y_a, Y_c;
  Real rho_i_a, rho_i_c;
  Real V_L_c, V_g_c, V_L_a, V_g_a;
  Real ptot_a_scale, ptot_c_scale;
  Real Volume_fraction_a,Volume_fraction_c;

  function HenryConstant
    input Real T_Henry "Temperature [K]";
    input Real H_max "Maximum Henry constant [atm]";
    input Real inv_T_max "Inverse Tmax [1/K]";
    input Real T_c "Critical temperature [K]";
    output Real H "Henry constant at temperature T [atm]";
  protected
    Real denom;
    Real T_bar;
    Real exponent;
    Real H_bar;
  algorithm
    denom := inv_T_max - (1/T_c);// ref Eq. 57
    T_bar := (1/T_Henry - 1/T_c)/denom;//ref Eq. 57
    exponent := 1.142 - 2.846*T_bar + 2.486*(T_bar^2) - 0.9761*(T_bar^3) + 0.2001*(T_bar^4);//ref Eq. 58
    H_bar := 10^(-exponent);//ref Eq. 58
    H := H_max*H_bar;//ref Eq. 56
  end HenryConstant;

  function VTFlashFunction
    input Real V;
    input Real T_flash;
    input Real n_H2O;
    input Real n_H2;
    input Real n_O2;
    output Real V_l; 
    output Real V_g;
    output Real x_a_H2O;
    output Real x_a_H2;
    output Real x_a_O2;
    output Real y_a_H2O;
    output Real y_a_H2;
    output Real y_a_O2;
    output Real p_a_tot;
    output Real V_molar_g;
    output Real V_molar_L;
  protected
    Real n_H2O_g, n_H2_g, n_O2_g, n_g_temp, n_H2_g_temp, n_O2_g_temp;
    Real n_H2O_l, n_H2_l, n_O2_l;
    Real n_tot, n_tot_g, n_tot_l;
    Real p_H2O, p_H2O_sat, p_H2, p_O2, p_tot;
    Real phi_H2O_l_sat;
    Real f_H2O_l_sat;
    Real y_a_n_H2O, y_a_n_H2, y_a_n_O2;
    Real H2_Henry, O2_Henry;
    Real H2O_molar_volume;
    // Constants (From table 5)
    constant Real A_H2O = 5.11564;//[K], Antoine constant
    constant Real B_H2O = 1687.537;//[K], Antoine constant
    constant Real C_H2O = 230.17;//[K], Antoine constant
    constant Real H_H2_max = 7.54e4*101325;//[Pa]
    constant Real H_O2_max = 7.08e4*101325;//[Pa]
    constant Real T_max_H2 = 3.09e-3;//1/K, inv max temp
    constant Real T_c_H2O = 641.7;//[K], water critical temp
    constant Real T_max_O2 = 2.73e-3;//1/K	inv max temp
    constant Real R = 8.314;//[J/mol K]	Gas constant

  algorithm
  // Step 1
    p_H2O_sat := 1e5*10^(A_H2O - B_H2O/(C_H2O + T_flash - 273.15));//ref Eq.54
    p_H2O := p_H2O_sat;//ref Eq.34
  // Step 2
    n_H2_g_temp := n_H2;//ref Eq.36
    n_O2_g_temp := n_O2;//ref Eq.36
  // Step 3
    phi_H2O_l_sat := 1.0012 - 1.6e-3*exp(8.7*(T_flash - 273.15)/373.15);//ref Eq.53
    H2O_molar_volume := (2.23 - 3.332e-3*T_flash + 6.421e-6*T_flash^2)/1e5;//ref Eq.55
    f_H2O_l_sat := phi_H2O_l_sat*p_H2O_sat;//ref Eq.52
    n_H2O_l := (f_H2O_l_sat*V - n_H2O*R*T_flash)/(f_H2O_l_sat*H2O_molar_volume - R*T_flash);//ref Eq.37
    n_H2O_g := n_H2O - n_H2O_l;//ref Eq.38
  // Step 4
    V_l := n_H2O_l*H2O_molar_volume;//ref Eq.39
    V_g := V - V_l;//ref Eq.40
    // Step 5
    p_H2 := (n_H2_g_temp*R*T_flash)/V_g;//ref Eq.41
    p_O2 := (n_O2_g_temp*R*T_flash)/V_g;//ref Eq.41
    p_tot := p_H2O + p_H2 + p_O2;//ref Eq.42
  // Step 6
    n_tot := n_H2O + n_H2 + n_O2;
    n_g_temp := n_H2O_g + n_H2_g_temp + n_O2_g_temp;//ref Eq.43
    y_a_n_H2O := n_H2O_g/n_g_temp;//ref Eq.44
    y_a_n_H2 := n_H2_g_temp/n_g_temp;//ref Eq.44
    y_a_n_O2 := n_O2_g_temp/n_g_temp;//ref Eq.44
    H2_Henry := HenryConstant(T_flash, H_H2_max, T_max_H2, T_c_H2O);//ref Eq.45
    O2_Henry := HenryConstant(T_flash, H_O2_max, T_max_O2, T_c_H2O);//ref Eq.45
    n_H2_l := n_H2O_l*((p_tot)*y_a_n_H2/H2_Henry);//ref Eq.45
    n_O2_l := n_H2O_l*((p_tot)*y_a_n_O2/O2_Henry);//ref Eq.45
    n_tot_l := n_H2O_l + n_H2_l + n_O2_l;//ref Eq.45
    x_a_H2O := n_H2O_l/n_tot_l;//ref Eq.46
    x_a_H2 := n_H2_l/n_tot_l;//ref Eq.46
    x_a_O2 := n_O2_l/n_tot_l;//ref Eq.46
    // Step 7
    n_H2_g := n_H2 - n_H2_l;//ref Eq.47
    n_O2_g := n_O2 - n_O2_l;//ref Eq.47
    n_tot_g := n_tot - n_tot_l;//ref Eq.47
    y_a_H2O := n_H2O_g/n_tot_g;//ref Eq.48
    y_a_H2 := n_H2_g/n_tot_g;//ref Eq.48
    y_a_O2 := n_O2_g/n_tot_g;//ref Eq.48
    p_a_tot := p_H2O + (n_H2_g*R*T_flash)/V_g + (n_O2_g*R*T_flash)/V_g;//ref Eq.49
    V_molar_g := V_g/n_tot_g;//ref Eq.50
    V_molar_L := V_l/n_tot_l;//ref Eq.51
  end VTFlashFunction;

algorithm
// VTFlash for anode
  (V_L_a, V_g_a, x_H2O_a, x_H2_a, x_O2_a, y_H2O_a, y_H2_a, y_O2_a, p_tot_a,V_molar_g_a, V_molar_L_a) := VTFlashFunction(Va, T, n_H2O_a, n_H2_a, n_O2_a);
// VTFlash for cathode
  (V_L_c, V_g_c, x_H2O_c, x_H2_c, x_O2_c, y_H2O_c, y_H2_c, y_O2_c, p_tot_c,V_molar_g_c, V_molar_L_c) := VTFlashFunction(Vc, T, n_H2O_c, n_H2_c, n_O2_c);

equation

// Differential equations
  der(n_H2O_a) = nd_H2O_i + nd_H2O_cr + nd_H2O_ai - nd_H2O_ae; 
  der(n_O2_a) = nd_O2_cr + nd_O2_ai - nd_O2_ae; //
  der(n_H2_a) = nd_H2_cr + nd_H2_ai - nd_H2_ae;
  der(n_H2O_c) = nd_H2O_ci - nd_H2O_ce - nd_H2O_cr;
  der(n_O2_c) = nd_O2_ci - nd_O2_ce - nd_O2_cr;   
  der(n_H2_c) = nd_H2_ci - nd_H2_ce - nd_H2_cr;  // Partial pressures
  

equation
// Anode
  p_H2_a = p_tot_a*y_H2_a;
  p_O2_a = p_tot_a*y_O2_a;
// Cathode
  p_H2_c = p_tot_c*y_H2_c;
  p_O2_c = p_tot_c*y_O2_c;
//
//
  nd_O2_delta = Nc*Perm_O2*A_m*(p_O2_c - p_O2_a)/delta_m;
//ref Eq.79
  nd_H2_delta = Nc*Perm_H2*A_m*(p_H2_c - p_H2_a)/delta_m;
//ref Eq.80

// Membrane transport
  nd_d_c = (1.34e-2*T + 0.03)*(A_m*i_pp/F);//ref Eq.73/74
  nd_d = Nc*nd_d_c;//ref Eq.75
  nd_H2O_g = (Nc*(v_H2O/v_e)*(A_m*i_pp/F));//ref Eq.71
  nd_O2_g = Nc*(v_O2/v_e)*(A_m*i_pp/F);//ref Eq.72
  nd_H2_g = Nc*(v_H2/v_e)*(A_m*i_pp/F);//ref Eq.86
  nd_H2O_d = x_H2O_a*nd_d;//ref Eq.76
  nd_O2_d = x_O2_a*nd_d;//ref Eq.78
  nd_H2_d = x_H2_a*nd_d;//ref Eq.77
  
// Anode equations
  V_l_a_ref = Va/3;//ref formula Table 6
  del_p_a = max(0,p_tot_a - p_a_e);//ref Eq.81
  Y_a = 1 - min(1, 2/3*(del_p_a/p_tot_a));//ref Eq.81
  rho_i_a = ((y_H2O_a*M_H2O) + (y_O2_a*M_O2) + (y_H2_a*M_H2))/V_molar_g_a;//ref Eq.82
  md_v_a_e = u_v_a*Y_a*sqrt(((del_p_a/p_scale)*(rho_i_a/(rho_scale))))*md_v_a_c; //  //ref Eq.81
  nd_v_a_e = md_v_a_e/((y_H2O_a*M_H2O) + (y_O2_a*M_O2) + (y_H2_a*M_H2));//ref Eq.84
  nd_H2O_ae = y_H2O_a*nd_v_a_e;//ref Eq.83
  nd_O2_ae = y_O2_a*nd_v_a_e;//ref Eq.83
  nd_H2_ae = y_H2_a*nd_v_a_e;//ref Eq.83
  nd_H2O_ai = nd_H2O_g - nd_H2O_d;//
  nd_O2_ai = nd_O2_g - nd_O2_d - nd_O2_delta;//
  nd_H2_ai = nd_H2_delta - nd_H2_d;//
  nd_H2O_i_ref = nd_H2O_g; 
  nd_H2O_i = max(0,nd_H2O_i_ref + K_p_a*(V_l_a_ref - V_L_a)); //ref Eq.85
// Cathode equations
  V_l_c_ref = Vc/3;//ref formula Table 6
  del_p_c = max(0,p_tot_c - p_c_e);//ref Eq.87
  Y_c = 1 - min(1, 2/3*(del_p_c/p_tot_c));//ref Eq.87
  rho_i_c = ((y_H2O_c*M_H2O) + (y_O2_c*M_O2) + (y_H2_c*M_H2))/V_molar_g_c;//ref Eq.88 05.11
  md_v_c_e = u_v_c*Y_c*sqrt((del_p_c/p_scale)*(rho_i_c/rho_scale))*md_v_c_c;//ref Eq.87
  nd_v_c_e = md_v_c_e/((y_H2O_c*M_H2O) + (y_O2_c*M_O2) + (y_H2_c*M_H2));//ref Eq.90
  nd_H2O_ce = y_H2O_c*nd_v_c_e;//ref Eq.89
  nd_O2_ce = y_O2_c*nd_v_c_e;//ref Eq.89
  nd_H2_ce = y_H2_c*nd_v_c_e;//ref Eq.89
  nd_H2O_ci = -nd_H2O_d;
  nd_O2_ci = nd_O2_d + nd_O2_delta;
  nd_H2_ci = nd_H2_g + nd_H2_d + nd_H2_delta;
// Cathode recirculation
  md_p_cr_ref = 5.13; //ref Table 6 Problem 5
  md_p_cr = max(0,md_p_cr_ref + K_p_c*(V_L_c - V_l_c_ref)); //ref Eq.91
  nd_H2O_cr = md_p_cr*x_H2O_c/((x_H2O_c*M_H2O) + (x_O2_c*M_O2) + (x_H2_c*M_H2));//ref Eq.92
  nd_O2_cr = md_p_cr*x_O2_c/((x_H2O_c*M_H2O) + (x_O2_c*M_O2) + (x_H2_c*M_H2));//ref Eq.92
  nd_H2_cr = md_p_cr*x_H2_c/((x_H2O_c*M_H2O) + (x_O2_c*M_O2) + (x_H2_c*M_H2));//ref Eq.92

  
  ptot_c_scale = p_tot_c/1e5;
  ptot_a_scale = p_tot_a/1e5;
  Volume_fraction_a = V_L_a/Va;
  Volume_fraction_c = V_L_c/Vc;

end ElectrolyzerStack;
