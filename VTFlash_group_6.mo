model VTFlash
	// Inputs 
	parameter Real V = 30e-3 "[L], volume in anode seperator";
	parameter Real T = 60 + 273.15 "[◦C], temperature in anode seperator";
	parameter Real n_H2O = 556;
	parameter Real n_H2 = 0;
	parameter Real n_O2 = 0.7224;
	
	// Variables
      Real n_H2O_g, n_H2_g, n_O2_g, n_g_temp, n_H2_g_temp, n_O2_g_temp;
	  Real n_H2O_l, n_H2_l, n_O2_l;
      Real n_tot, n_tot_g, n_tot_l;
      Real V_l, V_g;
	  Real V_molar_l, V_molar_g;
      Real p_H2O, p_H2O_sat, p_H2, p_O2, p_tot, p_a_tot;
	  Real phi_H2O_l_sat;
      Real f_H2O_l_sat;
      Real y_a_n_H2O, y_a_n_H2, y_a_n_O2; 
      Real y_a_H2O, y_a_H2, y_a_O2;
      Real x_a_H2O, x_a_H2, x_a_O2;
      Real H2_Henry, O2_Henry; 
	  Real H2O_molar_volume;
    
	// Constants
	constant Real A_H2O = 5.11564		"[K], Antoine correlation, valid [273.20, 473.20]K x [10^-2, 16]bar";
	constant Real B_H2O = 1687.537		"[K], Antoine correlation, valid [273.20, 473.20]K x [10^-2, 16]bar";
	constant Real C_H2O = 230.17		"[K], Antoine correlation, valid [273.20, 473.20]K x [10^-2, 16]bar";
	constant Real H_H2_max = 7.54e4 	"[atm], Hydrogen Maximum Henry Constant";
	constant Real H_O2_max = 7.08e4 	"[atm], Hydrogen Maximum Henry Constant";
	constant Real T_max_H2 = 3.09e-3	"[1/K), Hydrogen inverse maximum temperature"; 
	constant Real T_c_H2 = 33.2	   		"[K], Hydrogen critical temperature";
	constant Real T_max_O2 = 2.73e-3	"[1/K], Oxygen inverse maximum temperature"; 
	constant Real T_c_O2 = 154.6	   	"[K], Oxygen critical temperature";
	constant Real R = 8.314 			"[J/ mol K], ideal gas constant";
	constant Real F = 9.648e4			"[C/mol], Faraday's constant";
	constant Real Pa_to_atm = 101325;
	constant Real T_c_H2O = 641.7;
	
	// Functions
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
			denom := inv_T_max - (1/T_c);
			T_bar := (1/T_Henry - 1/T_c) /denom;
			exponent := 1.142 - 2.846 * T_bar + 2.486 * (T_bar^2) - 0.9761 * (T_bar^3) + 0.2001 * (T_bar^4);
			H_bar := 10^(-exponent);
			H := H_max * H_bar;
	end HenryConstant;

equation

	//Step 1, Approximate partial water pressure
	p_H2O_sat = 1e5*10^(A_H2O - B_H2O/(C_H2O + T - 273.15)); //[Pa],Re-order of Equation 54, to get p_sat_H2O(T) alone with regards to p_sigma being 1
	p_H2O = p_H2O_sat;
    
	//Step 2, assumptions
	n_H2_g_temp = n_H2;
	n_O2_g_temp = n_O2;

	//Step 3, Water in liquid and gas phase
	phi_H2O_l_sat = 1.0012 - 1.6e-3 * exp(8.7 * (T-273.15) /373.15);
	H2O_molar_volume = (2.23 - 3.332e-3 * T + 6.421e-6 * T^2) / 1e5;
	f_H2O_l_sat = phi_H2O_l_sat * p_H2O_sat;
	n_H2O_l = (f_H2O_l_sat * V - n_H2O*R*T)/(f_H2O_l_sat*H2O_molar_volume- R*T);
	n_H2O_g = n_H2O - n_H2O_l;

	//Step 4, Phase volumes"
	V_l = n_H2O_l *H2O_molar_volume;
	V_g = V - V_l;

	//Step 5, Partial Pressure and total pressure
	p_H2 = (n_H2_g_temp * R * T) / V_g;
	p_O2 = (n_O2_g_temp * R * T) / V_g;
    p_tot = p_H2O +  p_H2 +  p_O2;
    
    
	//Step 6, Temporary phase amounts"
	n_tot = n_H2O + n_H2 + n_O2 "moles";
	n_g_temp = n_H2O_g + n_H2 + n_O2 "moles";
	y_a_n_H2O = n_H2O_g/n_g_temp"mole fraction";
	y_a_n_H2 = n_H2_g_temp/n_g_temp"mole fraction";
	y_a_n_O2 = n_O2_g_temp/n_g_temp"mole fraction";
	H2_Henry = HenryConstant(T, H_H2_max, T_max_H2, T_c_H2O);
	O2_Henry = HenryConstant(T, H_O2_max, T_max_O2, T_c_H2O); 
	n_H2_l = n_H2O_l * ((p_tot/Pa_to_atm) * y_a_n_H2/H2_Henry);
	n_O2_l = n_H2O_l * ((p_tot/Pa_to_atm) * y_a_n_O2/O2_Henry);
    n_tot_l = n_H2O_l + n_H2_l + n_O2_l;
	x_a_H2O = n_H2O_l / (n_tot_l);
	x_a_H2 =  n_H2_l /  (n_tot_l);
	x_a_O2 = n_O2_l / (n_tot_l);

	// Step 7, correction"
	n_H2_g = n_H2 - n_H2_l;
	n_O2_g = n_O2 - n_O2_l;
	n_tot_g = n_tot - n_tot_l;
	y_a_H2O = n_H2O_g/n_tot_g;
	y_a_H2 = n_H2_g/n_tot_g;
	y_a_O2 = n_O2_g/n_tot_g;
	p_a_tot = p_H2O + (n_H2_g * R * T) / V_g + (n_O2_g * R * T) / V_g;
	V_molar_g = V_g /n_tot_g;
	V_molar_l = V_l / n_tot_l;
end VTFlash;
