function [V_c J_c] = get_V_c(T_e, phi_e, phi_c, d)

    A = 1.202e+6;
    m_e = 9.10938356e-31;
    k_b_ev = 8.6173303e-5;
    k_b_si = 1.38064852e-23;
    e_e = 1.60217662e-19;
    e_0 = 8.85418782e-12;

    [z_c_table gamma_c_table z_e_table gamma_e_table] = gamma_zeta_table;

    a = -1.497;
    b = 0.5125;
    c = 2.554;

    z_e = 2.554;

    J_c = sqrt(((z_e*T_e^(3/4))/d) * ((e_0^2 * k_b_si^3)/(2*pi*m_e*e_e^2))^(1/4));

    V_c = (-1) * ((-1 * ((k_b_ev*T_e))*log((J_c)/(A*T_e^2))) - phi_c);

    J_sat = A * T_e^2 * exp((-1 * phi_e)/(k_b_ev*T_e));

    J = [1e-6*J_c:0.01:1e3*J_c];

    gamma_e = log(J_sat./J);

    z_e = a * exp(-1*b*gamma_e) + c;

    for i = 1:length(z_e)
        if z_e(i) <= 2
            z_e(i) = interp1(gamma_e_table, z_e_table, gamma_e(i), 'spline');
        end
    end

    x_0 = (((e_0^2 * k_b_si^3)/(2*pi*m_e*e_e^2))^(1/4)) * ((T_e^(3/4))./(J.^(1/2)));       %It is the Debye length

    z_c = (d./x_0) - z_e;

    J_c_2 = J(find(z_c>0,1)); %%%%%
    
    if isempty(J_c_2)       %%%%%
        nothing = 1;        %%%%%
        J_c = J_sat;        %%%%%
    else                    %%%%%
        J_c = J_c_2;        %%%%%
    end                     %%%%%    
    
    if J_c > J_sat 
        J_c = J_sat;
    end

    V_c = (-1) *  ((-1 * ((k_b_ev*T_e))*log((J_c)/(A*T_e^2))) - phi_c);



end