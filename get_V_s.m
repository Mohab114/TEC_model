function [V_s J_s] = get_V_s(T_e, phi_e, phi_c, d)

    A = 1.202e+6;
    m_e = 9.10938356e-31;
    k_b_ev = 8.6173303e-5;
    k_b_si = 1.38064852e-23;
    e_e = 1.60217662e-19;
    e_0 = 8.85418782e-12;
    
    [z_c_table gamma_c_table z_e_table gamma_e_table] = gamma_zeta_table;
    
    a = 0.7388;
    b = 1.333;
    c = 2.198;
    dd = 0.6712;
    e = -1.476;

    J_s = A * T_e^2 * exp((-1 * phi_e)/(k_b_ev*T_e));

    z_c = (((2*pi*m_e*e_e^2)/(e_0^2 * k_b_si^3))^(1/4)) * ((J_s^(1/2) * d)/(T_e^(3/4)));

    if z_c <= 10
        gamma_c = interp1(z_c_table, gamma_c_table, z_c, 'spline');
    else
        gamma_c = a*z_c.^b + c + e.*z_c.^dd;
    end

    V_s = (-1) * (phi_e - phi_c - gamma_c * k_b_ev * T_e);
    



end