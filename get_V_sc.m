function [V_sc J_sc] = get_V_sc(T_e, phi_e, phi_c, d, J_c, J_s)

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

    J_sc = [J_c:0.1:J_s];

    gamma_e = log(J_s./J_sc);

    z_e = a * exp(-1*b*gamma_e) + c;

    for i = 1:length(z_e)
        if z_e(i) <= 2
            z_e(i) = interp1(gamma_e_table, z_e_table, gamma_e(i), 'spline');
        end
    end

    x_0 = (((e_0^2 * k_b_si^3)/(2*pi*m_e*e_e^2))^(1/4)) * ((T_e^(3/4))./(J_sc.^(1/2)));       %It is the Debye length

    z_c = (d./x_0) - z_e;

    a = 0.7388;
    b = 1.333;
    c = 2.198;
    dd = 0.6712;
    e = -1.476;

    for i = 1:length(z_c)
        if z_c(i) <= 10
            gamma_c(i) = interp1(z_c_table, gamma_c_table, z_c(i), 'spline');
        else
            gamma_c(i) = a*z_c(i).^b + c + e.*z_c(i).^dd;
        end
    end

    V_sc = (-1) * (phi_e - phi_c + (k_b_ev * T_e * (gamma_e - gamma_c)));


end