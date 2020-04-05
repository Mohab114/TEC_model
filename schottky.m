function [J_acc, bb] = schottky(T_e, phi_e, V_s, d, V)

A = 1.202e+6;
m_e = 9.10938356e-31;
k_b_ev = 8.6173303e-5;
k_b_si = 1.38064852e-23;
e_e = 1.60217662e-19;
e_0 = 8.85418782e-12;

bb = 60;

J_acc = A .* T_e^2 .* exp(-1.*phi_e ./ (k_b_ev*T_e)) .* exp((e_e .*(bb .* e_e .* (V - V_s) ./ (16*pi*e_0*d)).^(1/2)) ./ (k_b_si*T_e));

end

