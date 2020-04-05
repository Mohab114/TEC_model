A = 1.202e+6;
m_e = 9.10938356e-31;
k_b_ev = 8.6173303e-5;
k_b_si = 1.38064852e-23;
e_e = 1.60217662e-19;
e_0 = 8.85418782e-12;

phi_e = 2.95;
phi_c = 4.0;
phi_c = phi_c + 0.85;

T_e_grad = [1223.9,1247.2,1236.6,1247.4,1256.5,1252.9,1254.7,1240.6,1264,1275,1278.6,1298.8,1324.6,1341.4,1356.5,1386.6,1539,1574.1,1532.5,1517.7,1504.3,1486.6,1484.9,1479.1,1463.7,1473.4];
T_e_grad = flip(T_e_grad);
T_e_grad = T_e_grad./T_e_grad(1);
T_e_grad = T_e_grad .* 1700;

T_e_fit = [1700,1700,1700,1690,1680,1660,1645,1625,1600,1590,1575,1570,1545,1540,1530,1520,1515,1505,1500,1490,1485,1485,1480,1480,1475,1470];

T_e = 1475;
T_e_range = linspace(1700,1480,26);

nn = 25;
T_e_fun = (linspace(1,(T_e_range(end)/T_e_range(1))^nn,26)).^(1/nn);
T_e_fun = T_e_fun * 1700;

%r = 42e-6;
%S = pi * r^2;
S = 3e-9;

gaps = [50e-6:20e-6:550e-6];
%gaps = [0 10e-6 gaps(2:end)];
%gaps = [530e-6];



f = figure;
hold on;
for i = 1:length(gaps)
    d = gaps(i);
    T_e = T_e_fit(i);


    %%% Saturation Voltage
    [V_s J_s] = get_V_s(T_e, phi_e, phi_c, d);
    %V_sat_range = [V_s:0.01:V_s+1];
    V_sat_range = [V_s:0.01:100];
    %%%

    
    %%% Schottky effect
    [J_acc, bb] = schottky(T_e, phi_e, V_s, d, V_sat_range);
    V_acc_range = V_sat_range;
    %%%


    %%% Critical Voltage
    [V_c J_c] = get_V_c(T_e, phi_e, phi_c, d);
    %V_retard_range = [V_c-1:0.01:V_c-0.01];
    V_retard_range = [-5:0.01:V_c-0.01];
    %%%


    %%% Space Charge Voltages
    [V_sc J_sc] = get_V_sc(T_e, phi_e, phi_c, d, J_c, J_s);
    V_sc_range = V_sc;
    %%%


    J_r = A * T_e^2 * exp(-1* (phi_c + (-1)* V_retard_range)./(k_b_ev * T_e));
    J_sc = J_sc;
    J_s = (A * T_e^2 * exp(-1* phi_e/(k_b_ev * T_e))) * ones(1, length(V_sat_range));


    V_sat_range = []; %%%%%%%
    J_s = []; %%%%%%%%%%%%%%%
    V_total = [V_retard_range V_sc_range V_sat_range V_acc_range];
    J_total = [J_r J_sc J_s J_acc];
    I_total = J_total * S;

    %f = figure;
    %plot(V_total, log10(I_total), '-o', 'linewidth', 2);
    semilogy(V_total, I_total, 'linewidth', 2);

    xlabel('Output voltage (V)');
    ylabel('Output current (A)');
    set(gca,'fontweight','bold', 'fontsize', 14);
    grid ON;

    %savefig([num2str(d*1e6)]);
    %print([num2str(d*1e6)],'-dpng', '-r300');
    %close(f);

end
set(gca,'yscale','log');
dim = [.5 .15 .3 .3];
str = {['Dia = ' num2str(2*sqrt(S/pi)*1e+6) ' \mum'], ...
    ['Gaps = ' num2str(gaps(1)*1e+6) ' to ' num2str(gaps(end)*1e+6) ' \mum'], ...
    ['\phi_{E} = ' num2str(phi_e) ' eV'], ...
    ['\phi_{C} = ' num2str(phi_c) ' eV'], ...
    ['T_{E} = ' num2str(T_e) ' K'], ...
    ['V_{S} = ' num2str(V_s) ' V'], ...
    ['V_{C} = ' num2str(V_c) ' V'], ...
    ['\beta = ' num2str(bb)]};
ann = annotation('textbox',dim,'String',str,'FitBoxToText','on');
ann.Color = 'black';
ann.FontSize = 14;
ylim([1e-12 1e-2]);
xlim([-5 100]);










hold on;
f2 = openfig('exp_IV2.fig');
axesObjs = get(f2, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes

xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');

close(f2);

semilogy(xdata-10, ydata);


