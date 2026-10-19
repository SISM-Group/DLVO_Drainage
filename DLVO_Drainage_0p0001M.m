clc;
clear all;

n1=1.411;                           %oil refractive indices
n3=1.33;                            %water
yo=-62*10^(-3);                     %oil/water potential, 0.0001 M KCl, pH 5.8   V

% substrateName = 'Magnetite';
% n2 = 2.42;
% ys = -11*10^(-3);

% substrateName = 'Gold';
% n2 = 3;
% ys = -21*10^(-3);

substrateName = 'Silica';
n2 = 1.448;
ys = -88*10^(-3);


Kb=1.381*10^(-23);                  %Boltzmann constant                 j/kelvin
T=273.15+20;                        %Kelvin

% Parameters for distance-dependent Hamaker calculation
c_light = 2.998e8;
q_ret = 1.185;
hp = 1.055e-34;
Ve = 2.068e16;

e=80;                               %dielectric constant of the aqueous solution
ee=8.854*10^(-12);                  %vacuum permittivity                F/m or C2/(N.m2)
DL=30.3*10^(-9);                    %Debye length, 0.0001 M KCl, pH 5.8     m
k=1/DL;                             %Debye parameter                         1/m

h_start = 0.5;
h_end = 500;
h_step = 0.1;
h = h_start:h_step:h_end;

h_m_all = h*1e-9;

% ============================================================
% Distance-dependent Hamaker calculation
% ============================================================

landa1 = c_light/((pi/2)*Ve)*sqrt(2/(n3^2*(n1^2+n3^2)));
landa2 = c_light/((pi/2)*Ve)*sqrt(2/(n3^2*(n2^2+n3^2)));

I1 = (1 + (h_m_all./landa1).^q_ret).^(-1/q_ret);
I2 = (1 + (h_m_all./landa2).^q_ret).^(-1/q_ret);

A_non = 3*hp*Ve/(8*sqrt(2)) ...
    *(n1^2-n3^2)*(n2^2-n3^2)/(n1^2-n2^2) ...
    .*(I2./sqrt(n2^2+n3^2) - I1./sqrt(n1^2+n3^2));

A_z = 3*Kb*T/4*(1+2*k*h_m_all).*exp(-2*k*h_m_all);

A_t = A_z + A_non;

A_used = A_t;


% EDL parameters
qe = 1.602176634e-19;
z  = 1;
kappa = k;

psi1 = yo;
psi2 = ys;


%% ============================================================
% Pressure calculations
% ============================================================

h_m = h*1e-9;
kh  = kappa*h_m;
% ============================================================
% EDL pressure HHF
% ============================================================

Pi_HHF = e*ee*kappa^2 .* ...
    (2*psi1*psi2.*cosh(kh) - psi1^2 - psi2^2) ./ ...
    (2*sinh(kh).^2);     % Pa
% ============================================================
% EDL pressure SA
% ============================================================

Pi_SA = 2*e*ee*kappa^2*psi1*psi2.*exp(-kh);   % Pa

% ============================================================
% vdW pressure using A(h)
% ============================================================

dA_dh = gradient(A_used, h_m);

DPr_vdW = -(A_used ./ (6*pi*h_m.^3)) ...
          + (dA_dh ./ (12*pi*h_m.^2));   % Pa


% ============================================================
% Gregory compression approximation pressure
% ============================================================

N_A = 6.02214076e23;
C = 0.0001*1000;
ni = C*N_A;

Yeo = z*qe*yo/(Kb*T);
Yes = z*qe*ys/(Kb*T);

PedlComp = zeros(size(h));

for i = 1:length(h)

    kh_i = k*h(i)*1e-9;

    PedlComp(i) = ni*Kb*T * ...
        ( ...
        2*sqrt(1 + 0.25*(Yeo + Yes)^2 * csch(kh_i/2)^2) ...
        - ((Yeo - Yes)^2 * exp(-kh_i)) / ...
          (1 + 0.25*(Yeo + Yes)^2 * csch(kh_i/2)^2) ...
        - 2 ...
        );

end



DPr_dlvoDrainage = PedlComp + DPr_vdW;   % Pi_DLVO_used_for_drainage




%% ============================================================
% Thin-film drainage calculation
% Using Gregory EDL compression approximation + vdW pressure A(h)
% ============================================================

r_drop = 185e-6;
sigma = 51.5e-3;
rho_water = 998.2;
rho_dodecane = 749.5;
mu_water = 1.002e-3;
g = 9.81;

delta_rho = rho_water - rho_dodecane;

c_g = r_drop^2 * sqrt(2*delta_rho*g/(3*sigma));

F_b = (4/3)*pi*r_drop^3*delta_rho*g;
P_b = F_b/(pi*c_g^2);

h_drain_nm = fliplr(h);
Pi_drain = fliplr(DPr_dlvoDrainage);

P_net = P_b - Pi_drain;

idx_eq = find(P_net <= 0, 1, 'first');
equilibrium_found = false;

if isempty(idx_eq)

    h_calc_nm = h_drain_nm;

else

    h_eq_nm = interp1([P_net(idx_eq-1), P_net(idx_eq)], ...
                      [h_drain_nm(idx_eq-1), h_drain_nm(idx_eq)], ...
                      0, 'linear');

    h_calc_nm = h_drain_nm(1:idx_eq-1);
    equilibrium_found = true;

end

time_s = zeros(size(h_calc_nm));

for j = 1:length(h_calc_nm)-1

    h_mid_nm = (h_calc_nm(j) + h_calc_nm(j+1))/2;
    h_mid_m = h_mid_nm*1e-9;

    Pi_mid = interp1(h, DPr_dlvoDrainage, h_mid_nm, 'linear');

    P_net_mid = P_b - Pi_mid;

    V_drain = (2*h_mid_m^3/(3*mu_water*c_g^2)) * P_net_mid;

    dh_m = (h_calc_nm(j) - h_calc_nm(j+1))*1e-9;
    dt_s = dh_m/V_drain;

    time_s(j+1) = time_s(j) + dt_s;

end


figure;

plot(time_s, h_calc_nm, '-b', 'LineWidth', 1.5);
grid on;
hold on;

if equilibrium_found

    yline(h_eq_nm, '--k', sprintf('h_{eq} = %.2f nm', h_eq_nm), ...
        'LabelHorizontalAlignment', 'left', 'LineWidth', 1.1);

    title(['Dodecane-' substrateName '-0.0001 M KCl: predicted stable film thickness']);

    fprintf('Predicted equilibrium film thickness: %.4f nm\n', h_eq_nm);

else

    title(['Dodecane-' substrateName '-0.0001 M KCl: film thinning versus time']);

    fprintf('No stable equilibrium was found within %.1f to %.1f nm.\n', min(h), max(h));
    fprintf('The curve stops at the lower Pi(h) calculation limit (%.1f nm); this is not assigned as rupture.\n', min(h));

end

xlabel('Time (s)');
ylabel('Film thickness, h (nm)');

fprintf('Film radius, c_g = %.4f micrometre\n', c_g*1e6);
fprintf('Buoyancy driving pressure, P_b = %.4f Pa\n', P_b);

%% ============================================================
% Export Excel file
% ============================================================

conditionName = '0p0001M';

output_excel_file = ['DLVO_Drainage_dodecane_' substrateName '_' conditionName '.xlsx'];

% ------- Pressure sheet --------

DPr_DLVO_HHF = Pi_HHF + DPr_vdW;
DPr_DLVO_SA  = Pi_SA + DPr_vdW;
DPr_DLVO_CA  = PedlComp + DPr_vdW;

T_pressure = table( ...
    h(:), ...
    Pi_HHF(:), ...
    Pi_SA(:), ...
    PedlComp(:), ...
    DPr_vdW(:), ...
    DPr_DLVO_HHF(:), ...
    DPr_DLVO_SA(:), ...
    DPr_DLVO_CA(:), ...
    DPr_dlvoDrainage(:), ...
    'VariableNames', { ...
    'h_nm', ...
    'Pi_EDL_HHF', ...
    'Pi_EDL_SA', ...
    'Pi_EDL_CA', ...
    'Pi_vdW', ...
    'Pi_DLVO_HHF', ...
    'Pi_DLVO_SA', ...
    'Pi_DLVO_CA', ...
    'Pi_DLVO_used_for_drainage'} ...
    );

writetable(T_pressure, output_excel_file, 'Sheet', 'Pressure');


% -------- Drainage sheet: final thickness vs time  -------

T_drainage = table( ...
    time_s(:), ...
    h_calc_nm(:), ...
    'VariableNames', { ...
    'Time_s', ...
    'FilmThickness_nm'} ...
    );

writetable(T_drainage, output_excel_file, 'Sheet', 'Drainage');

fprintf('Pressure and drainage data exported to: %s\n', output_excel_file);