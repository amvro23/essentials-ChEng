
% Define constant values
Cto = 0.1;
To = 373;
Ta = 373;
Ua = 4000;
dH_ra = -20000;
dH_r2a = -60000;
Cpa = 90;
Cpb = 90;
Cpc = 180;
A1 = 10;
Ea1 = 33256;
A2 = 0.09;
Ea2 = 74826;

% Discretize the volume
V = 1;

% Independent variable
tspan = [0, V];
t = linspace(tspan(1), tspan(2), 100);

% Initial conditions for the dependent variables (4)
y0 = [100; 0.1; 0.1; 450];

[v, y] = ode45(@(v,y) f_ODE(v, y, Cto, To, Ta, Ua, A1, Ea1, A2, Ea2, dH_ra, dH_r2a, Cpa, Cpb, Cpc), tspan, y0);

% Extract the solution for Fa, Fb, Fc, Fd
F = y(:, 1:3);
T = y(:, end);

figure;
% First subplot: Concentrations
subplot(2, 1, 1);
hold on;
plot(v, F(:,1), 'r-', 'LineWidth', 1.5);
plot(v, F(:,2), 'g-', 'LineWidth', 1.5);
plot(v, F(:,3), 'k-', 'LineWidth', 1.5);
hold off;
xlabel('Volume (V)');
ylabel('Concentration (Fi)');
legend('Fa', 'Fb', 'Fc');
grid on;

% Second subplot: Temperature
subplot(2, 1, 2);
plot(v, T, 'b-', 'LineWidth', 1.5);
xlabel('Volume (V)');
ylabel('Temperature (T)');
legend('T');
grid on;


function dydt = f_ODE(~, y, Cto, To, Ta, Ua, A1, Ea1, A2, Ea2, dH_ra, dH_r2a, Cpa, Cpb, Cpc)
    F = y(1:3);
    T = y(4);

    C = Cto * (F / sum(F) * (To / T));

    % Rate laws or reaction rate equations
    R = 8.314;
    r1a = -A1 * exp((Ea1 / R) * (1/300 - 1/T)) * C(1);
    r2a = -A2 * exp((Ea2 / R) * (1/300 - 1/T)) * C(1)^2;

    % Net rates or rates of species
    ra = r1a + r2a;
    rb = -r1a;
    rc = -r2a / 2;
    
    dF = [ra; rb; rc];

    dT = (Ua * (Ta - T) + (r1a * dH_ra + r2a * dH_r2a)) / (F(1) * Cpa + F(2) * Cpb + F(3) * Cpc);

    dydt = [dF; dT];
end

