
Fao = 0.0376;
Cpa = 163;
delCp = -9;
Cao = 18.8;
Cpc = 34.5;
mc = 0.111;
T0 = 1035;
X0 = 0;
Ta0 = 1250;
Ua = 0;

Vspan = [0 0.005];
y0 = [X0; T0; Ta0];
[v, y] = ode45(@(v,y) f_ODE(v, y, Fao, Cpa, delCp, Cao, Cpc, mc, T0, Ua), Vspan, y0);

% Extract the solution for Fa, Fb, Fc, Fd
X = y(:, 1);
T = y(:, 2);
Ta = y(:, 3);

k = 3.58*exp(34222*(1/T0-1./T));

% when I multiply or divide Nx1 lists
ra = -Cao * k .* (1-X) .* (T0./T) ./ (1+X);

figure;
% First subplot: Fa and T
subplot(2, 1, 1);
yyaxis left
plot(v, X, 'b-');
ylabel('Concentration (F)');
grid on;

yyaxis right
plot(v, T, 'g-');
ylabel('Temperature (T)');
grid on;

legend('Fa', 'T');
xlabel('Volume (V)');

% Second subplot: -ra
subplot(2, 1, 2);
plot(v, -ra, 'm-');
legend('-ra');
xlabel('Volume (V)');
ylabel('-ra');
grid on;


function dYdV = f_ODE(~, y, Fao, Cpa, delCp, Cao, Cpc, mc, T0, Ua)
  X = y(1); 
  T = y(2); 
  Ta = y(3);
  
  k = 3.58*exp(34222*(1/T0-1/T));
  ra = -Cao*k*(1-X)*(T0/T)/(1+X);
  deltaH = 80770 + delCp*(T-298);
  
  dXdV = -ra/Fao;
  dTdV = (Ua*(Ta-T)+ra*deltaH)/(Fao*(Cpa+X*delCp));
  dTadV = Ua*(T-Ta)/mc/Cpc;
  dYdV = [dXdV; dTdV; dTadV];
end