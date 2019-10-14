# Cholera_Model_code
% syms S I R B delta1 gamma1 mu1 betaH1 betaL1 N1 r1 xi1 sigma1 kappa1
% 
% 
% X(1,1) = mu1*N1 - betaH1*S*I - betaL1*(S*B/(kappa1 + B)) - mu1*S +sigma1*R;
% 
% X(1,2) = betaH1*S*I + betaL1*(S*B/(kappa1 + B)) - (gamma1 + mu1)*I;
% 
% X(1,3) = gamma1*I - (mu1 + sigma1)*R;
% 
% X(1,4) = xi1*I - delta1*B;
% 
% var  = {mu1,            betaL1,             gamma1,   delta1,   N1,      betaH1,               xi1,      sigma1,        kappa1 };
% dat1 = [1/(43.5*365),   0.011/(10^3),       1/5,      1/30,     12347,   (1.57)*(10^(-3)),      10,      10^(-3),      10^(-3) ];
% dat2 = [1/43.5,         0.011*365/(10^5),   365/5,    365/30,   12347,   (1.57)*365*(10^(-3)),  ];
% dat3 = [1/43.5,         0.011*365/(10^5),   365/5,    365/30,   12347,   (1.57)*365*(10^(-3)),  ];
% dat4 = [1/(43.5*365),   0.011/(10^6),       1/5,      1/30,     12347,   (1.57)*(10^(-5)),      ];
% dat5 = [1/(43.5*365),   0.011/(10^6),       1/4,      1/30,     12347,   (1.57)*(10^(-5)),      ];
% 
% dat = dat1;
% 
% eqn1 = subs(X,var,dat) == 0;
% range1 = [10^(-8) Inf; 10^(-8) Inf; 10^(-8) Inf; 10^(-8) Inf;]; % excludes boundary sols. only EE %
% range2 = [0 Inf; 0 Inf; 0 Inf; 0 Inf;]; % includes boundary solutions %
% EE = vpasolve(eqn1, [S I R B], range1);
% 
% Se = EE.S; Ie = EE.I; Re = EE.R; Be = EE.B;
% 
% format long
% 
% EEtemp = [Se Ie Re Be];
% EE1 = double(EEtemp)
% % R_0 = vpa(subs((N1*betaL1*xi1/(gamma1 + mu1)*()*(betaH1 + ((betaL1*10)/delta1)), var, dat))
% 
% global  delta gamma mu betaH betaL N xi sigma kappa;
% datcell = num2cell(dat);
% [mu, betaL, gamma, delta, N, betaH, xi, sigma, kappa] = datcell{:};

% y0 = [EE1(1,1) EE1(1,2) EE1(1,3) EE1(1,4)];



x0 = [12347-100 100 0 2000]; % Setting initial conditions]
tspan = 0:.1:200;


[t,y] = ode45(@G,tspan,x0);


figure
subplot(2,2,1)
plot(t,y(:,1),'b')
subplot(2,2,3)
plot(t,y(:,3),'g')
subplot(2,2,4)
plot(t,y(:,4),'c')
subplot(2,2,2)
plot(t,y(:,2),'r')


function dXdt = G(~,y)

global  delta gamma mu betaH betaL N xi sigma kappa
          
dXdt = zeros(4,1);

dXdt(1) = mu*N - betaH*y(1)*y(2) - betaL*(y(1)*y(4)/(kappa + y(4))) - mu*y(1) + sigma*y(3);
dXdt(2) = betaH*y(1)*y(2) + betaL*(y(1)*y(4)/(kappa + y(4))) - (gamma + mu)*y(2);
dXdt(3) = gamma*y(2) - (mu + sigma)*y(3);
dXdt(4) = xi*y(2) - delta*y(4);

end

