%% DOCUMENT TITLE: Explicit solution of the one-dimensional diffusion equation
% this program will be able to consider pressure distribution with explicit formulations
%%
clear
close all
clc

%%% Initialize & Define Parameters
fi=0.2; % porosity
k=9.86923*10^(-13); % permeability (m^2)
N=10; % number of layer
PL=2; % pressure of left (atm)
Pr=1; % pressure of right (atm)
P0=1; % pressure initial (atm)
mio=0.001; % viscosity (pa*s)
ct=(10^(-4))./(1.013*10^(5)); % Compressibility (pa^-1)
L=0.1; % lenght(m)
tmax=0.2; % maximum time (sec)
deltat=.0000025; % Time step (sec)
t=0:deltat:tmax;

eta=k./(mio*fi*ct); % hydraulic diffusivity
x=linspace(0,L,10);
deltax=x(1,2)-x(1,1); % grid size
P=zeros(length(t),length(x));

%%% Initial Condition & Boundary Condition
for n=1
    for i=1:length(x)
        
   P(n,i)=P0;
        
    end
end

for n=1:length(t)-1
    
    for i=2:length(x)-1
        
P(n+1,1)=Pr;

P(n+1,end)=PL;

P(n+1,i)=P(n,i)+(deltat.*eta./(deltax.^2)).*(P(n,i+1)+P(n,i-1)-2.*P(n,i));

    end
end

%%% plot section
[T,X]=meshgrid(0:deltax:L,0:deltat:tmax);
Z=P;

figure 
mesh(T,X,Z)
colorbar
xlabel('x(m)')
ylabel('t(second)')
zlabel('pressure(atm)')
title('Pressure distribution');

figure
plot(X,Z)
xlabel('t(second)')
ylabel('pressure(atm)')
title('Pressure distribution');

figure
plot(T,Z)
xlabel('x(m)')
ylabel('pressure(atm)')
title('Pressure distribution');
grid on

figure
plot3(T,X,Z)
xlabel('x(m)')
ylabel('t(second)')
zlabel('pressure(atm)')
title('Pressure distribution');
grid on

figure
surf(T,X,Z)
colorbar
xlabel('x(m)')
ylabel('t(second)')
zlabel('pressure(atm)')
title('Pressure distribution');
