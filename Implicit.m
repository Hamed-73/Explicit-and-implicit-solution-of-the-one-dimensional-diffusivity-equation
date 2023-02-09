% IN THE NAME OF GOD
%% DOCUMENT TITLE: Implicit solution of the one-dimensional diffusion equation
% this program will be able to consider pressure distribution with implicit formulations
%%
clear 
close all
clc
clf
%% Initialize & Define Parameters
fi=0.2; % porosity
k=9.86923*10^(-13); % permeability (m^2)
N=10; 
PL=2; % pressure of left (atm)
Pr=1; % pressure of right (atm)
P0=1; % pressure initial (atm)
mio=0.001; % viscosity (pa*s)
ct=(10^(-4))./(1.013*10^(5)); % Compressibility (pa^-1)
L=0.1; % lenght(m)
tmax=0.2; % maximum time(sec)
deltat=.0000025; % Time step(sec)
t=0:deltat:tmax; 
%%
eta=k./(mio*fi*ct); % hydraulic diffusivity
x=linspace(0,L,10);
deltax=x(1,2)-x(1,1); % grid size
P=zeros(length(t),length(x));
%%
Pprior=ones(1,length(x)).*P0;
P(1,:)=Pprior;

for ii=2:length(t)

%%% create Coefficient matrix
h=(deltax.^2)./(eta.*deltat);
A=sparse(length(x),length(x));
for i=2:length(x)-1
A(i,i-1)=-h;
A(i,i)=(1+2*h);
A(i,i+1)=-h;
end
A(1,1)=1;
A(length(x),length(x))=1;

b=zeros(length(x),1);
b(2:length(x)-1)=Pprior(2:length(x)-1);
if ii==1
    
b(1)=P0;b(length(x))=P0;
else
b(1)=Pr;b(length(x))=PL; 

x=A\b;
P(ii,:)=x;
x=Pprior;

end
end

%% Final Plots
[T,X]=meshgrid(0:deltax:L,0:deltat:tmax);
Z=P;
surf(T,X,Z)
colorbar
xlabel('x(m)')
ylabel('t(second)')
zlabel('pressure(atm)')
title('Pressure distribution');

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
