tspan = [0:1:7200];
%y1=5.636457E-3; %mol P1/RT
%y5=1.252546E-3; %mol P6/RT
pch4=0.45*101325;
par=pch4;
pna=0.1*101325;
N=6.023E17;
v0=0.001; %m^3
y1=pch4*v0/(8.314*973); %mol P1V0/RT
y5=pna*v0/(8.314*973); %mol P6V0/RT
yar=par*v0/(8.314*973); %mol ParV0/RT
n0=y1+y5+yar;
s=[y1,0,0,0,y5,0,0,0,0,0,0,0,0,0,0,0];
[T Y] = ode15s(@(t,a)drcrates(t,a),tspan,s);
%for volume calculations
nt1 = sum(Y,2)+yar;
for i=1:length(T)
    v(i)=v0*nt1(i)./n0;% Calculate volume over time
end
kf=[1.1*10^-8 2.2*10^9 1.5*10^11 N*2.2*10^-21 N*3.2*10^-22 N*10^-23 N*9.4*10^-27 N*2.3*10^-23 N*3.1*10^-21 N*6.9*10^-13 N*9.9*10^-24 9.5*10^9 8.7*10^10 6.8*10^3 8.4*10^5 7.1*10^3 5*10^9 5.8*10^8 4.4*10^-1 1.4*10^8 2.4*10^1 4.7*10^6 N*1.9*10^-10 N*2.8*10^-10];

kb=[N*4.1*10^-10 N*1.5*10^-9  N*4.4*10^-9 N*1.9*10^-10 4.1*10^9 N*1.9*10^-12 N*6.6*10^-11 3.9*10^2 4.7*10^4 N*5.7*10^-15 N*3.8*10^-19 N*3.3*10^-9 N*7.3*10^-9 N*3*10^-9  N*2.3*10^-10 N*1.4*10^-9 N*1.5*10^-9 N*3.1*10^-9 N*3.7*10^-12 N*5.1*10^-10 N*5*10^-10 N*3.9*10^-11 1.8*10^-10 2.7*10^-4];
v=v'

b = v.*(-kf(1)*(Y(:,1)./v)+kb(1)*(Y(:,2)./v).*(Y(:,3)./v)-kf(4)*(Y(:,1)./v).*(Y(:,5)./v)+kb(4)*(Y(:,7)./v).*(Y(:,2)./v)-kf(5)*(Y(:,5)./v).*(Y(:,1)./v)+kb(5)*(Y(:,8)./v)-kf(6)*(Y(:,1)./v).*(Y(:,4)./v)+kb(6)*(Y(:,9)./v).*(Y(:,2)./v)-kf(7)*(Y(:,1)./v).*(Y(:,4)./v)+kb(7)*(Y(:,10)./v).*(Y(:,3)./v)-kf(8)*(Y(:,1)./v).*(Y(:,4)./v)+kb(8)*(Y(:,11)./v)-kf(9)*(Y(:,1)./v).*(Y(:,6)./v)+kb(9)*(Y(:,12)./v)-kf(10)*(Y(:,1)./v).*(Y(:,3)./v)+kb(10)*(Y(:,2)./v).*(Y(:,13)./v)-kf(11)*(Y(:,1)./v).*(Y(:,2)./v)+kb(11)*(Y(:,14)./v).*(Y(:,3)./v));

kfn = 1.01*kf;             % Updated reaction rate constants
keq = kf ./ kb;       % Equilibrium constants
kbn = kfn ./ keq;     % Modified reaction rate constants
%% correct one

x1 = v.*((-kfn(1)*(Y(:,1)./v)./b) + (kbn(1)*((Y(:,2)./v).*(Y(:,3)./v))./b));

x4  = v.*((-kfn(4)*(Y(:,1)./v).*(Y(:,5)./v)./b) + (kbn(4)*(Y(:,7)./v).*(Y(:,2)./v)./b));

x5  = v.*((-kfn(5)*(Y(:,5)./v).*(Y(:,1)./v)./b) + (kbn(5)*(Y(:,8)./v)./b)); 
 
x6 = v.*((-kfn(6)*(Y(:,1)./v).*(Y(:,4)./v)./b)+(kbn(6)*(Y(:,9)./v).*(Y(:,2)./v)./b));
 
x7  = v.*((-kfn(7)*(Y(:,1)./v).*(Y(:,4)./v)./b)+(kbn(7)*(Y(:,10)./v).*(Y(:,3)./v)./b));
 
x8  = v.*((-kfn(8)*(Y(:,1)./v).*(Y(:,4)./v)./b)+(kbn(8)*(Y(:,11)./v)./b));
 
x9  = v.*((-kfn(9)*(Y(:,1)./v).*(Y(:,6)./v)./b)+(kbn(9)*(Y(:,12)./v)./b));
 
x10 = v.*((-kfn(10)*(Y(:,1)./v).*(Y(:,3)./v)./b)+(kbn(10)*(Y(:,2)./v).*(Y(:,13)./v)./b));
 
x11 = v.*((-kfn(11)*(Y(:,1)./v).*(Y(:,2)./v)./b)+(kbn(11)*(Y(:,14)./v).*(Y(:,3)./v)./b));

% for i=1:length(T)
%      x1 = v(i).*((-kfn(1)*(Y(i,1)./v(i))/b(i)) + (kbn(1)*(Y(i,2)./v(i))*(Y(i,3)/v(i))/b(i)));
% 
%      x4(i)  = v(i).*((-kfn(4)*(Y(i,1)./v(i))*(Y(i,6)/v(i))/b(i)) + (kbn(4)*(Y(i,8)./v(i))*(Y(i,2)/v(i))/b(i)));
% 
%      x5(i)  = v(i).*((-kfn(5)*(Y(i,6)./v(i))*(Y(i,1)/v(i))/b(i)) + (kbn(5)*(Y(i,11)./v(i))/b(i))); 
% 
%      x6(i)  = v(i).*((-kfn(6)*(Y(i,1)./v(i))*(Y(i,5)/v(i))/b(i))+(kbn(6)*(Y(i,10)./v(i))*(Y(i,2)/v(i))/b(i)));
% 
%      x7(i)  = v(i).*((-kfn(7)*(Y(i,1)./v(i))*(Y(i,5)/v(i))/b(i))+(kbn(7)*(Y(i,16)/v(i))*(Y(i,3)./v(i))/b(i)));
% 
%      x8(i)  = v(i).*((-kfn(8)*(Y(i,1)./v(i))*(Y(i,5)/v(i))/b(i))+(kbn(8)*(Y(i,12)./v(i))/b(i)));
% 
%      x9(i)  = v(i).*((-kfn(9)*(Y(i,1)./v(i))*(Y(i,7)/v(i))/b(i))+(kbn(9)*(Y(i,13)./v(i))/b(i)));
% 
%      x10(i) = v(i).*((-kfn(10)*(Y(i,1)/v(i))*(Y(i,3)/v(i))/b(i))+(kbn(10)*(Y(i,2)./v(i))*(Y(i,4)/v(i))/b(i)));
% 
%      x11(i) = v(i).*((-kfn(11)*(Y(i,1)/v(i))*(Y(i,2)./v(i))/b(i))+(kbn(11)*(Y(i,14)./v(i))*(Y(i,3)/v(i))/b(i)));
% 
%  end
%plot(T,x1);
%% % Plot results
figure
plot(T,Y(:,1),'-')
title('solution using ode15s');
xlabel('time_t');
ylabel('Mol');
figure
plot(T(1:10),x1(1:10),'DisplayName','x1','Color','red')
hold on
plot(T(1:10),x4(1:10),'DisplayName','x4','Color','blue')
plot(T(1:10),x5(1:10),'DisplayName','x5','Color','black')
plot(T(1:10),x6(1:10),'DisplayName','x6','Color','green')
plot(T(1:10),x7(1:10),'DisplayName','x7','Color','magenta')
plot(T(1:10),x8(1:10),'DisplayName','x8','Color',"#A2142F")
plot(T(1:10),x9(1:10),'DisplayName','x9','Color','yellow')
plot(T(1:10),x10(1:10),'DisplayName','x10','Color',"#EDB120")
plot(T(1:10),x11(1:10),'DisplayName','x11','Color','cyan')
title('X\(_i\) VS Time','Interpreter','latex')
xlabel('Time (sec)','Interpreter','latex')
ylabel('X\(_i\)',"Interpreter","latex")
legend("Interpreter",'latex');