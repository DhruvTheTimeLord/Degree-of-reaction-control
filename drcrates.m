function [dydt,x]  = drcrates(t,a,x_matrix)
% drcrates - Calculate reaction rates and degree of rate control
%
% Inputs:
%   t - time
%   a - concentration of species
%   x_matrix - matrix for the degree of rate control (optional)
%
% Outputs:
%   dydt - array of reaction rates
%   x - array of the degree of rate control
%

%a is concentration of species and t is time
% b = d[CH4]/dt , a1=[CH4]
% c = d[CH3]/dt , a2=[CH3]
% d = d[H]/dt , a3=[H]
% e = d[H2]/dt , a4=[H2]
% f = d[(na)2]/dt , a5=[(na)2]
% g = d[(na)]/dt , a6=[(na)]
% h = d[(na)3]/dt , a7=[(na)3]
% i = d[(na)H]/dt , a8=[(na)H]
% j = d[(na)3H]/dt , a9=[(na)3H]
% k = d[(na)2H]/dt , a10=[(na)2H]
% l = d[H(na)CH3]/dt , a11=[H(na)CH3]
% m = d[H(na)2CH3]/dt , a12=[H(na)2CH3]
% n = d[H(na)3CH3]/dt , a13=[H(na)3CH3]
% o = d[C2H6]/dt , a14=[C2H6]
% p = d[(na)CH3]/dt , a15=[(na)CH3]
% q = d[(na)2CH3]/dt , a16=[(na)2CH3]

% Constants and initial values
dydt = zeros(16,1);
pch4=0.45*101325;
par=pch4;
pna=0.1*101325;
N=6.023E17;
v0=0.001; %m^3
y1=pch4*v0/(8.314*973); %mol P1V0/RT
y5=pna*v0/(8.314*973); %mol P6V0/RT
yar=par*v0/(8.314*973); %mol ParV0/RT
n0=y1+y5+yar;
nt=a(1)+a(2)+a(3)+a(4)+a(5)+a(6)+a(7)+a(8)+a(9)+a(10)+a(11)+a(12)+a(13)+a(14)+a(15)+a(16)+yar;
v=nt*v0/n0;

%Rate constants
kf=[1.1*10^-8 2.2*10^9 1.5*10^11 N*2.2*10^-21 N*3.2*10^-22 N*10^-23 N*9.4*10^-27 N*2.3*10^-23 N*3.1*10^-21 N*6.9*10^-13 N*9.9*10^-24 9.5*10^9 8.7*10^10 6.8*10^3 8.4*10^5 7.1*10^3 5*10^9 5.8*10^8 4.4*10^-1 1.4*10^8 2.4*10^1 4.7*10^6 N*1.9*10^-10 N*2.8*10^-10];

kb=[N*4.1*10^-10 N*1.5*10^-9  N*4.4*10^-9 N*1.9*10^-10 4.1*10^9 N*1.9*10^-12 N*6.6*10^-11 3.9*10^2 4.7*10^4 N*5.7*10^-15 N*3.8*10^-19 N*3.3*10^-9 N*7.3*10^-9 N*3*10^-9  N*2.3*10^-10 N*1.4*10^-9 N*1.5*10^-9 N*3.1*10^-9 N*3.7*10^-12 N*5.1*10^-10 N*5*10^-10 N*3.9*10^-11 1.8*10^-10 2.7*10^-4];

%Modified rate constants
kfn=1.01*kf; %kf is changed(increased) by 1%
keq=kf./kb;
kbn=kfn./keq;

%timespan
tspan = [0 7200];
x_matrix = [];
% Calculate reaction rates

r1=-kf(1)*(a(1)/v)+kb(1)*(a(2)/v)*(a(3)/v);

r2=-kf(2)*(a(4)/v)+kb(2)*(a(5)/v)*(a(5)/v);

r3=-kf(3)*(a(6)/v)+kb(3)*(a(4)/v)*(a(5)/v);

r4=-kf(4)*(a(1)/v)*(a(5)/v)+kb(4)*(a(7)/v)*(a(2)/v);

r5=-kf(5)*(a(5)/v)*(a(1)/v)+kb(5)*(a(8)/v);

r6=-kf(6)*(a(1)/v)*(a(4)/v)+kb(6)*(a(9)/v)*(a(2)/v);

r7=-kf(7)*(a(1)/v)*(a(4)/v)+kb(7)*(a(10)/v)*(a(3)/v);

r8=-kf(8)*(a(1)/v)*(a(4)/v)+kb(8)*(a(11)/v);

r9=-kf(9)*(a(1)/v)*(a(6)/v)+kb(9)*(a(12)/v);

r10=-kf(10)*(a(1)/v)*(a(3)/v)+kb(10)*(a(2)/v)*(a(13)/v);

r11=-kf(11)*(a(1)/v)*(a(2)/v)+kb(11)*(a(14)/v)*(a(3)/v);

r12=-kf(12)*(a(8)/v)+kb(12)*(a(15)/v)*(a(3)/v);

r13=-kf(13)*(a(8)/v)+kb(13)*(a(7)/v)*(a(2)/v);

r14=-kf(14)*(a(7)/v)+kb(14)*(a(3)/v)*(a(5)/v);

r15=-kf(15)*(a(15)/v)+kb(15)*(a(5)/v)*(a(2)/v);

r16=-kf(16)*(a(9)/v)+kb(16)*(a(4)/v)*(a(3)/v);

r17=-kf(17)*(a(9)/v)+kb(17)*(a(5)/v)*(a(7)/v);

r18=-kf(18)*(a(10)/v)+kb(18)*(a(4)/v)*(a(2)/v);

r19=-kf(19)*(a(10)/v)+kb(19)*(a(5)/v)*(a(15)/v);

r20=-kf(20)*(a(16)/v)+kb(20)*(a(4)/v)*(a(7)/v);

r21=-kf(21)*(a(16)/v)+kb(21)*(a(6)/v)*(a(3)/v);

r22=-kf(22)*(a(16)/v)+kb(22)*(a(9)/v)*(a(5)/v);

r23=-kf(23)*(a(3)/v)*(a(3)/v)+kb(23)*(a(13)/v);

r24=-kf(24)*(a(2)/v)*(a(2)/v)+kb(24)*(a(14)/v);

b=v*(r1+r4+r5+r6+r7+r8+r9+r10+r11);

c=v*(-r1-r4-r6-r10+r11-r13-r15-r18+2*r24);

d=v*(-r1-r7+r10-r11-r12-r14-r16-r21+2*r23);

e=v*(r2-r3+r6+r7+r8-r16-r18-r20);   

f=v*(-2*r2-r3+r4+r5-r14-r15-r17-r19-r22);

g=v*(r3+r9-r21);

h=v*(-r4-r13+r14-r17-r20);

i=v*(-r5+r12+r13);

j=v*(-r6+r16+r17-r22);   

k=v*(-r7+r18+r19);

l=v*(-r8);

m=v*(-r9);

n=v*(-r10-r23);

o=v*(-r11-r24);

p=v*(-r12+r15-r19);

q=v*(r20+r21+r22);

%x=degree of rate control
%x=ki/r*(dr/dki)
%ki over here will be the new ki that is the changed ki. 
%xi=kfn*(dr/dki)/r

%b=v*(-kf(1)*(a(1)/v)+kb(1)*(a(2)/v)*(a(3)/v)-kf(4)*(a(1)/v)*(a(6)/v)+kb(4)*(a(8)/v)*(a(2)/v)-kf(5)*(a(6)/v)*(a(1)/v)+kb(5)*(a(11)/v)-kf(6)*(a(1)/v)*(a(5)/v)+kb(6)*(a(10)/v)*(a(2)/v)-kf(7)*(a(1)/v)*(a(5)/v)+kb(7)*(a(16)/v)*(a(3)/v)-kf(8)*(a(1)/v)*(a(5)/v)+kb(8)*(a(12)/v)-kf(9)*(a(1)/v)*(a(7)/v)+kb(9)*(a(13)/v)-kf(10)*(a(1)/v)*(a(3)/v)+kb(10)*(a(2)/v)*(a(4)/v)-kf(11)*(a(1)/v)*(a(2)/v)+kb(11)*(a(14)/v)*(a(3)/v));

x1  = (-kfn(1)*a(1)/b) + (kbn(1)*a(2)*(a(3)/v)/b);

x2  = 0;
 
x3  = 0;

x4  = (-kfn(4)*a(1)*(a(6)/v)/b) + (kbn(4)*a(8)*(a(2)/v)/b);

x5  = (-kfn(5)*a(6)*(a(1)/v)/b) + (kbn(5)*a(11)/b);

x6  = (-kfn(6)*a(1)*(a(5)/v)/b)+(kbn(6)*(a(10))*(a(2)/v)/b);

x7  = (-kfn(7)*(a(1))*(a(5)/v)/b)+(kbn(7)*(a(16)/v)*(a(3))/b);

x8  = (-kfn(8)*(a(1))*(a(5)/v)/b)+(kbn(8)*(a(12))/b);

x9  = (-kfn(9)*(a(1))*(a(7)/v)/b)+(kbn(9)*(a(13))/b);

x10 = (-kfn(10)*(a(1)/v)*(a(3))/b)+(kbn(10)*(a(2))*(a(4)/v)/b);

x11 = (-kfn(11)*(a(1)/v)*(a(2))/b)+(kbn(11)*(a(14))*(a(3)/v)/b);

x12 = 0;

x13 = 0;

x14 = 0;

x15 = 0;

x16 = 0;

x17 = 0;

x18 = 0;

x19 = 0;

x20 = 0;

x21 = 0;

x22 = 0;

x23 = 0;

x24 = 0;

x = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24]

dydt=[b;c;d;e;f;g;h;i;j;k;l;m;n;o;p;q];
