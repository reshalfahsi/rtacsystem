clear all;
close all;
 
M = 1.3608;
m = 0.096;
R = 0.0592;
I = 2.175*10^-4;
k = 186.3;
 
%pembagi di persamaan
h = (M+m)*(I+0.5*m*R^2);
 
%nilai di matriks A
a11 = 0;
a12 = 1;
a13 = 0;
a14 = 0;
a21 = (-k*(I+m*R^2))/h;
a22 = 0;
a23 = -(((m^2)*(R^3)*(0.555*((M+m)*I+(M*m*R^2))+m*(0.616*I+0.555*m*R^2)))/(h^2));
a24 = (1.57*I*m*R+1.11*(m^2)*(R^3))/h;
 
a31 = 0;
a32 = 0;
a33 = 0;
a34 = 1;
 
a41 = (0.707*m*R*k)/h;
a42 = 0;
a43 = 0;
a44 = -((0.785*(m^2)*(R^2))/h);
 
%nilai di matriks B
b11 = 0;
b12 = 0;
b21 = (I+m*R^2)/h;
b22 = (-0.707*m*R/h);
b31 = 0;
b32 = 0;
b41 = (-0.707*m*R/h);
b42 = (M+m)/h;
 


%matriks A,B,C, dan D
A = [a11 a12 a13 a14; a21 a22 a23 a24; a31 a32 a33 a34; a41 a42 a43 a44];
B = [0 0; b21 b22; 0 0; b41 b42];
C = eye(4);
D = zeros(4,2);
 
sys_ss = ss(A,B,C,D); % representasi state space
 
%Diskretisasi Sistem
fs = 250; %Hz
Ts = 1/fs; % periode sampling (sekon)
 
sys_ss_d = c2d(sys_ss, Ts)
[Ad, Bd, Cd, Dd] = ssdata(sys_ss_d)
Ob=obsv(Ad,Cd); 
[n,m]=size(Ob); 
unob = m-rank(Ob); %observable
if(unob==0)
   disp('Given System is Observable.');
else
   disp('Given System is Unobservable');
end 
 
Co = ctrb(Ad,Bd);
[n,m]=size(Co);
unco=n-rank(Co); %controllable
if(unco==0)
   disp('Given System is Controllable.');
else
   disp('Given System is Uncontrollable');
end
 
%Mengecek kestabilan (melihat apakah pole berada di dalam unit circle)
pole = eig(Ad)
%Membuat stabil dengan menggunakan LQR
Qd=Cd'*Cd;
R=99999999999999999999*eye(2);
p=care(Ad,Bd,Qd,R);
K=inv(R)*Bd'*p
Anew=Ad-Bd*K; Bnew=Bd; Cnew=Cd; Dnew=Dd;
polenew=eig(Anew)
