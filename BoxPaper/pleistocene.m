% Tziperman-Gildor Pleistocene climate model
%  with sea-ice

% Parameters


 A0=20e6*1e6 %ocean km2 to m2
 A1=20e6*1e6 %land km2 to m2

 H=350          %W/m2
 Cp=4000        %J/kg C
 VOL=21.6e6*1e9 %km3 to m3
 rho=1000.      %kg/m3
 C=Cp*VOL*rho

 P0=0.03*1e6 %m3/s  Precip
 P1=20e6     %m3/s (/Pa deleted)

 Aq=2.53e11  % Pa
 Bq=5.42e3   % K
 qr=0.7      % no units
 eq=0.622    % no units
 Ps=1e5      % Pa

 LEW=4000e3 %4000 km to m
 lam=10   %m

 I0S=0.3*A0 %max ice cover

 T0=-10+273 %degrees C
 T1=-2 +273

 S0=0.08e6    %m3/s
 SM=0.02e6    %m3/s
 ST=0.00075e6 %m3/s/K

 e=0.6        %no units
 sig=5.67e-8  %W/m2 K^-4
 Tscale=(H/(e*sig))^.25-273;  %scale

 al=0.8  %albedos
 as=0.6
 ac=0.3

%time steps
 N=2000000;  %20000  %number of steps
 kyr=365*86400*1e3;  %s
 dt=0.001*kyr;       %time step 1kyr
 Nkyr=kyr/dt;        %number of steps/kyr
 Tint=N*dt;          %integration time

%  save values every 1kyrs 

Nsave=N/Nkyr;  %number of values saved
Tsave=zeros(1,Nsave);
Vsave=zeros(1,Nsave);


% Initial conditions

 t=0;
 T=273.;
 V=1e16;  %1km thick * LEW^2

%Plot set-up

subplot(211)
PLOT1=plot(0,0,'.',...
'EraseMode','none','MarkerSize',2)
axis([0 N*dt/kyr -1 1])
xlabel('time (kyrs)')
ylabel('Land Ice and T')
hold on
% P and S
subplot(212)
PLOT2=plot(0,0,'.',...
'EraseMode','none','MarkerSize',2)
axis([0 N*dt/kyr 0 5])
xlabel('time (kyrs)')
ylabel('P and S')
hold on



% Iterate model for T(t),V(t)

for kk=1:Nsave,
 Vsave(kk)=V;
 Tsave(kk)=T;

for k=1:Nkyr,

 t=t+dt; %s

 Tf=T0+(T1-T0)*t/(N*dt);

 AS=I0S;
 if(T>Tf) 
  AS=0.;
 end;

 AL=LEW*2*(V/(LEW*2*lam^0.5))^.66;

 q = qr*eq*Aq*exp(-Bq/T)/Ps;

 P = (P0+P1*q)*(1-AS/A0);

 S = S0+SM*sin(2*pi*t/(41*kyr))+ST*(T-273);

 Ice= P - S ;

 Heat=(A0/C)*(-e*sig*T^4 + H*(1-as*AS/A0)...
                           *(1-al*AL/A0)...
                           *(1-ac) );

 V = V + Ice*dt;

 T = T + Heat*dt;

end

% Plot next point
 TC=T-273;
 tkyr=t/kyr;
subplot(211)
 set(PLOT1,'XData',tkyr,'YData',AL/A1)
 drawnow
 set(PLOT1,'XData',tkyr,'YData',TC/Tscale)
 drawnow
subplot(212)
 set(PLOT2,'XData',tkyr,'YData',P/P0)
 drawnow
 set(PLOT2,'XData',tkyr,'YData',S/S0)
 drawnow


end
