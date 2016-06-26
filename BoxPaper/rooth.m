% Rooth model for S only, fixed T, V=1
%
%  ------   -------------    ------ 
% |      | |      2      |  |      |
% |  3   |  -------------   |   1  |
% |      |                  |      |
%  ------                    ------

% Parameters

 Pn=1.e-9
 Ps=1.e-9

 k=1e-6
 al=1.e-4
 be=1.e-3

% Initial conditions

 S1=35.
 S2=35.
 S3=35.
 T1=1.
 T3=0.

 q=k*(al*(T3-T1) - be*(S3-S1));
 dt=0.1/k;

%P=plot3(S1,S2,S3,'.',...
P=plot3(S1,q/(k*be),S3,'.',...
'EraseMode','none','MarkerSize',2)
%axis([34 36 34 36 34 36])
axis([34 36 -1 1 34 36])
xlabel('S1')
ylabel('q/(k*be)')
zlabel('S3')
hold on

for i=1:10000

 if q > 0
  ds = (  q*(S2-S1)-Pn )*dt;
 else 
  ds = ( -q*(S3-S1) - Pn)*dt;
 end

 S1=S1+ds

 if q > 0
  ds = (  q*(S3-S2)-(Pn+Ps) )*dt;
 else 
  ds = ( -q*(S1-S2)-(Pn+Ps) )*dt;
 end

 S2=S3+ds

 if q > 0
  ds = ( q*(S1-S3)-Ps )*dt;
 else 
  ds = ( -q*(S2-S3) - Ps)*dt;
 end

 S3=S3+ds;

% set(P,'XData',S1,'YData',S3,'XData',S3)
 set(P,'XData',S1,'YData',q/(k*be),'XData',S3)
 drawnow

 q=k*(al*(T3-T1) - be*(S3-S1));
 dt=0.1/k;

end
