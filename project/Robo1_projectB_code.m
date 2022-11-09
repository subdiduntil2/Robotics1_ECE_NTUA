%% *** Robot (kinematic) model parameters *** 
clear all; 
close all; 
l0= 10.0;  %% in cm 
l1= 15.0;
l2= 30.0;
l3= 30.0;
R=5;

%% *** sampling period *** 
%% *** for the robot motion, kinematic simulation: 
dt = 0.001; %dt = 0.001; i.e. 1 msec)   

%% *** Create (or load from file) reference signals *** 
%% *** DESIRED MOTION PROFILE - TASK SPACE *** 
T=2; 
t=0:dt:T;



%xd0,td0,yd1: initial/final end-point position --> desired task-space trajectory  
ydA = -30.0;
ydB = -30.0;
ydC = -30.0;
yd_top = -25.0;
zdA = 40.00; 
zdB = 50.00;  
zdC = 45.0000001;
zd_top= 45.00;
xdA = 15.0;

% Example of desired trajectory : linear segment (x0,y0)-->(x1,y1); Time duration: Tf; 
disp('Initialising Desired Task-Space Trajectory (Motion Profile) ...'); %% 
disp(' ');   
xd(1) = xdA;
zd(1) = zdA; 
yd(1) = ydA;
xd_dot(1)=0;
yd_dot(1)=0;
zd_dot(1)=0;
a_z=45;
lambda_z=15;




for i=2:335
   t1=i/1002;
   zd(i) = zd(i-1)+ 0.5*lambda_z*dt;    
   zd_dot(i) = zd_dot(i-1)+a_z*dt;
   yd(i)=sqrt(R^2-(zd(i)-zdC)^2)+ydC;
   yd_dot(i)=(-(zd(i)-zdC).*zd_dot(i))./(sqrt(R^2-(zd(i)-zdC)^2));
   xd(i)=xdA;
   xd_dot(i)=0;
end  

for i=336:669
   t2=i/1002;
   zd(i) =  zd(i-1)+lambda_z*dt;    
   zd_dot(i)= lambda_z;
   yd(i)=sqrt(R.^2-(zd(i)-zdC).^2)+ydC;
   yd_dot(i)=(-(zd(i)-zdC).*zd_dot(i))./(sqrt(R^2-(zd(i)-zdC)^2));
   xd(i)=xdA;
   xd_dot(i)=0;
end 

for i=670:1002
   t3=i/1002;
   zd(i) = zd(i-1) + lambda_z*dt - 0.5*lambda_z*dt; 
   if zd(i)>50.0
      zd(i)=50;
   end 
   zd_dot(i)=zd_dot(i-1)-a_z*dt;
   yd(i)=sqrt(R.^2-(zd(i)-zdC).^2)+ydC;
   yd_dot(i)=(-(zd(i)-zdC).*zd_dot(i))./(sqrt(R^2-(zd(i)-zdC)^2));
   if i==1000 
       yd_dot(i)=-1.546;
   end
   if i==1001
       yd_dot(i)=-0.773;
   end
   if i==1002
       yd_dot(i)=0;
   end
   xd(i)=xdA;
   xd_dot(i)=0;
end  

for i=1003:1335
   zd(i) = zd(i-1) - 0.5*lambda_z*(dt); 
   if(zd(i)>50)
       zd(i)=50;
   end
   zd_dot(i)=zd_dot(i-1)+a_z*dt;
   yd(i)=yd(i-1);
   yd_dot(i)=0;
   xd(i)=xdA;
   xd_dot(i)=0;
end  

for i=1336:1668
   zd(i) = zd(i-1) - lambda_z*dt;    
   zd_dot(i)=lambda_z;
   yd(i)=yd(i-1);
   yd_dot(i)=0;
   xd(i)=xdA;
   xd_dot(i)=0;
end  

for i=1669:2001
   zd(i) = zd(i-1) - lambda_z*dt + 0.5*lambda_z*(dt);    
   zd_dot(i)=zd_dot(i-1)-a_z*dt;
   yd(i)=yd(i-1);
   yd_dot(i)=0;
   xd(i)=xdA;
   xd_dot(i)=0;
end  



save;  %% --> save data to 'matlab.mat' file   


fig1 = figure;  
subplot(2,2,1); 
plot(t,zd); 
ylabel('zd (cm)'); 
xlabel('time t (sec)');  

subplot(2,2,2); 
plot(t,yd); 
ylabel('yd (cm)'); 
xlabel('time t (sec)');  

subplot(2,2,[3:4]); 
plot(t,xd); 
ylabel('xd (cm)'); 
xlabel('time t (sec)'); 

fig2 = figure;  
subplot(2,2,1); 
plot(t,zd_dot); 
ylabel('uyd (cm)'); 
xlabel('time t (sec)');  

subplot(2,2,2); 
axis([0 2 -10 10])
axis on 
hold on 
plot(t,yd_dot); 
ylabel('uyd (cm)'); 
xlabel('time t (sec)');  

subplot(2,2,[3:4]); 
plot(t,xd_dot); 
ylabel('uxd (cm)'); 
xlabel('time t (sec)');

rd1 = xd(:).^2 + yd(:).^2;
A1=yd(:)+sqrt((rd1.^2)-(l1.^2));
A2=xd(:)+l1;
q1=2*atan(A1./A2);

A4=xd(:).*sin(q1)-cos(q1).*yd(:); %A0
A5=zd(:)-l0; %A1
s3=(A4.^2 +A5.^2 -l2.^2-l3.^2)./(2*l2*l3);
c3=sqrt(1-s3.^2);
q3=atan2(s3,c3);
A6=l3.*cos(q3); %A2
A7=l3.*sin(q3); %A3
c2=(A4.*A6+A5.*(A7+l2))./(A6.^2 +(A7.^2+l2).^2);
s2=(A5.*A6-A4.*(A7+l2))./(A6.^2 +(A7.^2+l2).^2);
q2=atan2(s2,c2);

%xd1=cos(q1).*l1;
%yd1=l1.*sin(q1);
%zd1=l0.*ones(2001,1);
%xd2=cos(q1).*l1-l2.*sin(q1).*sin(q2);
%yd2=l1.*sin(q1)+cos(q1).*l2.*sin(q2);
%zd2=l0+cos(q2).*l2;
%xd3=cos(q1).*l1+sin(q1).*(cos(q2+q3).*l3-l2.*sin(q2));
%yd3=l1.*sin(q1)+cos(q1).*(l2.*sin(q2)-l3.*cos(q2+q3));
%zd3=l0+cos(q2).*l2+l3.*sin(q2+q3);

xd1=zeros(2001,1);
yd1=zeros(2001,1);
zd1=l0.*ones(2001,1);
xd2=cos(q1).*l1;
yd2=l1.*sin(q1);
zd2=l0.*ones(2001,1);
xd3=cos(q1).*l1-l2.*sin(q1).*sin(q2);
yd3=l1.*sin(q1)+l2.*cos(q1).*sin(q2);
zd3=l0+cos(q2).*l2;

fig3 = figure; 
tdiff=(1:2000);
subplot(2,2,1); 
plot(tdiff./1000,1000*diff(q1)); 
ylabel('q1 dot (cm)'); 
xlabel('time t (sec)');  

subplot(2,2,2); 
plot(tdiff./1000,1000*diff(q2)); 
ylabel('q2 dot (cm)'); 
xlabel('time t (sec)');  ;  

subplot(2,2,[3:4]); 
plot(tdiff./1000,1000*diff(q3)); 
ylabel('q3 dot (cm)'); 
xlabel('time t (sec)');  






fig4=figure;
subplot(2,2,1); 
plot(t,q1); 
ylabel('q1 (rad)'); 
xlabel('time t (sec)');  

subplot(2,2,2);
plot(t,q2); 
ylabel('q2 (rad)'); 
xlabel('time t (sec)');  

subplot(2,2,[3:4]);
plot(t,q3); 
ylabel('q3 (rad)'); 
xlabel('time t (sec)'); 

fig5 = figure; 
axis([-50 50 -20 70 -40 40]) %%set zy plot axes (caution: square axes, i.e. dz=dy) 
axis on 
hold on 
%xlabel('z (cm)'); 
%ylabel('y (cm)'); 
loops=5;
plot3(xd,zd,yd,'rs'); 
dtk=200; %% plot robot position every dtk samples, to animate its motion 
plot3([0],[0],[0],'o'); 
for n=1:loops; 
  for tk=1:dtk:2001;    %%% 	
   pause(0.5);	%% pause motion to view successive robot configurations    
   plot3([0,xd1(tk)],[0,zd1(tk)],[0,yd1(tk)]);					
   plot3([xd1(tk)],[zd1(tk)],[yd1(tk)],'o'); 
   plot3([xd1(tk),xd2(tk)],[zd1(tk),zd2(tk)],[yd1(tk),yd2(tk)]);	
   plot3([xd2(tk)],[zd2(tk)],[yd2(tk)],'o');
   plot3([xd2(tk),xd3(tk)],[zd2(tk),zd3(tk)],[yd2(tk),yd3(tk)]);	
   plot3([xd3(tk)],[zd3(tk)],[yd3(tk)],'o');
   plot3([xd3(tk),xd(tk)],[zd3(tk),zd(tk)],[yd3(tk),yd(tk)]);
   plot3([xd(tk)],[zd(tk)],[yd(tk)],'g+');  
  end
end

fig6 = figure; 
axis([-20 70 -40 40]) %%set zy plot axes (caution: square axes, i.e. dz=dy) 
axis on 
hold on 
xlabel('z (cm)'); 
ylabel('y (cm)'); 
plot(zd,yd,'rs'); 
dtk=200; %% plot robot position every dtk samples, to animate its motion 
plot([0],[0],'o'); 
for n=1:loops;  
  for tk=1:dtk:2001;    %%% 	
    pause(0.5);	%% pause motion to view successive robot configurations    
    plot([0,zd1(tk)],[0,yd1(tk)]);					
    plot([zd1(tk)],[yd1(tk)],'o');    
    plot([zd1(tk),zd2(tk)],[yd1(tk),yd2(tk)]);	
    plot([zd2(tk)],[yd2(tk)],'o'); 
    plot([zd2(tk),zd3(tk)],[yd2(tk),yd3(tk)]);	
    plot([zd3(tk)],[yd3(tk)],'o');
    plot([zd3(tk),zd(tk)],[yd3(tk),yd(tk)]);
    plot([zd(tk)],[yd(tk)],'g+');  
  end
end



