%{
figure
[X,Y,Z] = sphere(16);
x = [0.5*X(:); 0.75*X(:); X(:)];
y = [0.5*Y(:); 0.75*Y(:); Y(:)];
z = [0.5*Z(:); 0.75*Z(:); Z(:)];
scatter3(x,y,z)
%}
%% Correlation growth
figure
mxc1=sort(mxc066_386L(1,:)');
mxc2=sort(mxc066_386L(2,:)');
mxc3=sort(mxc066_386R(1,:)');
mxc4=sort(mxc066_386R(2,:)');
h1=plot(mxc1(end-80:end))
hold on
legend('Left v=6.0 m/s','Left v=0.2 m/s','Right v=6.0 m/s','Right v=0.2 m/s')
ylabel('correlation level')
xlabel('80 maximal correlation coefficients')

h2=plot(mxc2(end-80:end))
h3=plot(mxc3(end-80:end))
h4=plot(mxc4(end-80:end))
%% Vertices of cortex old
LeftH=1;
if LeftH==1 
    v=VVLcortex.Vertices;
else
    v=VVRcortex.Vertices;
end
x=v(:,1);
y=v(:,2);
z=v(:,3);
scatter3(x,y,z)
%% Vertices of cortex
hf=figure;
if LeftH==1 
    vv=v(mxv066_386L(1,:),:);
else
    vv=v(mxv066_386R(1,:),:);
end
x=vv(:,1);  
y=vv(:,2);
z=vv(:,3);

SizeContour=0.1;
scatter3(x,y,z,SizeContour,'k','.')
grid off
%title('Jumping epicenters of traveling waves for the mesoscopic model')
title('Jumping epicenters of traveling waves for the macroscopic model')
set (hf, 'Position', [1200 500 750 500])
%axis off
hold on
if LeftH==1
  view(0,0)
else
  view(180,0)  
end    
%% if wave speed ==6.0 m/s vel=1 else if wave speed==0.2 m/s vel=2

vel=1;
%vel=2;
tresh=0.7;
vv=[];
vt=vv;
vvv=vv;
dtt=vv;
j=1;
if LeftH==1  
for i=1:size(mxv066_386L,2)
    if mxc066_386L(vel,i) > tresh
      vv(j,:)=v(mxv066_386L(vel,i),:);       
      vt(j)=i;
      j=j+1
    end
end
else
for i=1:size(mxv066_386L,2)
    if mxc066_386R(vel,i) > tresh
      vv(j,:)=v(mxv066_386R(vel,i),:);
      vt(j)=i;
      j=j+1
    end
end    
end

j=1;
vvv(j,:)=vv(j,:);
for i=1:length(vt)-1
   dt=(vt(i+1)-vt(i))*2;
   if dt>50
      dtt(j)=dt; 
      j=j+1; 
      vvv(j,:)=vv(i,:); 
   end
end

x=vvv(:,1);
y=vvv(:,2);
z=vvv(:,3);
SizeMarker=60;
scatter3(x,y,z,SizeMarker,'r','*')
if LeftH==1
  view(0,0)
else
  view(180,0)  
end    
hold off
%% The epicentrs traectory
%figure
h=line(x,y,z); 
h.Color='b';
h.LineWidth=1;
%% The epicentrs traectory animation
%figure
i=1;
h=line(x(1:2),y(1:2),z(1:2));
for i=3:length(x)
  h=line(x(1:i),y(1:i),z(1:i)); 
  h.Color='r';
  h.LineWidth=1;
  pause(dtt(i-1)/1000)
end
%% The epicentr time
for i=1:length(vt)-1
   vvt(i)=(vt(i+1)-vt(i))*0.002; 
end
bar(vvt')
%% The epicentr distace 
figure(1)
dd=zeros(length(x),1);
%j=1;
for i=1:length(x)-1
  dd(i)=dist3(x(i),y(i),z(i),x(i+1),y(i+1),z(i+1));
  %{
  if dd>0.005 
    ddd(j)=dd;
    j=j+1;
  end
  %}
end
plot(sort(dd))
%{
p=bar(dd);
p.BarWidth=0.5;
p.FaceColor=[0 0 0];
%}
%%
v = VideoWriter('peaks.avi');
open(v);
Z = peaks;
surf(Z); 
axis tight manual 
set(gca,'nextplot','replacechildren'); 
for k = 1:20 
   surf(sin(2*pi*k/20)*Z,Z)
   frame = getframe(gcf);
   writeVideo(v,frame);
end

close(v);
%%
v = VideoWriter('epicentr.avi');
open(v);
line(x(1:2),y(1:2),z(1:2)); 
axis tight manual 
set(gca,'nextplot','replacechildren'); 
for i=3:length(x)
   i=i 
   h=line(x(1:i),y(1:i),z(1:i)); 
   h.Color='b';
   h.LineWidth=1;
   
   frame = getframe(gcf);
   writeVideo(v,frame);
   
   pause(dtt(i-2)*2/10000)
end
close(v);
%%
D=pdist(vvv);
%%
SD=squareform(D)
%%
%[idx,C,sumd,D]=kmeans(xyz,Kclastrs);
% C - returns the k cluster centroid locations
% sumd - returns the within-cluster sums of point-to-centroid distances in the k-by-1 vector sumd
% D - returns distances from each point to every centroid in the n-by-k matrix D.
%%
Ld=[];
j=1;
i=j;
while j==1
   i=i+1; 
   [idx,C,sumd,D]=kmeans(vvv,i);
   j=all(sumd);
   if j==1
       Ld(i-1)=sum(sumd)/i;
   end   
end
plot(Ld)
%%
ylabel('correlation level')
xlabel('wave velocity')
paL=anova1(mxc066_386L');

paR=anova1(mxc066_386R');