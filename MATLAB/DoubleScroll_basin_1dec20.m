
%%
% PA & JN
% Matlab script for
% Double Scroll basin images
%
% Dec 2020
%
%
%% Set seed
tic;
clear variables;
rng(1811);

plotts=false;
printfigs=true;

%% params for functions
p.a=9;
p.b=14;

%% noise
p.sigma=0;

%% no shift
p.shift=0;
p.rate=0;
% p.rate varies

%% total duration
total=20;

%% timestep
h=1e-2;
ni=total/h;

%% set up figure
fignum=0;
%fignum=fignum+1;
%f1=figure(fignum);
%f1.PaperUnits='centimeters';
%f1.PaperSize=[12 12];
%f1.Units='centimeters';
%f1.InnerPosition=[1 1 12 12];

%% phase space scan
npts=21;
ny1=npts;
ny2=npts;
ny3=npts;
%
yis1=linspace(-5,5,ny1);
yis2=linspace(-2,2,ny2);
yis3=linspace(-5,5,ny3);


%%
basin=zeros(ny1,ny2,ny3);
Xg=zeros(ny1,ny2,ny3);
Yg=zeros(ny1,ny2,ny3);
Zg=zeros(ny1,ny2,ny3);

for i1=1:ny1
    %% report outer loop
    fprintf('i1=%d, ',i1);
    for i2=1:ny2
        for i3=1:ny3
            clf;
            % set initial condition
            yi=[yis1(i1),yis2(i2),yis3(i3)];
            ti=-total/2;
            [y,tt]=runmodel(yi,ti,total/2,ni,p);
            %[tt,y]=runmodel(ti,yi,total/2,h,p);
            yy=y(end,:)';
            if plotts==true
                %% plot the time series
                for ii=1:3
                    subplot(3,1,ii)
                    plot(tt,y(:,ii));
                    hold on;
                    xlabel('t');
                    ylabel(sprintf('x_%d',ii));
                    xlim([tt(1) tt(end)])
                    %ylim([0.45 0.55]);
                    subplot(3,1,1);
                    title(sprintf('(i1,i2,i3)=%d %d %d',i1,i2,i3));
                    drawnow();
                    hold off;
                end
            end
            endnorm=norm([3*yy(2),yy(3)]);
            % store grid and endnorm
            basin(i1,i2,i3)=endnorm;
            Xg(i1,i2,i3)=yi(1);
            Yg(i1,i2,i3)=yi(2);
            Zg(i1,i2,i3)=yi(3);
        end
    end
end

fignum=fignum+1;
%% figs of alpha, epsilon and Gamma
f2=figure(fignum);
clf;
f2.PaperUnits='centimeters';
f2.PaperSize=[12 12];
f2.Units='centimeters';
f2.InnerPosition=[1 1 12 12];

clf;
hold off;

%% Threshold for attractor
ythresh=8;

%% smooth and interpolate
b2=smooth3(basin);
it=1;
data = interp3(b2,it,'cubic');
Xgi=interp3(Xg,it,'linear');
Ygi=interp3(Yg,it,'linear');
Zgi=interp3(Zg,it,'linear');
fv = isosurface(Xgi,Ygi,Zgi,data,ythresh);
p2 = patch(fv,'FaceColor','red','FaceAlpha',0.5,'EdgeColor','none');
view(3);
%axis tight;
xlabel('x_1');
ylabel('x_2');
zlabel('x_3');
%axis tight;
camlight;
camlight(-80,-10);
lighting gouraud;
%title('Basin boundary');
hold on


%% initial value: chaos
[y,tt]=runmodel([0,0.6,0],0,200,5000,p);
ti=1000;
%plot3(y(1:ti,1),y(1:ti,2),y(1:ti,3),'y');
plot3(y(ti:end,1),y(ti:end,2),y(ti:end,3),'b');
%% initial value: periodic
[y,tt]=runmodel([0,1.0,0],0,60,2000,p);
plot3(y(1:ti,1),y(1:ti,2),y(1:ti,3),'c');
plot3(y(ti:end,1),y(ti:end,2),y(ti:end,3),'k');
toc

%keyboard

fname=sprintf('ds-basin-%s',date);
savefigure(fname,fignum,printfigs);

%% rotate around: thanks to https://uk.mathworks.com/matlabcentral/answers/86940-animate-3d-plot-view
az = 0;
el = 90;
view([az,el]);
degStep = 2.5;
detlaT = 0.1;
fCount = 141;
f = getframe(gcf);
[im,map] = rgb2ind(f.cdata,256,'nodither');
im(1,1,1,fCount) = 0;
k = 1;
% spin 45°
for i = 0:-degStep:-45
  az = i;
  view([az,el]);
  f = getframe(gcf);
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
  k = k + 1;
end
% tilt down
for i = 90:-degStep:15
  el = i;
  view([az,el]);
 f = getframe(gcf);
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
  k = k + 1;
end
% spin left
for i = az:-degStep:-120
  az = i;
  view([az,el]);
  f = getframe(gcf);
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
  k = k + 1;
end
% spin right
for i = az:degStep:0
  az = i;
  view([az,el]);
  f = getframe(gcf);
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
  k = k + 1;
end
% tilt up to original
for i = el:degStep:90
  el = i;
  view([az,el])
  f = getframe(gcf);
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
  k = k + 1;
end
imwrite(im,map,sprintf('%s_anim.gif',fname),'DelayTime',detlaT,'LoopCount',inf)

keyboard;


function [yo,tt]=runmodel(yi,ti,te,ni,p)
% runs from (yi,ti) to te in ni steps
ts=linspace(ti,te,ni);
[tt,yo]=ode45(@(t,y) FF(y,t,p),ts,yi);

end


function dy=FF(y,t,p)
%% RHS of ODE
% double scroll RHS
% with parameter shift in x

% shift y1
temp=p.sigma*Lambda(t,p.rate);
xx=y(1)-temp;
yy=y(2)-temp;
zz=y(3);

dy(1,1) = p.a*(yy-Phi(xx));
dy(2,1) = xx-yy+zz;
dy(3,1) = -p.b*yy;

end

function ll=Lambda(t,r)
temp=tanh(r*t);
ll=0.5*(1+temp);
end

function tt=Phi(xx)
% Hirsch Smale Devaney nonlinearity
%

tt=xx.^3/16-xx/6;

end

function savefigure(name,number,printfigs)
%
ffname = sprintf('./figs/%s-fig%i', name, number);

h=figure(number);
set(h,'Units','centimeters');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

if(printfigs==true)
    %savefig(ffname);
    print(sprintf('%s.pdf',ffname),'-dpdf','-r0');
end

end