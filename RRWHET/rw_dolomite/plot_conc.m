clear all
%geom=zeros(5,6);
%concfile=fopen('dolomite.con');
concfile='dolomite.con';
geom=dlmread(concfile);
times=(size(geom,1)-4)/3;
mom=zeros(times,4);
nx=geom(1,1);ny=geom(1,2);nz=geom(1,3);
dx=geom(1,4);dy=geom(1,5);dz=geom(1,6);
dirs=0;
for i=1:3
    if geom(1,i)>1.5
        dirs=dirs+1;
    end
end
nxy=nx*ny;
maxc=max(geom(7,:))
%domain=squeeze(zeros(nx,ny,nz));
domain=zeros(nx,ny,nz);

for l=1:times
    adv=3*l-3;
    t=geom(adv+5,1);
    nboxes=geom(adv+5,2);
% convert box numbers to i,j,k
    for n=1:nboxes
        ijk=geom(adv+6,n);
        kbox=floor((ijk-1)/nxy)+1;
        jbox=floor((ijk-1-(kbox-1)*nxy)/nx)+1;  
        ibox=ijk-(kbox-1)*nxy-(jbox-1)*nx;
        domain(ibox,jbox,kbox)=geom(adv+7,n);   % puts conc in right box
    end
% Plot in the right number of dimensions
  if dirs==1
    xpos=dx*geom(adv+6,1:nboxes);  xpos=xpos';
    cpos=domain(domain>0);
    figure(11);plot(xpos,cpos);
    axis([0 nx*dx 0 maxc])
    mean=sum(xpos.*cpos)/sum(cpos);
    var=sum(xpos.*xpos.*cpos)/sum(cpos)-mean.^2
    title(['Time = ',num2str(t),', Mean = ',num2str(mean),', Var = ',num2str(var)])
    domain=zeros(size(domain));
    outfile=['plots/plot',num2str(t)]
    print('-f11',outfile,'-dpdf');
    mom(l,1)=t; mom(l,2)=mean; mom(l,3)=var; mom(l,4)=2*1.2e-7*t;
    mom(l,5)=sum(cpos)*dx;
%    pause
  end
  if dirs==2
    contourf(squeeze(domain));  % not sure which dimension is singleton
  end
  if dirs==3  
  end
end   % time loop

figure; plot(mom(:,1),mom(:,2),'-o');
figure; plot(mom(:,1),mom(:,3),'-o'); hold on
plot(mom(:,1),mom(:,4));
figure(); plot(mom(:,1),mom(:,5));
    