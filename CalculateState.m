%% Extract the final snapshot and calculate the ratio sigma

function qcplot=CalculateState(irunfoldername)
disp('CalculateState...')
cd(irunfoldername)
defprm=importdata('00defprm.dat');  
npm=defprm(4);   %number of particles      
vol=defprm(5);   %volume fraction 
boxside=sqrt(2./3.*pi*npm/vol); Lx=boxside; Ly=Lx; % size of simulation box
pos=importdata('00pos.dat');            %position
ori0=importdata('00orient0.dat');       %orientation    
ori1=importdata('00orient1.dat');       %orientation    
ori2=importdata('00orient2.dat');       %orientation    
r_lim=2.5; %1st neighbour limit

istepmax=size(pos,1)/npm;
pbcx=importdata('00PBC.dat');
istepplot=istepmax;

%% extract final snapshot   
    coords =pos(end-npm+1:end,1:2); %now pos=(x,y)
    save('01LastPos.dat','coords','-ascii')
    coords0 =ori0(end-npm+1:end,1:3); %now pos=(x,y)
    coords1 =ori1(end-npm+1:end,1:3); %now pos=(x,y)
    coords2 =ori2(end-npm+1:end,1:3); %now pos=(x,y)
    save('01LastOri0.dat','coords0','-ascii')
    save('01LastOri1.dat','coords1','-ascii')
    save('01LastOri2.dat','coords2','-ascii')

%% Determine the local structure for the last snapshot
qcplot=zeros(length(istepplot),7); %matrix of istep|sigma|H|Z|U|  D1|D2
% IND=zeros(npm,length(istepplot)); %store ind at each timestep
islot=0;
for istep=istepplot 
    nmin=(istep-1)*npm+1; nmax=(istep-1)*npm+npm;
    coords =pos(nmin:nmax,1:2); 
    sh=pbcx(istep);
    
    [~,~,mnb]=neighbour_list_sh(coords, r_lim, Lx, Ly,sh);
    ind=zeros(npm,1); %local indicator of particle
    for ipm=1:npm   %the common of ipm and ipmneighbour
        mi=mnb(ipm,:); mi(mi==0)=[]; mi2=mi;
        jj=0;
        for jpm=mi            
            mj=mnb(jpm,:); mj(mj==0)=[];
            mij=intersect(mi,mj); %intersection
            jj=jj+1; mi2(jj)=length(mij);
        end
        msort=sort(mi2,'descend');
        strname='';
        for k=1:length(msort)
            strname=strcat(strname,num2str(msort(k)));
        end
        ind(ipm)=str2double(strname);
    end   
    
 
    % sigma = 21111 ; H=22110; Z=222222; D1= 22211; D2=222211
    islot=islot+1;
    qcplot(islot,1)=istep; 
    qcplot(islot,2)=sum(ind==21111)/npm; %sigma
    qcplot(islot,3)=sum(ind==22110)/npm; %H
    qcplot(islot,4)=sum(ind==222222)/npm;%Z
    qcplot(islot,6)=sum(ind==22211)/npm; %Undefined 1
    qcplot(islot,7)=sum(ind==222211)/npm;%Undefined 2
    qcplot(islot,5)=1-sum(qcplot(islot,[2:4 6 7]))   ; %Undefined 3  

end

end


function [listnb,rnb,mnb]=neighbour_list_sh(coords, r_lim, Lx, Ly,sh)
npm=size(coords,1); 
rmat=zeros(npm); % assign distant matrix     
for jpm = 1:(npm-1)                    
for ipm = (jpm+1):npm                     
    dr = coords(ipm,:) - coords(jpm,:)      ; 
    
    dr(1)= dr(1) - round(dr(1)/Lx)*Lx;          %if cross x  
    dr(1)= dr(1) - round(dr(2)/Ly)*Ly*sh;    %if cross y 
    dr(2)= dr(2) - round(dr(2)/Ly)*Ly ;      %if cross y      
    r=sqrt(dot(dr,dr));
    rmat(ipm,jpm)=r ;
    rmat(jpm,ipm)=r;
end
end

rnb=zeros(npm);% neighbour list matrix
mnb=zeros(npm);% neighbour list matrix using ipm, jpm 
for ipm=1:npm-1
for jpm=ipm+1:npm
    if (rmat(ipm,jpm)<r_lim) && (rmat(ipm,jpm)>1.5)
        rnb(ipm,jpm)=1; %register if the pair is in range
        rnb(jpm,ipm)=1;  
        mnb(ipm,jpm)=jpm;
        mnb(jpm,ipm)=ipm;
    end
end
end
listnb=sum(rnb,2); % number of neighbours for ipm
end