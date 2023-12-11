%% Reinforcement Learning 
% for learning the temperature policy of the self-assemblies of patchy particles
%List of related codes:
% 1. RL_patchy.m (main code)
% 2. CallSimulation.m (function used in (1)) 
% 3. CalculateState.m (function used in (1))
% 4. PatchyRL_step1.f90 (used in (2))
% 5. PathcyRL.f90 (used in (2))

clearvars
workingfolder=pwd;%the folder contain the basic files
cd(workingfolder)

%% Set parameters for RL
    epoch=2;   % number of epochs
    eps=linspace(1,0,epoch);  %e psilon-greedy    
    nstep=2;   % number of update steps in each epoch                
    
    Tmin=0.2;Tmax=1.3;  % temperature range       
    T_int=0.1; %temperature mesh size
    s_int=0.1; %sigma mesh size
    A=[-0.05 0 0.05]; %actions on temperature
    s_tar=0.91 ;  %sigma value of target
    alpha=0.7;  gamma=0.9; %learning rate and discount factor
    Rmode=2; %reward function = -|s1-s_tar|^Rmode
    
%% Set parameters for Brownian Dynamics simulations:
    npm=256; %number of particles 
    N_BD=20000; %number of BD steps at each update (fixed T) 
    frac=0.5; %volume fraction in confined space, equivalant to area fraction of frac*3/2 
    ioutstep=10000; % write output after ioutstep steps

    
%% Others
    irng=randi(1000); rng(irng); %control the seed for random number
    timeBD=N_BD*1e-4;
%% Save the input:
    inputfile = fopen('00input.dat','wt');
    fprintf(inputfile,'%-10s\t  %u\n','epoch',epoch);
    fprintf(inputfile,'%-10s\t','epsilon'); fprintf(inputfile,'%.3f\t',eps);
    fprintf(inputfile,'\n%-10s\t  %u\n','nstep',nstep);
    fprintf(inputfile,'%-10s\t  %.3f\n','Tmin', Tmin);
    fprintf(inputfile,'%-10s\t  %.3f\n','Tmax', Tmax);
    fprintf(inputfile,'%-10s\t  %.3f\n','s_interval', s_int);
    fprintf(inputfile,'%-10s\t  %.3f\n','T_interval', T_int);    
    fprintf(inputfile,'%-10s\t', 'A'); fprintf(inputfile,'%.3f\t',A);
    fprintf(inputfile,'\n%-10s\t  %.3f\n','s_target',s_tar);
    fprintf(inputfile,'%-10s\t  %.3f\n','Reward mode ',Rmode);
    fprintf(inputfile,'%-10s\t  %.3f\n','alpha',alpha);
    fprintf(inputfile,'%-10s\t  %.3f\n','gamma',gamma);    
    fprintf(inputfile,'%-10s\t  %u\n','rng',irng);
    fprintf(inputfile,'%-10s\t  %u\n','npm',npm);
    fprintf(inputfile,'%-10s\t  %u\n','v.frac',frac);
    fprintf(inputfile,'%-10s\t  %u\n','ioutstep',ioutstep);
    fprintf(inputfile,'%-10s\t  %u\n','BDstep',N_BD);
    fprintf(inputfile,'%-10s\t  %.2f\n','timeBD',timeBD);
    fclose(inputfile);

%% Training process:    
Ss=(s_int/2):s_int:1;            
ST=(Tmin+T_int/2):T_int:Tmax;
Q=zeros(length(Ss),length(ST),length(A));   
for iepoch=1:epoch     
    ieps=eps(iepoch); 
    Statedata=zeros(nstep+1,3); %[istep sigma T]
    Actiondata=zeros(nstep,1);% action index

    s1=0.1*rand; T1=round(Tmin+(Tmax-Tmin)*rand,2); 
    s=s1; T=T1;
    
    xx1=abs(Ss-s1); xx2=abs(ST-T1);
    iss1=find(xx1==min(xx1)) ;    iss2=find(xx2==min(xx2)) ;
    iss1=iss1(randi(length(iss1))); iss2=iss2(randi(length(iss2))); 
    istep=0; 
    Statedata(istep+1,1:3)=[istep s1 T1]; %[istep sigma T]

    for istep=1:nstep
        eps_rand=rand; 
        if eps_rand >= ieps % choose action at maxQ    
            ia=find(Q(iss1,iss2,:)==max(Q(iss1,iss2,:)));
            if length(ia) ~= 1 
                ia=ia(randi(length(ia)));
            end                  
            T1=T+A(ia); %stay there if out of range of T
            if T1 < Tmin-0.001;     T1=T;
            elseif T1 > Tmax+0.001; T1=T;
            end
            T1=round(T1,2);        
        else % choose random action
            ia=randi(length(A)); T1=T+A(ia); 
            if T1 < Tmin-0.001;     T1=T;
            elseif T1 > Tmax+0.001; T1=T;
            end
            T1=round(T1,2);
        end

        % call the simulation at fixed T1:
        [irunfoldername]=CallSimulation(T1,istep,workingfolder,iepoch,npm,frac,ioutstep,N_BD,timeBD);    
        % calculate sigma of the simulation above:
        qcplot=CalculateState(irunfoldername); 
        cd(workingfolder)
        s1=qcplot(end,2); %sigma ratio of last step of the simulation above
        
        % calculate the reward:
        if Rmode==1;       rwd_1=-abs(s1-s_tar);
        elseif Rmode==2;   rwd_1=-(s1-s_tar)^2;
        end

        %update Q table by s_t, a_t, R_(t+1), s_(t+1) 
        xx1=abs(Ss-s1); xx2=abs(ST-T1);
        iss1_1=find(xx1==min(xx1)) ;    iss2_1=find(xx2==min(xx2)) ; 
        iss1_1=iss1_1(randi(length(iss1_1)));  iss2_1=iss2_1(randi(length(iss2_1)));     
        Q(iss1,iss2,ia)=  Q(iss1,iss2,ia)...
                    + alpha*(rwd_1 + gamma*max(Q(iss1_1,iss2_1,:))-Q(iss1,iss2,ia)) ;
        Statedata(istep+1,1:3)=[istep s1 T1]; %[istep sigma T]
        Actiondata(istep)=ia;% action index

       iss1=iss1_1;  iss2=iss2_1; T=T1; 
    end
%% save data each epoch
reshapeQ=reshape(Q,[],1);
save(strcat('train_Q_epoch',num2str(iepoch,'%u'),'.dat'),'reshapeQ','-ascii')

%% plot data s and T each epoch
figure; 
subplot(2,1,1);
plot(Statedata(:,1),Statedata(:,3),'.-k','LineWidth',1)
ylabel 'T'; ylim([Tmin-0.1 Tmax+0.1]) ;ytickformat('%.2f')
set(gca,'FontSize',14)
title(strcat('Training data, epoch',num2str(iepoch),', eps=',num2str(ieps,'%.2f')))

subplot(2,1,2);
plot(Statedata(:,1),Statedata(:,2),'.-m','LineWidth',1)
ylabel '\sigma'; ylim([0 1]); ytickformat('%.2f')
xlabel 'Update step'; 
set(gca,'FontSize',14) ;
savefig(strcat('fig_epoch',num2str(iepoch,'%u'),'_sT.fig')); close

%% plot policy each epoch
col(1,:)=[0 0 1]; col(3,:)=[1 0 0];col(2,:)=[0.5 0.5 0.5];
figure; hold on
for iSp=1:length(Ss)
for iSv=1:length(ST) 
    ia=find(Q(iSp,iSv,:)==max(Q(iSp,iSv,:)));
    if length(ia)==1
        scatter(Ss(iSp),ST(iSv),550,col(ia,:),'s','filled','MarkerEdgeColor',[0 0 0])
    end
end
end
axis equal; box on;  grid on;
xlim([0 1]);ylim([Tmin-0.1 Tmax+0.1]); 
xlabel '\sigma'; ylabel 'T';
title(strcat('Policy after epoch',num2str(iepoch,'%u')))
hold off
savefig(strcat('fig_epoch',num2str(iepoch,'%u'),'_policy.fig')); close
end







