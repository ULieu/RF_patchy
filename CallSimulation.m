%% Call the simulation in Fortran90
% In this function, the folder "iepochXrunY" containing the data of 
% the Y^th step in X^th epoch created and the simulation is executed. 
% X and Y are iepoch and irun, respectively

function [irunfoldername]=CallSimulation(T,irun,folder,iepoch,npm,frac,ioutstep,istep_interval,timeMD)
cd(folder) 

%---No need to change these two lines:
    tinterval = 0.01; tmax = T; ibigemax =1; omg=2*pi/100; Pe0=0;
    lrank=5; mrank=5; ncode=4; rlist=7; 
    
%---Create folder "iepochXrunY" contain .f90 file
    system(strcat('mkdir iepoch',num2str(iepoch,'%03u'),'run', num2str(irun,'%03u'))) 
    fname=dir(strcat('iepoch',num2str(iepoch,'%03u'),'run',num2str(irun,'%03u'),'*')); %find the directory 
    %--for linux:
%         irunfoldername=strcat(fname.folder,'/',fname.name); 
%         cd(fname.name);    
%         if irun==1; copyfile ../PatchyRF_step1.f90             
%         else;       copyfile ../PatchyRF.f90     
%         end
    %--for windows:
        irunfoldername=strcat(fname.folder,'\',fname.name) ; 
        cd(fname.name);     
        if irun==1; copyfile ..\PatchyRF_step1.f90             
        else;       copyfile ..\PatchyRF.f90     
        end
    
    if irun~=1 %copy the last snapshot of previous folder to current one
        cd(folder)
        fnamep=dir(strcat('iepoch',num2str(iepoch,'%03u'),'run',num2str(irun-1,'%03u'),'*')); %find the dir of previous run 
        cd(fnamep.name); 
        posfile=dir('01LastPos.dat'); copyfile(posfile.name,irunfoldername)  %copy the latest snapshot to irun
        ori0file=dir('01LastOri0.dat'); copyfile(ori0file.name,irunfoldername)  %copy the latest snapshot to irun
        ori1file=dir('01LastOri1.dat'); copyfile(ori1file.name,irunfoldername)  %copy the latest snapshot to irun
        ori2file=dir('01LastOri2.dat'); copyfile(ori2file.name,irunfoldername)  %copy the latest snapshot to irun

        cd(irunfoldername)
        movefile '01LastPos.dat' '01LastPos_previous.dat'
        movefile '01LastOri0.dat' '01LastOri0_previous.dat'
        movefile '01LastOri1.dat' '01LastOri1_previous.dat'
        movefile '01LastOri2.dat' '01LastOri2_previous.dat'
    end
    
    cd(irunfoldername);    
    fprintf               ('%u\t  %u\t  %u\t  %u\t  %.5f\t  %u\t  %.2f\t %u\t %u\t %.5f\t %.5f\t %.5f\n', lrank, mrank, ncode, npm, frac, ioutstep,timeMD,ibigemax, istep_interval,tinterval, tmax, rlist)
    parainputfile = fopen('00defprm.dat','w');
    fprintf(parainputfile,'%4u\t  %4u\t  %u\t  %u\t  %.5f\t  %u\t  %.2f\t %u\t %u\t %.5f\t %.5f\t %.5f\n', lrank, mrank, ncode, npm, frac, ioutstep,timeMD,ibigemax, istep_interval,tinterval, tmax, rlist);
    fclose(parainputfile);    
    
    seedfile = fopen('00seedparfor.dat','w');  
    fprintf(seedfile,'%u\n', irun);             
    fclose(seedfile);                           
    
    ratefile = fopen('00Pe_osci.dat','w');         
    fprintf(ratefile,'%.5f\t %.5f\n', Pe0,omg);             
    fclose(ratefile);                           

%---Run simulation 
    %--in linux:   
%         if irun==1
%         system('/opt/intel/bin/ifort -O3 -qopt-matmul PatchyRL_step1.f90 -o PatchyRL_step1')
%         system('./PatchyRL_step1') 
%         system('rm ./PatchyRL_step1')    
%         else
%         system('/opt/intel/bin/ifort -O3 -qopt-matmul PatchyRF.f90 -o PatchyRF')
%         system('./PatchyRF') 
%         system('rm ./PatchyRF') 
%         end
%         system('rm ./*.mod') 
%         system('rm ./*.f90') 

    %--in windows:  
        if irun==1
        system('ifort /O3 /Qopt-matmul PatchyRF_step1.f90')
        system('PatchyRF_step1') 
        else
        system('ifort /O3 /Qopt-matmul PatchyRF.f90')
        system('PatchyRF') 
        end
       
cd(folder)    
end


