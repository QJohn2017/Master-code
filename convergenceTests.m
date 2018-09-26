%Convergence testing for the GMRES propagator in dt. 


clear all

foldername='convergenceTests';
cd(foldername)
cd 'krylovConvergence'
load('main_parameters.mat')


%Allocate
%convParmVector=convParmVector(1:end-1);
%convParmVector = [convParmVector, 0.005, 0.001];

numOfFiles = 8;
for i=1:numOfFiles
    
   clear Dh_corr W current spin_corr
   loadFile=strcat('Results_dt=',sprintf('%.3f',convParmVector(i)),'.mat');
   load(loadFile)
   
   Dh_corr_load{i}=Dh_corr;
   W_load{i}=W;
   current_load{i}=current;
   spin_corr_load{i}=spin_corr;
   
   
end

cd ..
cd ..


%%%%%%%%%%
% BE CAREFUL HERE
%%%%%%%%%



%Calculate the convergence in some parameter.
Dh_converged(1:numOfFiles)=0;
spin_converged(1:numOfFiles)=0;
current_converged(1:numOfFiles)=0;
gs_converged(1:numOfFiles)=0;



for i=1:numOfFiles
    clear domain acceleration
    domain=linspace(t_0,t_1,length(Dh_corr_load{i}));

    
    Dh_converged(i) = trapz(domain,Dh_corr_load{i});
    spin_converged(i) = trapz(domain,spin_corr_load{i});
    
    %omega_vect(index)-omega_vect(1)
    %Calculate the acceleration
    acceleration(1:length(current_load{i}),1)=0;
    acceleration(1:length(current_load{i})-1)=current_load{i}(2:end,1)-current_load{i}(1:end-1,1);
    acceleration(length(current_load{i}))=acceleration(end);
    acceleration(:,1)=acceleration(:,1)/(time_vector(2)-time_vector(1));
    
    tau=(time_vector(1) + time_vector(end))/2;
    FWHM=(130*10^(-15))/atomic_sec;
    sigma_wide=FWHM/2.35482;
    J_Gabor(1:length(acceleration),1)=exp(-((linspace(t_0,t_1,length(acceleration))'-tau).^2)/(sigma_wide^2));
    
    current_converged(i) = trapz(domain,acceleration(:,1).*J_Gabor);
    gs_converged(i)=trapz(domain,W_load{i});  
    
    figure
    plot(linspace(t_0,t_1,length(acceleration)),real(current_load{i}))
    hold off
    

end



figure

subtr = 5;

subplot(2,2,1)
plot(convParmVector(end-subtr:end),Dh_converged(end-subtr:end)/max(Dh_converged(end-subtr:end)),'r*')
legend({'$ \int D(t) dt $'},'Interpreter','latex','FontSize',16)
xlabel('dt','FontSize',18)
title('Convergence test for N=6 and U/t=8','FontSize',16)

subplot(2,2,2)
plot(convParmVector(end-subtr:end),spin_converged(end-subtr:end)/max(spin_converged(end-subtr:end)),'b*')
legend({'$\int \eta(t) dt$'},'Interpreter','latex','FontSize',16)
xlabel('dt','FontSize',18)

subplot(2,2,3)
plot(convParmVector(end-subtr:end),current_converged(end-subtr:end)/max(current_converged(end-subtr:end)),'g*')
legend({'$\int j(t) dt$'},'Interpreter','latex','FontSize',16)
xlabel('dt','FontSize',18)

subplot(2,2,4)
plot(convParmVector(end-subtr:end),gs_converged(end-subtr:end)/max(gs_converged(end-subtr:end)),'k*')
legend({'$\int W(t) dt$'},'Interpreter','latex','FontSize',16)
xlabel('dt','FontSize',18)





HHGConverged(1:numOfFiles)=0;
%fixedCutOff=end;

eps=10^(-4);
harmOrder=100;
figure
for i=1:numOfFiles

    clear omega_vect index Fs time_vector freqat freqHz index HHG acceleration J_Gabor ...
        index tau
    
    %DOMAIN IN FREQ. SPACE
    %We pick a domain taylored to each specific simulation. 
    time_vector=linspace(t_0,t_1,length(current_load{i}));
    
    Fs=length(time_vector)/((time_vector(end)-time_vector(1))*atomic_sec);
    freqHz=(0:1:length(time_vector)-1)*Fs/(length(time_vector)); %Fs is the sampling rate (1/s)
    freqat=freqHz*atomic_sec;
    omega_vect=2*pi*freqat/omega;
    
    %omega_vect(end)-omega_vect(1)
    
    % We should go up to a certain energy... since the domain is not the
    % same..
    %omega_vect(end)-omega_vect(1) is not constant for the diff. simulations!
    
    %pick omega/omega_0 = 50;
    index = find(abs(omega_vect -harmOrder)<eps); 
    
    if length(index)>1
        index=index(1);
    elseif isempty(index)
        j=1;
        while isempty(index)
           j=j+1;
           index=find(abs(omega_vect -harmOrder)<j*eps);
        end
        
        if length(index)>1
            index=index(1);
        end
    end
    
    %omega_vect(index)-omega_vect(1)
    %Calculate the acceleration
    acceleration(1:length(current_load{i}),1)=0;
    acceleration(1:length(current_load{i})-1)=current_load{i}(2:end,1)-current_load{i}(1:end-1,1);
    acceleration(length(current_load{i}))=acceleration(end);
    acceleration(:,1)=acceleration(:,1)/(time_vector(2)-time_vector(1));
    
    tau=(time_vector(1) + time_vector(end))/2;
    FWHM=(130*10^(-15))/atomic_sec;
    sigma_wide=FWHM/2.35482;
    J_Gabor(1:length(acceleration),1)=exp(-((linspace(t_0,t_1,length(acceleration))'-tau).^2)/(sigma_wide^2));

    HHG=abs(fft(acceleration.*J_Gabor)).^2;%'.*fft(current_load{1})');
    HHGConverged(i)=trapz(omega_vect(1:index),log10(HHG(1:index)));
    
    plot(omega_vect(1:index),log10(HHG(1:index)))
    hold on
    
end
subtr = 5;

figure
plot(convParmVector(end-subtr:end),HHGConverged(end-subtr:end)/max(HHGConverged(end-subtr:end)),'-r*')
legend({'$ \int I(\omega) d\omega $'},'Interpreter','latex','FontSize',16)
xlabel('dt','FontSize',18)
title('Convergence test for N=6 and U/t=8','FontSize',16)
hold off