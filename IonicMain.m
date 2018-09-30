%follow conventions of the article cited above. 
%Number of sites enumerated 0,1,2,...,N-1
clear all

%atomic units:
atomic_sec=2.418*10^(-17); %secs
atomic_freq = 1/(2.418*10^(-17));
atomic_length=5.292*10^(-11); %m
atomic_energy=27.21; %eV
atomic_field=5.142*10^(11); %V/m
hbar=1;

a=4*10^(-10)/(5.29*10^(-11));
e=1;
t_hopping=0.52/atomic_energy;%hopping parm 


E_0=(10^9)/atomic_field; %omega*A_0
freq=(32.9*10^12)*atomic_sec;
omega=2*pi*freq;  
CEP=0;
M=10;
T=2*pi/omega;

%pulse and time-specs
dt=2;  %0.05
duration=T*M; %300*10^(-15)/atomic_sec; 300*10^(-15)/atomic_sec; %
t_0=0;

%The actual duration of the time propagation (not of the pulse)
t_1=duration;

time_vector=t_0:dt:t_1;
Nost=length(time_vector);



%OBS! If you have a hamiltonian whose ground state is located right near 0
%energy, then the ground state will not converge...

%----------------------
% USER PROMPTS 
%----------------------

N =input('Number of sites:  ');

enableKrylovGroundState = input('Calculate the ground state via the Krylov method? (true/false): ');

enableKrylovTimePropagation = input('Do time propagation with the Krylov method? (true/false): ');


%Number of configurations of one spin-species of electrons. 
sz=(factorial(N)/(factorial(N/2)*factorial(N/2)));


%Specifications for finding the ground state via the Krylov method.
%Krylov
L_time=15; %dimension of the Krylov space in time propagation
L_flux=15; %dimension of the Krylov space for obtaining the ground state



%----------------------
% DEFINING THE PULSE 
%----------------------

average_time=M*T/2; %(time_vector(end)+time_vector(1))/2;
A_pulse = PulseRoutine(average_time,M,T,E_0,0);
f_tilde=@(t) exp(-1i*e*a*(A_pulse(t)));





%----------------------
%   USER DEFINED INPUTS
%----------------------


enableParmSweep=input('enableParmSweep? (true/false): ');
if enableParmSweep
   U_array = input('U-array: ');
   Q_array = input('Corresponding Q-array: ');
   numberOfSimulations=length(U_array);

   %Note that we need to diagonalize for each time. 
else

    Q_factor = input('Value of Q/t_0? (Size of the gap) ');
    Q = Q_factor * t_hopping;

    U_factor = input('Value of U/t_0? ');
    U = U_factor * t_hopping;
end


enableConvergenceTest=input('enableConvergenceTest? (true/false): ');
if enableConvergenceTest
    convParmVector = input('Vector of dts: ');
end
% 
% 
% 
% 
%--------------------------------------------
%   GENERATING THE BASIS AND HAMTILTONIAN
%--------------------------------------------



[~,H_kin_left,H_kin_right,H_U_eff,H_U_ones,Q_matrix_ones,H_spin]=...
ionicMatrixGen(sz,t_hopping,1,0,N,1);



%--------------------------------------------
% SAVE WORKSPACE VARIABLES
%--------------------------------------------


foldername=num2str(yyyymmdd(datetime));
date=datestr(now,'HH:MM:SS');
foldername=strcat(foldername,date);
mkdir(foldername)

cd(foldername)
save('main_parameters.mat')
cd ..



interactionQuench =@(t) U;


if enableParmSweep
    
    %---------------------------------------
    % PARALLELIZATION FOR PARAMETER SWEEPS
    %---------------------------------------
    
    
    
    %Specify the required files to be made available to the parallell pool
    poolobj = parpool([2 8]); % "get current pool"
    addAttachedFiles(poolobj,{'PulseRoutine.m','heaviside.m','gmres.m',...
        'GMRES_timePropIonic.m','iterchk.m','iterapp.m'})
    
    
    %Preallocation
    W_array(1:length(time_vector),1:numberOfSimulations)=0;
    current_array(1:length(time_vector),1:numberOfSimulations)=0;
    Dh_corr_array(1:length(time_vector),1:numberOfSimulations)=0;
    spin_corr_array(1:length(time_vector),1:numberOfSimulations)=0;


    parfor i=1:numberOfSimulations
    


        %-----------------------------
        %       GROUND STATE
        %-----------------------------

        if enableKrylovGroundState
           [V_stationary,E]=ground_state(H_kin_left,H_kin_right,U_array(i) * H_U_ones,Q_array(i) * Q_matrix_ones,...
       1,L_flux,true,30);

     
        else
            [V_stationary, DD ] = eigs(...
                H_kin_left + H_kin_right+ U_array(i)*H_U_ones + Q_matrix_ones*Q_array(i),1,'sa');
            D = diag(DD);
            [E,I]=sort(D);


            V_stationary= V_stationary(:,I);
        end



        %Parfor stuff. Restating f_tilde...
        A_pulse = PulseRoutine(average_time,M,T,E_0,0);
        f_tilde=@(t) exp(-1i*e*a*(A_pulse(t)));


        %Regular time propagation
        [W,current,Dh_corr,spin_corr] = GMRES_timePropIonic(V_stationary,time_vector,...
    H_kin_left,H_kin_right,H_U_ones, Q_array(i) * Q_matrix_ones,H_spin,U_array(i),foldername,f_tilde);


        %Append to the preallocated arrays
        W_array(:,i)=W(:,1);
        current_array(:,i)=current(:,1);
        Dh_corr_array(:,i)=Dh_corr(:,1);
        spin_corr_array(:,i)=spin_corr(:,1);

    end

    cd(foldername)
    format compact;
    savedFileName='Results.mat';
    save(savedFileName,'W_array','current_array','Dh_corr_array','spin_corr_array')
    cd ..
    
elseif enableConvergenceTest 
    
    
    %-------------------------------------
    %   CONVERGENCE TEST IN THE TIME STEP
    %-------------------------------------
    

    if enableKrylovGroundState
       [V_stationary,E]=ground_state(H_kin_left,H_kin_right,U * H_U_ones,Q * Q_matrix_ones,...
   1,L_flux,true,30);


    else
        [V_stationary, DD ] = eigs(...
            H_kin_left + H_kin_right+ U*H_U_ones + Q_matrix_ones*Q,1,'sa');
        D = diag(DD);
        [E,I]=sort(D);


        V_stationary= V_stationary(:,I);
    end
    
    
    
    
    for i=1:length(convParmVector)
        clear time_vector
        dt = convParmVector(i);
        Nost = (t_1-t_0)/dt + 1;
        time_vector =linspace(t_0,t_1,Nost);
        


        if enableKrylovTimePropagation

           [W,current,Dh_corr,spin_corr] = KrylovTimePropagation(V_stationary,...
            H_kin_left,H_kin_right,H_U_ones,H_spin,...
            L_time,Nost,t_0,t_1,dt,f_tilde,...
        e,a,N,interactionQuench);
    
        else %Use GMRES

            [W,current,Dh_corr,spin_corr] = GMRES_timePropIonic(V_stationary,time_vector,...
        H_kin_left,H_kin_right,H_U_ones, Q * Q_matrix_ones,H_spin,U,foldername,f_tilde);

        end
        
        
        cd(foldername)
        format compact;
        savedFileName=strcat('Results','_dt=',num2str(dt),'.mat');
        save(savedFileName,'W','current','Dh_corr','spin_corr')
        cd ..         
    end
    
else
    
    [V_stationary, DD ] = eigs(...
    H_kin_left+H_kin_right + U*H_U_ones + Q_matrix_ones*Q,1,'sa');
    D = diag(DD);
    [E,I]=sort(D);

    V_stationary= V_stationary(:,I);
    
    interactionQuench =@(t) U;

    
    if enableKrylovTimePropagation

       [W,current,Dh_corr,spin_corr] = KrylovTimePropagation(V_stationary,...
        H_kin_left,H_kin_right,H_U_ones,H_spin,...
        L_time,Nost,t_0,t_1,dt,f_tilde,...
    e,a,N,interactionQuench);

    else %Use GMRES

        [W,current,Dh_corr,spin_corr] = GMRES_timePropIonic(V_stationary,time_vector,...
    H_kin_left,H_kin_right,H_U_ones, Q * Q_matrix_ones,H_spin,U,foldername,f_tilde);

    end
    
    
    cd(foldername)
    format compact;
    savedFileName=strcat('Results_dt=',num2str(dt),'.mat');
    save(savedFileName,'W','current','Dh_corr','spin_corr')
    cd ..

end



if not(enableParmSweep) && not(enableConvergenceTest) %Since this case involves the data processing of a multi-dimensional array. 
    %It can however be done by generalizing the code snippet below. 

    %---------------------------
    % HHG SPECTRUM
    %---------------------------

    %DOMAIN IN FREQ. SPACE
    %atomic_sec=2.418*10^(-17);
    Fs=length(time_vector)/((time_vector(end)-time_vector(1))*atomic_sec);
    freqHz=(0:1:length(time_vector)-1)*Fs/(length(time_vector)); %Fs is the sampling rate (1/s)
    freqat=freqHz*atomic_sec;
    omega_vect=2*pi*freqat/omega;
    omega_vect=omega_vect(1:end/6);

    %Window function
    tau=(time_vector(1) + time_vector(end))/2;
    FWHM=(125*10^(-15))/atomic_sec;
    sigma_wide=FWHM/2.35482;
    J_Gabor(1:length(time_vector))=exp(-((time_vector-tau).^2)/(sigma_wide^2));


    acceleration(1,1:length(time_vector))=0;
    acceleration(1:length(time_vector)-1)=current(2:end,1)-current(1:end-1,1);
    acceleration(length(time_vector))=acceleration(end);
    acceleration(:)=acceleration/(time_vector(2)-time_vector(1));
    HHG=fft(acceleration.*J_Gabor);


    %High harmonic spectrum
    plot(omega_vect,log10(HHG(1:length(omega_vect))))
    xlabel('$\frac{\omega}{\omega_0}$','Interpreter','latex','FontSize',20)
    ylabel('Harmonic intensity (arb.u)','Interpreter','latex','FontSize',18)
    
end
