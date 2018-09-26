clear all

foldername='2018092311:39:33';
cd simulations
cd(foldername)
load('main_parameters.mat')
load('Results.mat')
cd ..
cd ..

%---------------------------
%----- Calculations --------
%---------------------------

E_pulse=@(t) -(E_0)*cos(omega*(t-average_time)+chirp_parm*(t-average_time).^2)...
    .*(cos(omega*(t-average_time)/(2*M)).^2).*...
    ((heaviside((t-average_time)+M*T/2))-heaviside((t-average_time)-M*T/2));


figure

%DOMAIN IN FREQ. SPACE
%We pick a domain taylored to each specific simulation. 

eps=10^(-4);
harmOrder=70;

Fs=length(time_vector)/((time_vector(end)-time_vector(1))*atomic_sec);
freqHz=(0:1:length(time_vector)-1)*Fs/(length(time_vector)); %Fs is the sampling rate (1/s)
freqat=freqHz*atomic_sec;
omega_vect=2*pi*freqat/omega;


k=5;

for i=k:k
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
    acceleration(1,1:Nost)=0;
    acceleration(1,1:Nost-1)=transpose(current_array(2:end,i))-transpose(current_array(1:end-1,i));
    acceleration(1,Nost)=acceleration(1,Nost-1);
    acceleration(1,:)=acceleration(1,:)/(time_vector(2)-time_vector(1));
    
    tau=(time_vector(1) + time_vector(end))/2;
    FWHM=(130*10^(-15))/atomic_sec;
    sigma_wide=FWHM/2.35482;
    J_Gabor(1,1:Nost)=exp(-((time_vector-tau).^2)/(sigma_wide^2));

    HHG(i,:)=abs(fft(acceleration.*J_Gabor)).^2;%'.*fft(current_load{1})');
    HHGConverged(i)=trapz(omega_vect(1:index),log10(HHG(1:index)));
    
    %plot(omega_vect(1:index),log10(HHG(1,1:index)))
    %hold on
    
    %Time frequency analysis: choose one k at a time...
    if i==k
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Perform a time frequency analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        sigma=1/(3*omega);
        window=@(t,tau_var) exp(-(t-tau_var).^2/(sigma^2));
        %figure
        %plot(time_vector,feval(window,time_vector,0))
        step=sigma; %(time_vector(end)-time_vector(1))/50;
        time_length=round((time_vector(end)-time_vector(1))/step);
        temp_time_array=linspace(time_vector(1),time_vector(end),time_length);


        %DOMAIN IN FREQ. SPACE
        %atomic_sec=2.418*10^(-17);
        Fs=length(time_vector)/((temp_time_array(end)-temp_time_array(1))*atomic_sec);
        freqHz=(0:1:length(time_vector)-1)*Fs/(length(time_vector)); %Fs is the sampling rate (1/s)
        freqat=freqHz*atomic_sec;
        omega_vect=2*pi*freqat/omega;
        omega_vect=omega_vect(1:index); % just pick a fixed value...


        %The matrix must be reduced in size, so pick a step for the
        %omega-dimension
        omegaStep = 5;

        %preallocate
        time_freq_array(1:length(omega_vect(1:omegaStep:end)),1:time_length)=0; %matrix containing the frequency spectrum for each time. 



        for j=1:time_length
            tau_var=temp_time_array(j);
            temp_accel=acceleration.*window(time_vector,tau_var);

            %FFT
            HHG_1=fft(temp_accel);
            HHG_1=abs(HHG_1).^2;
            HHG_1 = HHG_1(1:length(omega_vect(1:omegaStep:end)));
            time_freq_array(:,j)=log10(HHG_1');

        end 
    end
    
end


figure
surf(temp_time_array,omega_vect(1:omegaStep:end),time_freq_array)
xlabel('$t [a.u.]$','Interpreter','latex','FontSize',18)
%xlim([time_vector(1) time_vector(end)])
ylabel('$\frac{\omega}{\omega_0}$','Interpreter','latex','FontSize',18)
%xlim([0 omega_vect(end)])
title('$\frac{U}{t_0}=\frac{1}{2} \frac{Q}{t_0}=4$','Interpreter','latex','FontSize',18)
shading interp


STOP


figure 
plot(omega_vect(1:index+200),log10(HHG(6,1:index+200)))
hold off


STOP


%d_infinity




d_infinity(1:numberOfSimulations)=0;
indexArray(1:numberOfSimulations)=0;
eps=5*10^(-1);
for i=1:numberOfSimulations
    start_index=find(abs(time_vector-intervalVector(i))<eps);
    indexArray(i)=start_index;
    
    %calculate an average:
    dummy_vec=Dh_corr_array(start_index:end,i);
    avg = (1/length(dummy_vec))*sum(dummy_vec);
    
    d_infinity(i)=avg;
end


%We also need to calculate the slopes of D(t) at the various t_c (quench
%times)
S(1:numberOfSimulations)=0;

%The actual value of D(t) at t_c
V(1:numberOfSimulations)=0;
for i=1:numberOfSimulations
   if indexArray(i)>1
       S(i)=(Dh_corr_array(indexArray(i)-1,i)-Dh_corr_array(indexArray(i)-2,i) )/dt;
       V(i)=Dh_corr_array(indexArray(i),i);
   end
end


%-----------------------------
%----- Double axis plots -----
%-----------------------------

figure
yyaxis left
plot(intervalVector/T,d_infinity);
%title('$D(t), -\eta(t)$','Interpreter','latex','FontSize',18)
xlabel('Time [a.u.]','FontSize',18)
%xticks(avg_dur + [-4*T -3*T -2*T -T 0 T 2*T 3*T 4*T])
%xticklabels({'4','3','2','1','0','1','2','3','4'})
xlabel('Time [o.c.]','Interpreter','latex','FontSize',18)
xlim([intervalVector(1) intervalVector(end)]/T)


yyaxis right
plot(intervalVector/T,E_pulse(intervalVector)/E_0);
xlim([intervalVector(1) intervalVector(end)]/T)
leg1 = legend('$d_{\infty}(t)$','$E(t) (norm.)$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',17);
hold off






for i=1:numberOfSimulations
    if i==1
        [~,omega_vect,HHG_vector,time_freq_array,temp_time] = ...
            HHGDataProcessing(time_vector,current_array(:,i));
        
        %Preallocation:
        HHG_Array(1:length(HHG_vector),1:numberOfSimulations)=0;
        SZ=size(time_freq_array);
        time_freq_mega_array(1:SZ(1),1:SZ(2),1:numberOfSimulations)=0;
        
        
    else
        [~,~,HHG_vector,time_freq_array,temp_time] = ...
            HHGDataProcessing(time_vector,current_array(:,i));
        

    end
    HHG_Array(:,i)=HHG_vector;
    time_freq_mega_array(:,:,i)=time_freq_array;
end


%----------------------
%--HHG YIELD ---------
%---------------------

yieldArray(1:numberOfSimulations)=0;
for i=1:numberOfSimulations
   totYield=sum(HHG_Array(:,i));
   yieldArray(i)=totYield;
end

%Present the results
yyaxis left
plot(intervalVector,d_infinity);
%title('$D(t), -\eta(t)$','Interpreter','latex','FontSize',18)
xlabel('Time [a.u.]','FontSize',18)
%xticks(avg_dur + [-4*T -3*T -2*T -T 0 T 2*T 3*T 4*T])
%xticklabels({'4','3','2','1','0','1','2','3','4'})
xlabel('Time [a.u.]','Interpreter','latex','FontSize',18)
xlim([intervalVector(1) intervalVector(end)]/T)


yyaxis right
plot(intervalVector,yieldArray);
xlim([intervalVector(1) intervalVector(end)])
leg1 = legend('$d_{\infty}(t)$','$Total HHG-yield$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',17);
hold off


stop
%Plot the HHG-spectrum in subplot
figure
for i=1:numberOfSimulations

    subplot(5,4,i)
    plot(omega_vect,HHG_Array(:,i),'r')
    title('HHG spectrum, 10 sites','Interpreter','latex')
    xlabel('$\frac{\omega}{\omega_0}$','Interpreter','latex')
    %strLeg=strcat('critTime= ',num2str(intervalVector(i)),', U=',num2str(quenchVal));
    %legend(strLeg);
    %ylabel('Harmonic Intensity (log10) [arb.u]')
    axis([0 70 -7 2])
end





figure
for i=1:numberOfSimulations

    subplot(5,4,i)
    plot(time_vector,Dh_corr_array(:,i),'r')
    hold on;
    line([intervalVector(i) intervalVector(i)],[ min(Dh_corr_array(:,i)),max(Dh_corr_array(:,i))]);
    title('HHG spectrum, 10 sites','Interpreter','latex')
    xlabel('$\frac{\omega}{\omega_0}$','Interpreter','latex')
    strLeg=strcat('critTime= ',num2str(intervalVector(i)),', U=',num2str(quenchVal));
    legend(strLeg);

end

%The most clear ones. Figure 1
figure
subplot(1,2,1)
plot(time_vector/T,Dh_corr_array(:,1),'r')
xlabel('Time [o.c.]','FontSize',18);
ylabel('d(t)','FontSize',18);
legend({'$ t_c= 5.2893 [T]$'},'Interpreter','latex','FontSize',16)

subplot(1,2,2)
plot(time_vector/T,Dh_corr_array(:,2),'b')
hold off
xlabel('Time [o.c.]','FontSize',2);
legend({'$ t_c= 5.4999 [T]$'},'Interpreter','latex','FontSize',16)



figure
for i=1:numberOfSimulations

    subplot(3,2,i)
    plot(time_vector,W_array(:,i),'r')
    hold on;
    line([intervalVector(i) intervalVector(i)],[ min(Dh_corr_array(:,i)),max(Dh_corr_array(:,i))]);
    title('HHG spectrum, 10 sites','Interpreter','latex')
    xlabel('$\frac{\omega}{\omega_0}$','Interpreter','latex')
    strLeg=strcat('critTime= ',num2str(intervalVector(i)),', U=',num2str(quenchVal));
    legend(strLeg);

end



%time frequency


for i=16:16
    figure
    surf(temp_time,omega_vect,time_freq_mega_array(:,:,i))
    shading interp
    c=colorbar;
    c.Label.String='Harmonic intensity';
    c.Label.Interpreter='latex';
    c.Label.FontSize = 22;
    az=0;
    el=90;
    view(az,el);
    title('Undoped material','Interpreter','latex');
    xlabel('$\frac{\omega}{\omega_0}$','Interpreter','latex','FontSize', 18);
    ylabel('$time$','Interpreter','latex','FontSize', 24);
    axis([time_vector(1) time_vector(end)  0 omega_vect(end)])
    hold off

end