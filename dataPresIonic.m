clear all

cd 'Results'
cd 'FIG1'
cd 'FIG1A'
load('main_parameters.mat')
load('Results_dt=0.1.mat')
cd ..
cd ..
cd ..


%---------------------------
%----- Calculations -------
%---------------------------



E_pulse=@(t) -(E_0)*cos(omega*(t-average_time)+chirp_parm*(t-average_time).^2)...
    .*(cos(omega*(t-average_time)/(2*M)).^2).*...
    ((heaviside((t-average_time)+M*T/2))-heaviside((t-average_time)-M*T/2));



%DOMAIN IN FREQ. SPACE
%atomic_sec=2.418*10^(-17);
Fs=length(time_vector)/((time_vector(end)-time_vector(1))*atomic_sec);
freqHz=(0:1:length(time_vector)-1)*Fs/(length(time_vector)); %Fs is the sampling rate (1/s)
freqat=freqHz*atomic_sec;
omega_vect=2*pi*freqat/omega;
omega_vect=omega_vect(1:end/6);


tau=(time_vector(1) + time_vector(end))/2;
FWHM=(125*10^(-15))/atomic_sec;
sigma_wide=FWHM/2.35482;
J_Gabor(1:length(time_vector))=exp(-((time_vector-tau).^2)/(sigma_wide^2));


acceleration(1,1:length(time_vector))=0;
acceleration(1:length(time_vector)-1)=current(2:end,1)-current(1:end-1,1);
acceleration(length(time_vector))=acceleration(end);
acceleration(:)=acceleration/(time_vector(2)-time_vector(1));
HHG=fft(acceleration.*J_Gabor);


figure
plot(omega_vect,log10(HHG(1:length(omega_vect))),'r')
xlim([0 60])
hold on     
  

%current
current_fig1a=real(current);




cd 'Results'
cd 'FIG1'
cd 'FIG1B'
load('main_parameters.mat')
load('Results_dt=0.1.mat')
cd ..
cd ..
cd ..


%---------------------------
%----- Calculations -------
%---------------------------

E_pulse=@(t) -(E_0)*cos(omega*(t-average_time)+chirp_parm*(t-average_time).^2)...
    .*(cos(omega*(t-average_time)/(2*M)).^2).*...
    ((heaviside((t-average_time)+M*T/2))-heaviside((t-average_time)-M*T/2));



%DOMAIN IN FREQ. SPACE
%atomic_sec=2.418*10^(-17);
Fs=length(time_vector)/((time_vector(end)-time_vector(1))*atomic_sec);
freqHz=(0:1:length(time_vector)-1)*Fs/(length(time_vector)); %Fs is the sampling rate (1/s)
freqat=freqHz*atomic_sec;
omega_vect=2*pi*freqat/omega;
omega_vect=omega_vect(1:end/6);


tau=(time_vector(1) + time_vector(end))/2;
FWHM=(125*10^(-15))/atomic_sec;
sigma_wide=FWHM/2.35482;
J_Gabor(1:length(time_vector))=exp(-((time_vector-tau).^2)/(sigma_wide^2));


acceleration(1,1:length(time_vector))=0;
acceleration(1:length(time_vector)-1)=current(2:end,1)-current(1:end-1,1);
acceleration(length(time_vector))=acceleration(end);
acceleration(:)=acceleration/(time_vector(2)-time_vector(1));
HHG=fft(acceleration.*J_Gabor);


plot(omega_vect,log10(HHG(1:length(omega_vect))),'k')
xlim([0 60])
leg1 = legend('$Q=2 t_0, U=0$','$Q=0, U=2t_0$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',17);
xlabel('$\frac{\omega}{\omega_0}$','Interpreter','latex','fontSize',20)
ylabel('HHG yield [arb.u]')
hold off           



%current
current_fig1b=real(current);

figure
plot(time_vector,current_fig1a,'k')
hold on
plot(time_vector,current_fig1b,'r')
leg1 = legend('$Q=2 t_0, U=0$','$Q=0, U=2t_0$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',17);
xlabel('time [a.u.]','fontSize',18)
ylabel('J(t) [a.u.]','fontSize',18)
hold off 

%----------------------
%--FIGURE 2 ---------
%---------------------



cd 'Results'
cd 'FIG2'
cd 'FIG2A'
load('main_parameters.mat')
load('Results_dt=0.1.mat')
cd ..
cd ..
cd ..



%current
current_fig2a=real(current);

acceleration(1,1:length(time_vector))=0;
acceleration(1:length(time_vector)-1)=current(2:end,1)-current(1:end-1,1);
acceleration(length(time_vector))=acceleration(end);
acceleration(:)=acceleration/(time_vector(2)-time_vector(1));
HHG=fft(acceleration.*J_Gabor);


figure
plot(omega_vect,log10(HHG(1:length(omega_vect))),'r')
xlim([0 70])
hold on     



cd 'Results'
cd 'FIG2'
cd 'FIG2B'
load('main_parameters.mat')
load('Results_dt=0.1.mat')
cd ..
cd ..
cd ..



%current
current_fig2b=real(current);


acceleration(1,1:length(time_vector))=0;
acceleration(1:length(time_vector)-1)=current(2:end,1)-current(1:end-1,1);
acceleration(length(time_vector))=acceleration(end);
acceleration(:)=acceleration/(time_vector(2)-time_vector(1));
HHG=fft(acceleration.*J_Gabor);



plot(omega_vect,log10(HHG(1:length(omega_vect))),'b')
xlim([0 70])
hold on  


cd 'Results'
cd 'FIG2'
cd 'FIG2C'
load('main_parameters.mat')
load('Results_dt=0.1.mat')
cd ..
cd ..
cd ..

%current
current_fig2c=real(current);


acceleration(1,1:length(time_vector))=0;
acceleration(1:length(time_vector)-1)=current(2:end,1)-current(1:end-1,1);
acceleration(length(time_vector))=acceleration(end);
acceleration(:)=acceleration/(time_vector(2)-time_vector(1));
HHG=fft(acceleration.*J_Gabor);



plot(omega_vect,log10(HHG(1:length(omega_vect))),'k')
xlim([0 70])
leg1 = legend('$Q=t_0, U=4t_0$','$Q=t_0, U=2t_0$','$Q=t_0, U=0.5t_0$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',17);
xlabel('$\frac{\omega}{\omega_0}$','Interpreter','latex','fontSize',20)
ylabel('HHG yield [arb.u]')
hold off     



%current
figure
plot(time_vector,current_fig2a,'k')
hold on
plot(time_vector,current_fig2b,'r')
hold on
plot(time_vector,current_fig2c,'b')
leg1 = legend('$Q=t_0, U=4t_0$','$Q=t_0, U=2t_0$','$Q=t_0, U=0.5t_0$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',17);
xlabel('time [a.u.]','fontSize',18)
ylabel('J(t) [a.u.]','fontSize',18)
hold off 



STOP 


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