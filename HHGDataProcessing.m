%Function that takes a variable number of inputs, the most basic of which
%are time and acceleration. If the user inputs a 2D array, it is
%the data is processed by applying a "default" sliding window.


%Inputs
%time_vector: ...
%current: single current vector
%2Darray_x: this will belong to a parameter sweep. It is assumed that the
%time vector corresponding to this sweep is given as the first input argument. 

%varargout = [omega_vect,HHG_spectrum,HHGsurface] as an example
function varargout = HHGDataProcessing(varargin)
    
    %inputs
    [time_vector,J] = varargin{1:nargin};

    nOutputs = nargin + 3;
    varargout = cell(1,nOutputs);
    varargout{1}=nOutputs;
    
    
    %Atomic units
    atomic_sec=2.418*10^(-17); %secs
    freq=(32.9*10^12)*atomic_sec;
    omega=2*pi*freq; 


    %DOMAIN IN FREQ. SPACE
    %atomic_sec=2.418*10^(-17);
    Fs=length(time_vector)/((time_vector(end)-time_vector(1))*atomic_sec);
    freqHz=(0:1:length(time_vector)-1)*Fs/(length(time_vector)); %Fs is the sampling rate (1/s)
    freqat=freqHz*atomic_sec;
    omega_vect=2*pi*freqat/omega;
    fixed_cutOff=12000;
    omega_vect=omega_vect(1:fixed_cutOff);
    
    varargout{2}=omega_vect;

    switch nargin 
        case 2 %Processing a single current vector
            if 2+2==4
                tau=(time_vector(1) + time_vector(end))/2;
                FWHM=(90*10^(-15))/atomic_sec;%(125*10^(-15))/atomic_sec;
                sigma_wide=FWHM/2.35482;
                J_Gabor(1:length(time_vector))=exp(-((time_vector-tau).^2)/(sigma_wide^2));
            else
                J_Gabor(1:length(time_vector))=1;
            end


            %time frequency analysis
            %ASIZE=size(2Darray);

            %%%%%%%%%%%%%%%%%%%%%%%%
            %2D HHG SURFACE PLOT
            %%%%%%%%%%%%%%%%%%%%%%%%

            %compare with and without V-interaction
            HHG_surface(1:length(time_vector),1:length(omega_vect),1)=0;
            acceleration(1,1:length(time_vector))=0;
            acceleration(1:length(time_vector)-1)=J(2:end,1)-J(1:end-1,1);
            acceleration(length(time_vector))=acceleration(end);
            acceleration(:)=acceleration/(time_vector(2)-time_vector(1));
            HHG=fft(acceleration.*J_Gabor);
            HHG=(abs(HHG).^2);
            HHG_vector=transpose(log10(HHG(1:length(omega_vect))));
            
            varargout{3}=HHG_vector;
            
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
            omega_vect=omega_vect(1:fixed_cutOff); % just pick a fixed value...
            
            %preallocate
            time_freq_array(1:length(omega_vect),1:time_length)=0; %matrix containing the frequency spectrum for each time. 


            for i=1:time_length
                tau_var=temp_time_array(i);
                temp_accel=acceleration.*window(time_vector,tau_var);

                %FFT
                HHG_1=fft(temp_accel);
                HHG_1=abs(HHG_1).^2;
                HHG_1 = HHG_1(1:length(omega_vect));
                time_freq_array(:,i)=log10(HHG_1');

            end
            
            varargout{4} = time_freq_array;
            varargout{5}=temp_time_array;
            
        case 3 %This corresponds to a parameter sweep. 
            ASIZE=size(TWODIMarray);

            %%%%%%%%%%%%%%%%%%%%%%%%
            %2D HHG SURFACE PLOT
            %%%%%%%%%%%%%%%%%%%%%%%%

            %Pre-allocation
            HHG_surface(1:length(time_vector),1:ASIZE(2))=0;
            acceleration(1,1:length(time_vector))=0;

            for i=1:ASIZE(2)
                acceleration(1:length(time_vector)-1)=J(2:end,i)-J(1:end-1,i);
                acceleration(length(time_vector))=acceleration(end);
                acceleration(:)=acceleration/(time_vector(2)-time_vector(1));
                HHG=fft(acceleration.*J_Gabor);
                HHG=(abs(HHG).^2);
                HHG_surface(:,i)=transpose(log10(HHG(1:length(omega_vect))));
            end
            
            varargout{3} = HHG_surface;

        otherwise
            msg='Too many input arguments in function HHGsurface.m';
            error(msg)
    end
end




