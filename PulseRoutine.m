%EM pulse created by the following paramters (in this order)
%Avg_time: Average of the time interval
%M: Number of cycles
%T: Period of the carrier pulse
%E: Electric field magnitude

%For two pulses, specify the additional parameters
%M2: Number of cycles
%T2: Period
%E2: Electrif field magnitude
%Offs: Offset; time delay between each pulses maximum value

%Additional parameters
%phase1: CEP of the first pulse
%(note that Offs controls the second pulse)
function varargout = PulseRoutine(varargin)
    
    %Necessary calculations
    Avg_time=varargin{1};
    M = varargin{2};
    T = varargin{3};
    E_0 = varargin{4};
    chirp_parm = varargin{5};
    omega=(2*pi/T);
    
    
    switch nargin
        case 5
            A_pulse=@(t)  (E_0/omega)*sin(omega*(t-Avg_time) + chirp_parm*(t-Avg_time).^2).*(cos(omega*(t-Avg_time)/(2*M)).^2)...
            .*((heaviside((t-Avg_time) + M*T/2 ))-heaviside((t -Avg_time) - M*T/2 ));
            varargout{1} = A_pulse;
            %A_pulse=@(t)  (E_0/omega)*sin(omega*(t-Avg_time)).*(cos(omega*(t-Avg_time)/(2*M)).^2)...
            %.*(heaviside((t-Avg_time) + M*T/2 ))-heaviside((t -Avg_time) - M*T/2 );
            %varargout{1} = A_pulse;
        case 7
            A_pulse=@(t)  (E_0/omega)*sin(omega*(t-Avg_time) - pi/2).*(cos(omega*(t-Avg_time)/(2*M)).^2)...
            .*(heaviside((t-Avg_time) + M*T/2 ))-heaviside((t -Avg_time) - M*T/2 );
            varargout{1} = A_pulse;
            
            
            %second pulse
             A_pulse_sec=@(t)  (E_0/omega)*sin(omega*(t-Avg_time) + pi/2).*(cos(omega*(t-Avg_time)/(2*M)).^2)...
            .*(heaviside((t-Avg_time) + M*T/2 ))-heaviside((t -Avg_time) - M*T/2 );
            varargout{2} = A_pulse_sec;
            
        %case 8 to be implemented
            
        otherwise
            msg='Invalid number of input arguments';
            error(msg);
            
    end
end
            
        
