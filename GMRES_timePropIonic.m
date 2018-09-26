%Time propagation scheme for finding the generalized minimum residual of
%Ax=b. Note that because of the mere inclusion of the first order in the
%forward and backward Taylor expansion, you must have dt<1. Preferrably 0.1
%or less.. The upshot of this method is that it is implicit, and should be
%able to tackle rapid parameter quenches. 
function [W,current,Dh_corr,spin_corr] = GMRES_timePropIonic(varargin)

if nargin ==10
[startVector,time_vector,...
    H_kin_left,H_kin_right,H_U_ones,Q_matrix,H_spin,U,foldername,f_tilde]=varargin{1:10};

    %overwriting the old one.
    interactionQuench=@(t) U;
    
end

%load parameters
cd(foldername)
load('main_parameters.mat','e','a','N')
cd ..

%Observables and time dependent quantities:
Nost = length(time_vector);

t=time_vector(1);
dt=time_vector(2)-time_vector(1);

x_tilde=startVector(:,1)/norm(startVector(:,1));


n=size(H_kin_left);

%Observables and time dependent quantities:
W(1:Nost,1)=0;
Dh_corr(1:Nost,1)=0;
spin_corr(1:Nost,1)=0;
current(1:Nost,1)=0;

%Tolerance for the time propagation
tol = 1e-12;
for i=1:Nost
    
    
   %Create the time-dependent hamiltonian - something wrong here! Note that
   %the interaction quench is time-dependent. My own propagator
%    H_tempxdt= (dt/6)*(6*Q_matrix + ...
%        H_U_ones*interactionQuench(t) +H_kin_left*f_tilde(t) + H_kin_right*conj(f_tilde(t))...
%         +4*H_kin_left*f_tilde((2*t+dt)/2) + 4*H_kin_right*conj(f_tilde((2*t+dt)/2)) +...
%         4*H_U_ones*interactionQuench((2*t + dt)/2)+...
%         H_kin_left*f_tilde(t+dt) + H_kin_right*conj(f_tilde(t+dt)) + H_U_ones*interactionQuench(t+dt)) ;
%    
   %Crank nicolson: (see the book with Oka and Aoki) or this link:
   %https://scholarworks.umass.edu/cgi/viewcontent.cgi?article=1504&context=theses
   %for the modified midpoint rule..
   H_tempxdt=@(t) (Q_matrix + ...
       H_U_ones*interactionQuench(t) +H_kin_left*f_tilde(t) + H_kin_right*conj(f_tilde(t))) ;

    
   %On GMRES form
   A=(speye(n(1)) + (1i*dt/2)*H_tempxdt(t+dt/2));
   b = (speye(n(1)) - (1i*dt/2)*H_tempxdt(t+dt/2) )*x_tilde;
   
   
   %Call the built in method and normalize the vector
   [x_tilde,~] = gmres(A,b,10,tol);
   
   
    %Denominator to ensure that the expectation value do not blow up. 
    denom=x_tilde'*x_tilde;

   
   %Observables and time dependent quantities:
    W(i,1)=(startVector(:,1)'*x_tilde)*(x_tilde'*startVector(:,1))/denom;

    curr_op= -1i*e*a*H_kin_left*f_tilde(t) + 1i*e*a*H_kin_right*conj(f_tilde(t));
    current(i,1)=(x_tilde')*curr_op*x_tilde/denom;
    

    %doublon number
    Dh_corr(i,1)=(1/N)*x_tilde'*H_U_ones*x_tilde/denom;


    %spin correlation
    spin_corr(i,1)=(1/N)*x_tilde'*H_spin*x_tilde/denom;
    
    t=t+dt;
    
    %Progress
    if mod(round(t-time_vector(1)),round((time_vector(end)-time_vector(1))/20))==0
        display(num2str(dt),'..is the time step.')
        display(num2str(((t-time_vector(1))/(time_vector(end)-time_vector(1)))*100),'percent')
    end
end