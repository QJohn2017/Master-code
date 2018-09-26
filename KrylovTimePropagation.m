function [W,current,Dh_corr,spin_corr]=KrylovTimePropagation(varargin)


    if nargin==15
        
        [V_stationary,...
            H_kin_left,H_kin_right,H_U_ones,H_spin,...
            L_time,Nost,t_0,t_1,dt,f_tilde,...
        e,a,N,interactionQuench] = varargin{1:nargin};
    
    else 
        msg='Invalid number of input arguments';
        error(msg)
    end
    

    SIZE=size(H_U_ones);
    W(1:Nost,1)=0;
    running_wf=zeros(SIZE(1),1);
    Q_ref=zeros(SIZE(1),L_time);%(1:SIZE(1),1:L_time)=0;
    Dh_corr(1:Nost,1)=0;
    spin_corr(1:Nost,1)=0;
    current(1:Nost,1)=0;
    

    tn=t_0;
    i=1;
    
    tic
    while tn<t_1

        Q_time=Q_ref;

        if i==1
            norm_wf=norm(V_stationary(:,1));
            Q_time(:,1)=V_stationary(:,1)/norm_wf;

        else
            norm_wf=norm(running_wf);
            Q_time(:,1)=running_wf/norm(running_wf);
        end



        if mod(round(tn-t_0),round((t_1-t_0)/20))==0
            display(num2str(((tn-t_0)/(t_1-t_0))*100),'percent')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %SPECIFY THE HAMILTONIAN IN THIS TIME STEP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %HIGHER ORDER CORRECTION TO MAGNUS PROPAGATOR
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


       % Q(1:SIZE(1),1:L+1)=0;
       % Q(:,1)=prop_wf(:,i)/norm(prop_wf(:,i));   %PSI/norm(PSI);
        beta=zeros(L_time,L_time);%(1:L_time,1:L_time)=0; 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %SPECIFY THE HAMILTONIAN IN THIS TIME STEP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


       H_temp= (dt/6)*(6*H_U_ones*interactionQuench(tn) +H_kin_left*f_tilde(tn) + H_kin_right*conj(f_tilde(tn))...
        +4*H_kin_left*f_tilde((2*tn+dt)/2) + 4*H_kin_right*conj(f_tilde((2*tn+dt)/2)) +...
        H_kin_left*f_tilde(tn+dt) + H_kin_right*conj(f_tilde(tn+dt))) ;


        for k=1:L_time-1  %we let the loop decide when L is large enough... this is the time consuming?
            Q_time(:,k+1)=H_temp*Q_time(:,k);  %these represent first order vectors
            for j=max(1,k-10):k
                %create some matrix elements
                beta(j,k)=(Q_time(:,j))'*Q_time(:,k+1);
                Q_time(:,k+1)=Q_time(:,k+1)-Q_time(:,j)*beta(j,k);
            end
            beta(k+1,k)=norm(Q_time(:,k+1));
            Q_time(:,k+1)=Q_time(:,k+1)/beta(k+1,k);

        end



        %%%%%%%%%%%%
        [S,D]=eig(beta);   %note that this gives normalized vectors.
        S_min=inv(S); %conjugate transpose: see (4) in https://arxiv.org/pdf/1301.7596.pdf

        %S=sparse(S);
        %D=sparse(D);
        %S_min=sparse(S_min);
        
        %%%%%%%%DO not confuse beta with the hessenberg matrix. the latter is found
        %%%%%%%%like this
        %Hess=Q*beta*Q';

        %%%%%%%%%%%final expression
        %running_wf=0*Q_time(:,1);
        M=zeros(L_time,L_time);%M(1:L_time,1:L_time)=0;
        for j=1:L_time
            %prop_wf=prop_wf/norm(prop_wf); %normalize in each step
            for k=1:L_time
                adder=S(j,k)*exp(-1i*D(k,k))*S_min(k,1);
                M(j,1)=M(j,1)+adder; %note: transpose! (matlab documentation)
            end
        end

        
        %bottleneck. Make the matrix sparse before doing the product. 
        Q_time=sparse(Q_time);
        M=sparse(M);
        running_wf(:,1)=Q_time*M(:,1);%the same as running_wf(:,1)=Q_time*M*Q_time'*Q_time(:,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CURRENT EXPECTATION VALUE 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

        curr_op=-1i*e*a*H_kin_left*f_tilde(tn) + 1i*e*a*H_kin_right*conj(f_tilde(tn));
        %curr_op=(1/2)*curr_op'+(1/2)*curr_op; %hermitian
        denom=running_wf'*running_wf;
        current(i,1)=(running_wf')*curr_op*running_wf/denom;


            W(i,1)=(V_stationary(:,1)'*running_wf)*(running_wf'*V_stationary(:,1));
       

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   CORRELATION SECTION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %doublon correlation
        Dh_corr(i,1)=(1/N)*running_wf'*H_U_ones*running_wf/denom;
        
        
        %spin correlation sector
        spin_corr(i,1)=(1/N)*running_wf'*H_spin*running_wf/denom;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  END - CORRELATION SECTION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%
%         % SECURE SAVING 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%
%         if saving_enabler==1
%             if mod(i,round(Nost/5))==0
%                 
%                 %print to file
%                 filename=strcat('State at time; t=',num2str(tn));
%                 StandardFileWriting(transpose(running_wf),filename,foldername)
% 
%                 %current
%                 filename=strcat('Current at time; t=',num2str(tn));
%                 StandardFileWriting(transpose(current(1:i,1)),filename,foldername)
% 
%                 %Doublon correlation
%                 filename=strcat('Doublon at time; t=',num2str(tn));
%                 StandardFileWriting(transpose(Dh_corr(1:i,1)),filename,foldername)
% 
%                 %Spin correlation
%                 filename=strcat('Spin at time; t=',num2str(tn));
%                 StandardFileWriting(transpose(spin_corr(1:i,1)),filename,foldername)
% 
%                 %Ground state
%                 filename=strcat('Ground state proj. at time; t=',num2str(tn));
%                  StandardFileWriting(transpose(W(1:i,1)),filename,foldername)
% 
%             end
%         end
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%
%         % END - SECURE SAVING 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tn=tn+dt;
        i=i+1;
        
    end
    
    
    toc
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE THE RESULTS
    %%%%%%%%%%%%%%%%%%%%%%%%
    %current
%     filename=strcat('Current');
%     StandardFileWriting(transpose(current(:,1)),filename,foldername)
% 
%     %Doublon correlation
%     filename=strcat('Doublon');
%     StandardFileWriting(transpose(Dh_corr(:,1)),filename,foldername)
% 
%     %Spin correlation
%     filename=strcat('Spin');
%     StandardFileWriting(transpose(spin_corr(:,1)),filename,foldername)
% 
%     %Ground state
%     filename=strcat('Ground state projection');
%     StandardFileWriting(transpose(W(:,1)),filename,foldername)
%     

    %%%%%%%%%%%%%%%%%%%%%%%%
    % END - SAVE THE RESULTS
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    %toc
end