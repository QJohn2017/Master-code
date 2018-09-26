function [V_stationary,E]=ground_state(H_kin_left,H_kin_right,H_U,Q_matrix,...
NR_EXC_ST,L_flux,enableEnergyShift,energy_shift)


    tolerance=1e-05;


    %flux_steps=1;
    E(1:NR_EXC_ST,1)=0;
    SIZE=size(H_U);
    V_stationary=sparse(SIZE(1),NR_EXC_ST);%(1:SIZE(1),1:NR_EXC_ST)=0;


    diff_small=false;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %SPECIFY THE HAMILTONIAN IN THIS TIME STEP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if enableEnergyShift
        H_temp=energy_shift*speye(SIZE(1))  - (H_kin_left+H_kin_right+Q_matrix+H_U)  ;%DIAG(:,:,end);% +H_pol*A_pulse(1,i+1)/omega;
    else
        H_temp=H_kin_left + H_kin_right + H_U + Q_matrix;
    end


    max_index=min(4,L_flux);
    conv_energies(1:max_index,1)=102; %just a random number

    while diff_small==false

        Q_flux(1:SIZE(1),1)=rand(SIZE(1),1);
        Q_flux(:,1)=Q_flux(:,1)/norm(Q_flux(:,1));
        beta_flux(1:L_flux,1:L_flux)=0;


       for k=1:L_flux-1  %we let the loop decide when L is large enough... this is the time consuming?
            Q_flux(:,k+1)=H_temp*Q_flux(:,k);  %these represent first order vectors
            for j=max(1,k-2):k
                %create some matrix elements
                beta_flux(j,k)=(Q_flux(:,j))'*Q_flux(:,k+1);
                Q_flux(:,k+1)=Q_flux(:,k+1)-Q_flux(:,j)*beta_flux(j,k);
            end
            beta_flux(k+1,k)=norm(Q_flux(:,k+1));
            Q_flux(:,k+1)=Q_flux(:,k+1)/beta_flux(k+1,k);
        end


        %%%%%%%%%%%%
        [V_flux,D_flux]=eig(full(beta_flux));   %note that this gives normalized vectors.
        DD_flux=real(diag(D_flux));    %be careful here

        if enableEnergyShift
            [E_flux,I_flux]=sort(DD_flux,'descend'); %be aware of this
        else
            [E_flux,I_flux]=sort(DD_flux); 
        end
        V_flux=V_flux(:,I_flux);


        %tolerable?

        error=norm(E_flux(1:max_index)-conv_energies);
        conv_energies= E_flux(1:max_index);

        if error<tolerance || L_flux==SIZE(1)-1
            diff_small=true;
        else
            L_flux=L_flux+5;
        end
    end


%save an eigenvector of the ground state
    for p=1:NR_EXC_ST
        V_stationary(:,p)=Q_flux(:,:)*V_flux(:,p);
        if enableEnergyShift
            E(p,1)=-E_flux(p) + energy_shift;
        else
            E(p,1)=E_flux(p);
        end
    end

end

