function [BASIS,H_kin_left,H_kin_right,H_U_eff,H_U_ones,Q_matrix,H_spin]=ionicMatrixGen(sz,t_hopping,U,V,N,Q)

subb_up(1:sz)=0;
BASIS(1:sz^2)=0;


for g=1:length(subb_up) %looping through each basis vector
    
    flag=1;
    while 1   
        bin=0;
        occup_cnt=0;
        for i=0:N-1
            %make sure you dont have more than half-filling
            if not(gt(occup_cnt,N/2-1))
                
                binar=randi([0 1]);
                adder=binar*2^(i);
                if binar==1
                    occup_cnt=occup_cnt+1;
                else
                    continue
                end
                
            elseif i==N-1 && lt(occup_cnt,N/2)%this is the case where you have to insert ...
                binar=1;
                adder=binar*2^i;
                occup_cnt_=occup_cnt+1;
            elseif gt(occup_cnt,N/2-1)
                continue
            end
            
            bin=bin+adder;
        end
        
        %determine if you have already made this basis state
        if gt(g,1) && ismember(bin, subb_up)
            continue
        elseif not(ismember(bin,subb_up)) && occup_cnt==N/2
            flag=0;
        end
        
        if flag==0
            break
        end
        
    end
    subb_up(g)=bin;
end


%Sort the basis
%[B_sub_up,J]=sort(subb_up);
B_sub_up=subb_up;            
        
%%%%To complete the basis, create the spin-down sector. We are making unique configurations
%here.
M=length(subb_up);
%subb_dwn(1:sz)=0;
for i=1:M
    for j=1:M
        k=M*(i-1)+j;
        bin_up=B_sub_up(i);
        bin_dwn=B_sub_up(j); %why not?
        %associate a binary number to each configuration
        BASIS(k)=(2^N)*bin_dwn+bin_up;
    end
end
[BASIS,I]=sort(BASIS);





H_kin_left=sparse(length(BASIS),length(BASIS));
H_kin_right=sparse(length(BASIS),length(BASIS));
H_U=sparse(length(BASIS),length(BASIS));
H_U_ones=sparse(length(BASIS),length(BASIS));
H_V=sparse(length(BASIS),length(BASIS));
H_spin = sparse(length(BASIS),length(BASIS));
Q_matrix = sparse(length(BASIS),length(BASIS));

tic
for bas=1:length(BASIS)
    
    I_up=mod(BASIS(bas),2^N);
    I_dwn=(BASIS(bas)-I_up)/(2^N);
    state_up=dec2bin(I_up);
    state_dwn=dec2bin(I_dwn);

    if lt(length(state_dwn),N)
        diff=N-length(state_dwn);
        cat='';
        for irr=1:diff
            cat=strcat('0',cat);
        end
        state_dwn=strcat(cat,state_dwn);
    end
    
    if lt(length(state_up),N)
        diff=N-length(state_up);
        cat='';
        for irr=1:diff
            cat=strcat('0',cat);
        end
        state_up=strcat(cat,state_up);
    end

%     bas
%     state_up
%     state_dwn
%     
    
    %%%%%%%%%%%%%%%%%%%
    % HUBBARD SECTOR
    %%%%%%%%%%%%%%%%%%
    %SWEEP THROUGH POSITIONS AND COUNT THE DOUBLE OCCUPANCY
    %NOTE THAT WE DONT GET OFF-DIAGONAL ELEMENTS HERE. (BECAUSE OF MULTIPLICATION BY NUMBER OPERATORS ONLY)
   
    %counting the doublons, appending to each U-matrix.
    for us=1:1%NR_U
        for pos=1:N
            
            if pos<N
                modpos_right=pos+1;
            elseif pos==N
                modpos_right=1;
            end
            
            if str2double(state_dwn(pos))==str2double(state_up(pos)) && ...
                    not(str2double(state_dwn(pos))==0)
                H_U(bas,bas)=H_U(bas,bas)+U(us);
                H_U_ones(bas,bas)=H_U_ones(bas,bas)+1;
         
            end
            
            
            %Nearest Neighboor interaction    


            if str2double(state_dwn(pos))==str2double(state_dwn(modpos_right)) && ...
                not(str2double(state_dwn(pos))==0)
                H_V(bas,bas)=H_V(bas,bas)+V;
            end
                                
            if str2double(state_dwn(pos))==str2double(state_up(modpos_right)) && ...
                not(str2double(state_dwn(pos))==0)
                H_V(bas,bas)=H_V(bas,bas)+V;
            end
            
            if str2double(state_up(pos))==str2double(state_dwn(modpos_right)) && ...
                not(str2double(state_up(pos))==0)
                H_V(bas,bas)=H_V(bas,bas)+V;
            end    
            if str2double(state_up(pos))==str2double(state_up(modpos_right)) && ...
                not(str2double(state_up(pos))==0)
                H_V(bas,bas)=H_V(bas,bas)+V;
            end
            
            
        end
    end
    
    
    for i=1:N

        if i==1
            im=1;
        elseif not(i==1)
            im=mod(i,N);
        end
        
        %%%%%%%%%%%%%%%%%
        %SPIN UP SECTOR
        %%%%%%%%%%%%%%%%%
        %PICK OUT THE STATE FROM THE BASIS ABOVE

        new_st_up=state_up;
        %Vectors instead of strings

        new_st_up_v(1:N)=zeros(1,N);
        state_up_v(1:N)=zeros(1,N);
        for dum=1:N
            new_st_up_v(dum)=str2double(new_st_up(dum));
            state_up_v(dum)=str2double(state_up(dum));
        end  
        
%       action of a_dagg(i)a(i+1)
        prefac=2;
        
        if i<N
            if str2double(state_up(i))==0 && str2double(state_up(im+1))==0
                new_st_up_v(:)=0;
            elseif str2double(state_up(i))==1 && str2double(state_up(im+1))==1
                new_st_up_v(:)=0;
            elseif str2double(state_up(i))==0 && str2double(state_up(im+1))==1
                new_st_up_v(i)=1;
                new_st_up_v(im+1)=0;
            elseif str2double(state_up(i))==1 && str2double(state_up(im+1))==0
                new_st_up_v(:)=0;
            end
        else
            if str2double(state_up(i))==0 && str2double(state_up(im+1))==0
                new_st_up_v(:)=0;
            elseif str2double(state_up(i))==1 && str2double(state_up(im+1))==1
                new_st_up_v(:)=0;
            elseif str2double(state_up(i))==0 && str2double(state_up(im+1))==1
                new_st_up_v(i)=1;
                new_st_up_v(im+1)=0;
            elseif str2double(state_up(i))==1 && str2double(state_up(im+1))==0
                new_st_up_v(:)=0;
            end
            
            prefac=1;
        end
    
        %Create the corresponding matrix element
        new_st_dec=0;
        for dum=0:N-1
            adder=(new_st_up_v(N-dum))*2^(dum);
            new_st_dec=new_st_dec+adder;
        end
        %the spin down sector has remained the same, so its binary is still
        %I_dwn=I/2^N...
        
        %search in the basis for the corresponding element and its index. 
        for ind=1:length(B_sub_up)
            if new_st_dec==B_sub_up(ind)
                I_up_lft=B_sub_up(ind);
            end
        end
        
        
        %THE TOTAL BINARY VALUE:
        %I_up=mod(BASIS(bas),2^N);
        

        %this means...
        if not(all(new_st_up_v(:)==0))
            index=find(BASIS==(2^N)*I_dwn+I_up_lft);
            H_kin_left(index,bas)=H_kin_left(index,bas)-t_hopping*(-1)^prefac;
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %action of a_dagg(i+1)a(i)
        %%%%%%%%%%%%%%%%%%%%%%%%%

        %IMPORTANT TO INITIALIZE
        new_st_up=state_up;
        new_st_up_v=state_up_v;
        
        prefac=2;
        if i<N
            if str2double(state_up(i))==0 && str2double(state_up(im+1))==0
                new_st_up_v(:)=0;
            elseif str2double(state_up(i))==1 && str2double(state_up(im+1))==1
                new_st_up_v(:)=0;
            elseif str2double(state_up(i))==0 && str2double(state_up(im+1))==1
                new_st_up_v(:)=0;
            elseif str2double(state_up(i))==1 && str2double(state_up(im+1))==0
                new_st_up_v(i)=0;
                new_st_up_v(im+1)=1;
            end
        else
            if str2double(state_up(i))==0 && str2double(state_up(im+1))==0
                new_st_up_v(:)=0;
            elseif str2double(state_up(i))==1 && str2double(state_up(im+1))==1
                new_st_up_v(:)=0;
            elseif str2double(state_up(i))==0 && str2double(state_up(im+1))==1
                new_st_up_v(:)=0;
            elseif str2double(state_up(i))==1 && str2double(state_up(im+1))==0
                new_st_up_v(i)=0;
                new_st_up_v(im+1)=1;
            end
            prefac=1;
        end
            
        
        
        
        %Create the corresponding matrix element
        new_st_dec=0;
        for dum=0:N-1
            adder=(new_st_up_v(N-dum))*2^(dum);
            new_st_dec=new_st_dec+adder;
        end
        %search in the basis for the corresponding element and its index. 
        
        for ind=1:length(B_sub_up)
            if new_st_dec==B_sub_up(ind)
                I_up_lft=B_sub_up(ind);
            end
        end
        
        %THE TOTAL BINARY VALUE:
        %pick out the spin down state the binary above
        %I_up=mod(BASIS(bas),2^N);
        %I_dwn=(BASIS(bas)-I_up)/(2^N);
        %I_dwn=BASIS(bas)/(2^N);
        
        %this means...
        if not(all(new_st_up_v(:)==0))
            index=find(BASIS==(2^N)*I_dwn+I_up_lft);
            H_kin_right(index,bas)=H_kin_right(index,bas)-t_hopping*(-1)^prefac;
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CREATING THE Q-MATRIX
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if (str2double(state_up(i))==1)
            Q_matrix(bas,bas) = Q_matrix(bas,bas) + ((-1)^(i) *(Q/2));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % END - CREATING THE Q-MATRIX
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%
        %%SPIN DOWN SECTOR
        %%%%%%%%%%%%%%%%%%
        

        new_st_dwn=state_dwn;

        %Vectors instead of strings
        new_st_dwn_v(1:N)=zeros(1,N);
        state_dwn_v(1:N)=zeros(1,N);
        for dum=1:N
            new_st_dwn_v(dum)=str2double(new_st_dwn(dum));
            state_dwn_v(dum)=str2double(state_dwn(dum));
        end  
        
        %action of a_dagg(i)a(i+1)
        
        prefac=2;
        if i<N
            if str2double(state_dwn(i))==0 && str2double(state_dwn(im+1))==0
                new_st_dwn_v(:)=0;
            elseif str2double(state_dwn(i))==1 && str2double(state_dwn(im+1))==1
                new_st_dwn_v(:)=0;
            elseif str2double(state_dwn(i))==0 && str2double(state_dwn(im+1))==1
                new_st_dwn_v(i)=1;
                new_st_dwn_v(im+1)=0;
            elseif str2double(state_dwn(i))==1 && str2double(state_dwn(im+1))==0
                new_st_dwn_v(:)=0;
            end
        else
            if str2double(state_dwn(i))==0 && str2double(state_dwn(im+1))==0
                new_st_dwn_v(:)=0;
            elseif str2double(state_dwn(i))==1 && str2double(state_dwn(im+1))==1
                new_st_dwn_v(:)=0;
            elseif str2double(state_dwn(i))==0 && str2double(state_dwn(im+1))==1
                new_st_dwn_v(i)=1;
                new_st_dwn_v(im+1)=0;
            elseif str2double(state_dwn(i))==1 && str2double(state_dwn(im+1))==0
                new_st_dwn_v(:)=0;
            end
            
            prefac=1;
        end
        
        %search in the basis for the corresponding element and its index. 
        new_st_dec=0;
        for dum=0:N-1
            adder=(new_st_dwn_v(N-dum))*2^(dum);
            new_st_dec=new_st_dec+adder;
        end
        
        for ind=1:length(B_sub_up)
            if new_st_dec==B_sub_up(ind)
                I_dwn_lft=B_sub_up(ind);
            end
        end
        
        %THE TOTAL BINARY VALUE:
        %I_up=mod(BASIS(bas),2^N); %note; this is unchanged
        
        %this means...
        if not(all(new_st_dwn_v(:)==0))
            index=find(BASIS==(2^N)*I_dwn_lft+I_up);
            H_kin_left(index,bas)=H_kin_left(index,bas)-t_hopping*(-1)^prefac;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %action of a_dagg(i+1)a(i)
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        %IMPORTANT TO INITIALIZE
        new_st_dwn_v=state_dwn_v;
        new_st_dwn=state_dwn;
        
        prefac=2;
        
        
        if i<N
            if str2double(state_dwn(i))==0 && str2double(state_dwn(im+1))==0
                new_st_dwn_v(:)=0;
            elseif str2double(state_dwn(i))==1 && str2double(state_dwn(im+1))==1
                new_st_dwn_v(:)=0;
            elseif str2double(state_dwn(i))==0 && str2double(state_dwn(im+1))==1
                new_st_dwn_v(:)=0;
            elseif str2double(state_dwn(i))==1 && str2double(state_dwn(im+1))==0
                new_st_dwn_v(i)=0;
                new_st_dwn_v(im+1)=1;
            end
        else
            if str2double(state_dwn(i))==0 && str2double(state_dwn(im+1))==0
                new_st_dwn_v(:)=0;
            elseif str2double(state_dwn(i))==1 && str2double(state_dwn(im+1))==1
                new_st_dwn_v(:)=0;
            elseif str2double(state_dwn(i))==0 && str2double(state_dwn(im+1))==1
                new_st_dwn_v(:)=0;
            elseif str2double(state_dwn(i))==1 && str2double(state_dwn(im+1))==0
                new_st_dwn_v(i)=0;
                new_st_dwn_v(im+1)=1;
            end
            
            prefac=1;
        end
        
        
        %Create the corresponding matrix element
        new_st_dec=0;
        for dum=0:N-1
            adder=(new_st_dwn_v(N-dum))*2^(dum);
            new_st_dec=new_st_dec+adder;
        end
        %search in the basis for the corresponding element and its index. 
        
        for ind=1:length(B_sub_up)
            if new_st_dec==B_sub_up(ind)
                I_dwn_lft=B_sub_up(ind);   %the mistake was here!!!!!!!
            else
                continue
            end
        end
        %this means...
        
        if not(all(new_st_dwn_v(:)==0))
            index=find(BASIS==(2^N)*I_dwn_lft+I_up);
            H_kin_right(index,bas)=H_kin_right(index,bas)-t_hopping*(-1)^prefac;
        end
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CREATING THE Q-MATRIX
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if (str2double(state_dwn(i))==1)
            Q_matrix(bas,bas) = Q_matrix(bas,bas) + ((-1)^(i) *(Q/2));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % END - CREATING THE Q-MATRIX
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        
        
        

        
        
        
        
        
        
         
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % SPIN CORRELATION SECTOR
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
  
        
        
        %S_i^(+)S_j^(-)
        %IMPORTANT TO INITIALIZE
        new_st_up_v=state_up_v;
        new_st_dwn_v=state_dwn_v;
        
        if str2double(state_up(im+1))==0
            new_st_up_v(:)=0;
            new_st_dwn_v(:)=0;
        elseif str2double(state_up(im+1))==1
            if str2double(state_dwn(im+1))==1
                new_st_up_v(:)=0;
                new_st_dwn_v(:)=0;
            elseif str2double(state_dwn(im+1))==0
                new_st_up_v(im+1)=0;
                new_st_dwn_v(im+1)=1;
            end
        end
        
        
        %second step
        if str2double(state_dwn(i))==0
            new_st_up_v(:)=0;
            new_st_dwn_v(:)=0;
        elseif str2double(state_dwn(i))==1
            if str2double(state_up(i))==1
                new_st_up_v(:)=0;
                new_st_dwn_v(:)=0;
            elseif str2double(state_up(i))==0
                new_st_up_v(i)=1;
                new_st_dwn_v(i)=0;
            end
        end
        
        
        %Create the corresponding matrix element
        new_st_dec_dwn=0;
        new_st_dec_up=0;
        for dum=0:N-1
            adder=(new_st_dwn_v(N-dum))*2^(dum);
            new_st_dec_dwn=new_st_dec_dwn+adder;
            
            adder=(new_st_up_v(N-dum))*2^(dum);
            new_st_dec_up=new_st_dec_up+adder;
            
        end
        %search in the basis for the corresponding element and its index. 
        
        for ind=1:length(B_sub_up)
            if new_st_dec_dwn==B_sub_up(ind)
                I_dwn_lft=B_sub_up(ind);   %the mistake was here!!!!!!!
            else
                continue
            end
            
            if new_st_dec_up==B_sub_up(ind)
                I_up_lft=B_sub_up(ind);   %the mistake was here!!!!!!!
            else
                continue
            end
        end
        %this means...
        
        if not(all(new_st_dwn_v(:)==0)) %this is sufficient
            index=find(BASIS==(2^N)*I_dwn_lft+I_up_lft);
            H_spin(index,bas)=H_spin(index,bas) - 1/2; %always minus: see derivation
        end
        
        
        
        
        
        
        
        %S_i^(-)S_j^(+)
        %IMPORTANT TO INITIALIZE
        new_st_up_v=state_up_v;
        new_st_dwn_v=state_dwn_v;
        
        if str2double(state_dwn(im+1))==0
            new_st_up_v(:)=0;
            new_st_dwn_v(:)=0;
        elseif str2double(state_dwn(im+1))==1
            if str2double(state_up(im+1))==1
                new_st_up_v(:)=0;
                new_st_dwn_v(:)=0;
            elseif str2double(state_up(im+1))==0
                new_st_up_v(im+1)=1;
                new_st_dwn_v(im+1)=0;
            end
        end
        
        
        %second step
        if str2double(state_up(i))==0
            new_st_up_v(:)=0;
            new_st_dwn_v(:)=0;
        elseif str2double(state_up(i))==1
            if str2double(state_dwn(i))==1
                new_st_up_v(:)=0;
                new_st_dwn_v(:)=0;
            elseif str2double(state_dwn(i))==0
                new_st_up_v(i)=0;
                new_st_dwn_v(i)=1;
            end
        end
        
        
        %Create the corresponding matrix element
        new_st_dec_dwn=0;
        new_st_dec_up=0;
        for dum=0:N-1
            adder=(new_st_dwn_v(N-dum))*2^(dum);
            new_st_dec_dwn=new_st_dec_dwn+adder;
            
            adder=(new_st_up_v(N-dum))*2^(dum);
            new_st_dec_up=new_st_dec_up+adder;
            
        end
        %search in the basis for the corresponding element and its index. 
        
        for ind=1:length(B_sub_up)
            if new_st_dec_dwn==B_sub_up(ind)
                I_dwn_lft=B_sub_up(ind);   %the mistake was here!!!!!!!
            else
                continue
            end
            
            if new_st_dec_up==B_sub_up(ind)
                I_up_lft=B_sub_up(ind);   %the mistake was here!!!!!!!
            else
                continue
            end
        end
        %this means...
        
        if not(all(new_st_dwn_v(:)==0)) %this is sufficient
            index=find(BASIS==(2^N)*I_dwn_lft+I_up_lft);
            H_spin(index,bas)=H_spin(index,bas) - 1/2;
        end
        
        
        %THE Z-COMPONENT OF THE SPIN (atomic units)
        
        adder_i=(1/2) *(str2double(state_up(i)) - str2double(state_dwn(i)));
        adder_j=(1/2) *(str2double(state_up(im+1)) - str2double(state_dwn(im+1)));
        
        H_spin(bas,bas)=H_spin(bas,bas) + adder_i*adder_j;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % END - SPIN CORRELATION SECTOR
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
end


H_U_eff = H_U + H_V;



end