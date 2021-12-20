clc
%clear all

%=========================================================================
%-- This code simulate the Analytical expression for UL SE, in Theorem 1.--
%-- In order to simulate the analytical expressions for for the downlink SE
%-  in Theorem-2,  X_{pq,ij} and Kapa_{pq,ij} must be btained from the
%-  "OTFS_DD_channel_Xqij function", and accordingly the related variables, 
%-  i.e."var" and "Iq1" must be multipled by "Xpq(m,k)" and "Kpq(m,k)", 
%-  respectively.

% [Xpq,Kpq] = OTFS_DD_channel_Xqij(Lpq_vec);
%==========================================================================

global  M_a K_u M N deltaf elmax kmax Lpqmax
%% Inital parameters
uncor   = 0; % Set "uncor=0" to simulate correlated and "uncor=1" for uncrrelated

EVA     = 0;  %Set EVA =1 to model extended vehicular A model
EVB     = 1;  %Set EVB =1 to model extended vehicular B model
Derrick = 0;

%==========================================================================
%==========================================================================
%%%----------------------- OTFS parameters --------------------------------
%==========================================================================
%==========================================================================
% We consider an OTFS system with $\Delta f=15$ kHz, carrier frequency 
% $f_c=4$ GHz and vehicular speed of $500$ kmph. We consider two 3GPP 
% vehicular models, namely:
% (i) extended vehicular A (EVA) with $\ell_{pq}=9$ and $\tau_{max}=2.51~\mu sec$ 
% (ii) extended vehicular B (EVB) with  $\ell_{pq}=6$ and $\tau_{max}=20~\mu sec$.
% The block duration is assumed as $N=128$ which corresponds to $T_f=8.85~msec$,
% and the number of sub-carriers is $M=512$.

if EVA==1
    
    %channel parameters for vehicular speed of up to 500 kmph
    vspeed  = 300;            % vehicular speed in kmph
    clight  = 3*10^8;         % speed of light in m/s
    deltaf  = 15*10^3;        % Frequency spacing between adjacent sub-carriers
    f       = 4000;           % Carrier frequency in MHz (1900)
    N  = 128;                 % number of symbol N  = 128;
    M  = 512;                 % number of subcarriers
    MN = N*M;                 % number of symbols per frame
    
    taumax = 2.5*10^-6;  %maximum delay spreed according to EVA,
    elmax  = floor(taumax*deltaf*M);
    
    numax  = f *vspeed ./(3.6*clight/1000000);  %f is in kHz
    kmax   = floor (numax*N/deltaf);
    khat    = 1;   % reduced guard interval length
    
    Lpqmax  = 9;   % maximum number of path between the APs and uers
    
elseif EVB==1
    
    %channel parameters for vehicular speed of up to 500 kmph
    vspeed  = 300;           % vehicular speed in kmph
    clight  = 3*10^8;        % speed of light in m/s
    
    deltaf = 15*10^3;        % Frequency spacing between adjacent sub-carriers
    f  = 4000;               % Carrier frequency in MHz (1900)
    N  = 128;                % number of symbol
    M  = 512;                % number of subcarriers
    MN = N*M;                % number of symbols per frame
    
    taumax = 10*10^-6;  %maximum delay spreed according to EVA,
    elmax  = floor(taumax*deltaf*M);
    
    numax  = f *vspeed ./(3.6*clight./1000000);  %f is in kHz
    kmax   = floor (numax*N/deltaf);
    khat   = 1;   % reduced guard interval length
    
    Lpqmax = 6;   % maximum number of path between the APs and uers
else
    
%  TWC.2021.Performance Analysis of Coded OTFS Systems over High-Mobility Channels
%  Integer delay and Doppler case and set the maximum delay index
%  as lmax = 3 and the maximum Doppler index as kmax = 5, which is
%  corresponding to a relative speed around 250 km/h with 4 GHz
%  carrier frequency and 1.5 kHz sub-carrier spacing.
    
    %channel parameters for vehicular speed of up to 250 kmph
    vspeed  = 350;             % vehicular speed in kmph
    clight  = 3*10^8;          % speed of light in m/s
    deltaf  = 15*10^3;        % Frequency spacing between adjacent sub-carriers
    f       = 4000;            % Carrier frequency in MHz (1900)
    tumax   = 4.7*10^-6;
    
    N  = 30;                   % number of symbol
    M  = 40;                   % number of subcarriers
    MN = N*M;                  % number of symbols per frame
    
    numax  = f *vspeed ./(3.6*clight/1000000);  %f is in MHz
    kmax   = floor (numax*N/deltaf);
    
    elmax   = floor (tumax*M*deltaf);   % maximum delay index
    %kmax    = 2;   % maximum Doppler index
    khat    = 1;   % reduced guard interval length
    Lpqmax  = 5;   % maximum number of path between the APs and uers
end

% For each channel realization, we randomly select the delay and Doppler
% indices according to the uniform distribution, such that we have
% -kmax<ki<kmax and 0<li<lmax

Nguard = (4*kmax+4*khat)*(2*elmax+1);
MN./Nguard
%==========================================================================
%==========================================================================
%%%-------------------  Cell free Network parameters ----------------------
%==========================================================================
%==========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uplink
% Consider a square are of DxD m^2
% M_a distributed APs serves K_u terminals, randomly located in the area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_a   = 100;     %number of access points
K_uv  = [5:5:50];              %number of terminals

D =1;        %in kilometer

B  = deltaf*M;  
Hb = 15;      % Base station height in m
Hm = 1.65;    % Mobile height in m
%f  = 4000;    % Frequency in MHz (1900)
aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8);
L  = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL;


power_f     = 1; 

Pmax        = 1;
Paloc_sp    = 2.5*10^-3;
Paloc_emp   =0.25;
power_usp   = Paloc_sp*Pmax;    %uplink power: 200 mW
power_uemp  = Paloc_emp*Pmax;    %uplink power: 200 mW
power_p_sp  = (1-Paloc_sp)*Pmax; 
power_p_em  = 2*(1-Paloc_emp)*Pmax; 

noise_p    = 10^((-203.975+10*log10(B)+9)/10); %noise power

Pd         = power_f/noise_p;     %nomalized receive SNR
Pp_em      = power_p_em/noise_p;                  %pilot power
Pp_sp      = power_p_sp/noise_p;                  %pilot power
Pu_em      = power_uemp/noise_p;
Pu_sp      = power_usp/noise_p;

d0 = 0.01;%km
d1 = 0.05;%km


my_date   = date;
if EVA  ==1
    if uncor ==1
        save_file_name = ['UL_v_Ku_EVA_uncor_',my_date,'_N_',num2str(N),'_M_',num2str(M),'_Lpqmax_',num2str(Lpqmax),'_Ku_',num2str(K_u),'.mat'];
    else
        save_file_name = ['UL_v_Ku_EVA_cor_' ,my_date,'_N_',num2str(N),'_M_',num2str(M),'_Lpqmax_',num2str(Lpqmax),'_Ku_',num2str(K_u),'.mat'];
    end
elseif EVB==1
    if uncor ==1
        save_file_name = ['UL_v_Ku_EVB_uncor_',my_date,'_Speed_',num2str(vspeed),'kmph_N_',num2str(N),'_M_',num2str(M),'_Lpqmax_',num2str(Lpqmax),'_Ku_',num2str(K_u),'.mat'];
    else
        save_file_name = ['UL_v_Ku_EVB_cor_' ,my_date,'_Speed_',num2str(vspeed),'kmph_N_',num2str(N),'_M_',num2str(M),'_Lpqmax_',num2str(Lpqmax),'_Ku_',num2str(K_u),'.mat'];
    end
elseif Derrick==1
     if uncor ==1
        save_file_name = ['UL_v_Ku_Dre_uncor_',my_date,'_Speed_',num2str(vspeed),'kmph_N_',num2str(N),'_M_',num2str(M),'_Lpqmax_',num2str(Lpqmax),'_Ku_',num2str(K_u),'.mat'];
    else
        save_file_name = ['UL_v_Ku_Dre_cor_' ,my_date,'_Speed_',num2str(vspeed),'kmph_N_',num2str(N),'_M_',num2str(M),'_Lpqmax_',num2str(Lpqmax),'_Ku_',num2str(K_u),'.mat'];
    end
end

j = sqrt(-1);

R_cf_user_apx_emp_final   = zeros(1,length(K_uv));
R_cf_user_apx_sp_final    = zeros(1,length(K_uv));
R_cf_SE_apx_emp_final     = zeros(1,length(K_uv));
R_cf_SE_apx_sp_final      = zeros(1,length(K_uv));
   
ItR=5000;

for ijk =1:length(K_uv)
    
    ijk
    K_u    = K_uv(ijk);
    
    
    R_cf_user_exact    = zeros(ItR,K_u);
    R_cf_user_apxemp   = zeros(ItR,K_u);
    R_cf_user_apxsp    = zeros(ItR,K_u);
    
    R_cf_SE          = zeros(ItR,1);
    R_cf_SE_apxemp   = zeros(ItR,1);
    R_cf_SE_apxsp    = zeros(ItR,1);
    
    for itr=1:ItR
        
        %%%%%Randomly locations of M_a APs%%%%
        AP=zeros(M_a,2,9);
        AP(:,:,1)=unifrnd(-D/2,D/2,M_a,2);
        
        %Wrapped around (8 neighbor cells)
        D1=zeros(M_a,2);
        D1(:,1)=D1(:,1)+ D*ones(M_a,1);
        AP(:,:,2)=AP(:,:,1)+D1;
        
        D2=zeros(M_a,2);
        D2(:,2)=D2(:,2)+ D*ones(M_a,1);
        AP(:,:,3)=AP(:,:,1)+D2;
        
        D3=zeros(M_a,2);
        D3(:,1)=D3(:,1)- D*ones(M_a,1);
        AP(:,:,4)=AP(:,:,1)+D3;
        
        D4=zeros(M_a,2);
        D4(:,2)=D4(:,2)- D*ones(M_a,1);
        AP(:,:,5)=AP(:,:,1)+D4;
        
        D5=zeros(M_a,2);
        D5(:,1)=D5(:,1)+ D*ones(M_a,1);
        D5(:,2)=D5(:,2)- D*ones(M_a,1);
        AP(:,:,6)=AP(:,:,1)+D5;
        
        D6=zeros(M_a,2);
        D6(:,1)=D6(:,1)- D*ones(M_a,1);
        D6(:,2)=D6(:,2)+ D*ones(M_a,1);
        AP(:,:,7)=AP(:,:,1)+D6;
        
        D7=zeros(M_a,2);
        D7=D7+ D*ones(M_a,2);
        AP(:,:,8)=AP(:,:,1)+D7;
        
        D8=zeros(M_a,2);
        D8=D8- D*ones(M_a,2);
        AP(:,:,9)=AP(:,:,1)+D8;
        
        %Randomly locations of K_u terminals:
        Ter=zeros(K_u,2,9);
        Ter(:,:,1)=unifrnd(-D/2,D/2,K_u,2);
        
        %Wrapped around (8 neighbor cells)
        D1=zeros(K_u,2);
        D1(:,1)=D1(:,1)+ D*ones(K_u,1);
        Ter(:,:,2)=Ter(:,:,1)+D1;
        
        D2=zeros(K_u,2);
        D2(:,2)=D2(:,2)+ D*ones(K_u,1);
        Ter(:,:,3)=Ter(:,:,1)+D2;
        
        D3=zeros(K_u,2);
        D3(:,1)=D3(:,1)- D*ones(K_u,1);
        Ter(:,:,4)=Ter(:,:,1)+D3;
        
        D4=zeros(K_u,2);
        D4(:,2)=D4(:,2)- D*ones(K_u,1);
        Ter(:,:,5)=Ter(:,:,1)+D4;
        
        D5=zeros(K_u,2);
        D5(:,1)=D5(:,1)+ D*ones(K_u,1);
        D5(:,2)=D5(:,2)- D*ones(K_u,1);
        Ter(:,:,6)=Ter(:,:,1)+D5;
        
        D6=zeros(K_u,2);
        D6(:,1)=D6(:,1)- D*ones(K_u,1);
        D6(:,2)=D6(:,2)+ D*ones(K_u,1);
        Ter(:,:,7)=Ter(:,:,1)+D6;
        
        D7=zeros(K_u,2);
        D7=D7+ D*ones(K_u,2);
        Ter(:,:,8)=Ter(:,:,1)+D7;
        
        D8=zeros(K_u,2);
        D8=D8- D*ones(K_u,2);
        Ter(:,:,9)=Ter(:,:,1)+D8;
        
        sigma_shd=8; %in dB
        D_cor=0.1;
        
        %%%%%%Create the M_a x K_u correlated shadowing matrix %%%%%%%
        
        %%%%M correlated shadowing cofficients of M_a APs:
        Dist=zeros(M_a,M_a);%distance matrix
        Cor=zeros(M_a,M_a);%correlation matrix
        
        
        %distance between AP m1 and AP m2
        for m1=1:M_a
            for m2=1:M_a
                Dist(m1,m2) = min([norm(AP(m1,:,1)-AP(m2,:,1)), norm(AP(m1,:,1)-AP(m2,:,2)),...
                    norm(AP(m1,:,1)-AP(m2,:,3)),norm(AP(m1,:,1)-AP(m2,:,4)),norm(AP(m1,:,1)-AP(m2,:,5)),....
                    norm(AP(m1,:,1)-AP(m2,:,6)),norm(AP(m1,:,1)-AP(m2,:,7)),norm(AP(m1,:,1)-AP(m2,:,8)),...
                    norm(AP(m1,:,1)-AP(m2,:,9)) ]);
                Cor(m1,m2)=exp(-log(2)*Dist(m1,m2)/D_cor);
            end
        end
        A1 = chol(Cor,'lower');
        x1 = randn(M_a,1);
        sh_AP = A1*x1;
        for m=1:M_a
            sh_AP(m)=(1/sqrt(2))*sigma_shd*sh_AP(m)/norm(A1(m,:));
        end
        
        %%%K_u correlated shadowing matrix of K_u terminal:
        Dist=zeros(K_u,K_u);%distance matrix
        Cor=zeros(K_u,K_u);%correlation matrix
        
        for k1=1:K_u
            for k2=1:K_u
                Dist(k1,k2)=min([norm(Ter(k1,:,1)-Ter(k2,:,1)), norm(Ter(k1,:,1)-Ter(k2,:,2)),...
                    norm(Ter(k1,:,1)-Ter(k2,:,3)),norm(Ter(k1,:,1)-Ter(k2,:,4)),norm(Ter(k1,:,1)-Ter(k2,:,5)),...
                    norm(Ter(k1,:,1)-Ter(k2,:,6)),norm(Ter(k1,:,1)-Ter(k2,:,7)),norm(Ter(k1,:,1)-Ter(k2,:,8)),...
                    norm(Ter(k1,:,1)-Ter(k2,:,9)) ]); %distance between Terminal k1 and Terminal k2
                Cor(k1,k2)=exp(-log(2)*Dist(k1,k2)/D_cor);
            end
        end
        A2 = chol(Cor,'lower');
        x2 = randn(K_u,1);
        sh_Ter = A2*x2;
        for k=1:K_u
            sh_Ter(k)=(1/sqrt(2))*sigma_shd*sh_Ter(k)/norm(A2(k,:));
        end
        %
        % %%% The shadowing matrix:
        Z_shd=zeros(M_a,K_u);
        for m=1:M_a
            for k=1:K_u
                Z_shd(m,k)= sh_AP(m)+ sh_Ter(k);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Create an M_a x K_u large-scale coefficients beta_mk
        BETAA = zeros(M_a,K_u);
        dist=zeros(M_a,K_u);
        for m=1:M_a
            for k=1:K_u
                [dist(m,k),index] = min([norm(AP(m,:,1)-Ter(k,:,1)), norm(AP(m,:,2)-Ter(k,:,1)),...
                    norm(AP(m,:,3)-Ter(k,:,1)),norm(AP(m,:,4)-Ter(k,:,1)),norm(AP(m,:,5)-Ter(k,:,1)),...
                    norm(AP(m,:,6)-Ter(k,:,1)),norm(AP(m,:,7)-Ter(k,:,1)),norm(AP(m,:,8)-Ter(k,:,1)),...
                    norm(AP(m,:,9)-Ter(k,:,1)) ]); %distance between Terminal k and AP m
                if dist(m,k)<d0
                    betadB=-L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
                elseif ((dist(m,k)>=d0) && (dist(m,k)<=d1))
                    betadB= -L - 35*log10(d1) + 20*log10(d1) - 20*log10(dist(m,k));
                else
                    
                    if uncor ==1
                        betadB = -L - 35*log10(dist(m,k)) + sigma_shd*randn(1,1); %large-scale in dB (uncorrolated)
                    else
                        betadB = -L - 35*log10(dist(m,k)) + Z_shd(m,k); %large-scale in dB (corrolated)
                    end
                    
                end
                
                BETAA(m,k)=10^(betadB/10);
            end
            
        end
        
        Lpq_vec = randi([1,Lpqmax],M_a,K_u);       % Matrix contains nummber of paths between AP m and user k in (m,k)-th entry
        
        %% OTFS channel generation between APs and UEs %%%%
        %=============================================================================
        %%%---------------------- OTFS channel generation-----------------------------
        %=============================================================================
        %[Hpq_mat,ch_gains,ch_delay_indices,ch_Doppler_indices,ch_pDoppler_indices,Lpq_vec] =...
        %OTFS_DD_channel_gen(BETAA,Lpq_vec);
        
        %=============================================================================
        %---- Analytical expressions for X_{q,ij} and Kapa_{q,ij} -(Proposition 1)----
        %=============================================================================
        %[Xpq,Kpq] = OTFS_DD_channel_Xqij(Lpq_vec);
        
        
        %% Create Gamma matrix (variances of the channel estimates)
        Gammaa_emp = zeros(M_a,K_u);  % Embedded-pilot CHE
        Gammaa_sp  = zeros(M_a,K_u);  % Superimposed-pilot CHE
        
        mau=zeros(M_a,K_u);
        nu =zeros(M_a,K_u);
        for m=1:M_a
            for k=1:K_u
                mau(m,k)= (1/N)*sum(  Lpq_vec(m,:).*BETAA(m,:));
                nu(m,k) = ((4*kmax+4*khat+1)./(N^2))* Lpq_vec(m,k)* BETAA(m,k);
            end
        end
        
        for m=1:M_a
            for k=1:K_u
                Gammaa_emp(m,k)=Pp_em*BETAA(m,k)^2/( Pp_em*BETAA(m,k) + Pu_em*( mau(m,k)-nu(m,k))+1);
                Gammaa_sp(m,k) =Pp_sp*BETAA(m,k)^2/( N*Pu_sp*mau(m,k) + N*Pp_sp*mau(m,k) - Pp_sp*Lpq_vec(m,k).*BETAA(m,k)+ 1 );
            end
        end        
        
        %% Each AP has equal power allocations for K terminals
        
        % Compute etaa(m): (each AP transmits equal power to K_u terminals)
        etaa_emp = zeros(M_a,1);
        etaa_sp  = zeros(M_a,1);
        
        for m=1:M_a
            etaa_emp(m) = 1/(sum(Lpq_vec(m,:).*Gammaa_emp(m,:)));
            etaa_sp(m)  = 1/(sum(Lpq_vec(m,:).*Gammaa_sp(m,:)));
        end        
        
        %% Compute Rate
        %SINR=zeros(1,K_u);
        %SINR_exact=zeros(1,K_u);
        %R_cf_exact=zeros(1,K_u);
        
        SINR_apx_emp = zeros(1,K_u);
        SINR_apx_sp  = zeros(1,K_u);
        
        R_cf         = zeros(1,K_u);
        R_cf_apx_emp = zeros(1,K_u);
        R_cf_apx_sp  = zeros(1,K_u);
        
        % Note: assume that gamma_{pq,i} =\gamma_{pq} and \bata_{pq,i}=\bta_{pq}:        
        Gammaa_Lemp = Lpq_vec.* Gammaa_emp;
        Gammaa_Lsp  = Lpq_vec.* Gammaa_sp;
        BETAA_L     = Lpq_vec.* BETAA;
        
        for k=1:K_u            
   
            %denom =0;
            %var_exact =0;
            %Iq1_exact = 0;
            
            num_emp =0;
            Iq2_emp = 0;   
            var_apxemp =0;
            Iq1_apxemp =0;
            noise_emp =0;
            
            num_sp = 0;
            Iq2_sp = 0;   
            var_apxsp =0;
            Iq1_apxsp =0;
            noise_sp  =0;
            
            for m=1:M_a     % m==p and k ==q
                
                % embeded-pilot CHE
                num_emp       = num_emp + Lpq_vec(m,k)*Gammaa_emp(m,k);
                Iq1_apxemp    = Iq1_apxemp + (Lpq_vec(m,k)*BETAA(m,k))* (Gammaa_emp(m,k)*Lpq_vec(m,k));
                Iq2_emp       = Iq2_emp + BETAA_L(m,k)* (sum(Gammaa_Lemp(m,:))-Gammaa_Lemp(m,k));  
                noise_emp     = noise_emp + Gammaa_Lemp(m,k);
                
                % superimposed pilot CHE
                num_sp       = num_sp + Lpq_vec(m,k)*Gammaa_sp(m,k);
                Iq1_apxsp    = Iq1_apxsp + (Lpq_vec(m,k)*BETAA(m,k))* (Gammaa_sp(m,k)*Lpq_vec(m,k));
                Iq2_sp       = Iq2_sp + BETAA_L(m,k)* (sum(Gammaa_Lsp(m,:))-Gammaa_Lsp(m,k)); 
                noise_sp     = noise_sp + Gammaa_Lsp(m,k);
                
            end
            
            SINR_apx_emp(k)   = Pu_em*num_emp^2/(noise_emp + Pu_em*(Iq1_apxemp + Iq2_emp));
            SINR_apx_sp(k)    = Pu_sp*num_sp^2/(noise_sp + Pu_sp*( Iq1_apxsp + Iq2_sp));           
            
            % DL Rate of each user:
            R_cf_apx_emp(k)   = log2(1+ SINR_apx_emp(k));
            R_cf_apx_sp(k)    = log2(1+ SINR_apx_sp(k));
            %R_cf_exact(k)    = log2(1+ SINR_exact(k));            
        end
        
        R_cf_user_apxemp(itr,:)   = R_cf_apx_emp;
        R_cf_user_apxsp(itr,:)    = R_cf_apx_sp;
        
        R_cf_SE_apxemp(itr) = sum(R_cf_apx_emp);
        R_cf_SE_apxsp(itr)  = sum(R_cf_apx_sp);
        
    end
    if K_u> ceil(MN./Nguard)
        R_cf_user_apx_emp_final(ijk) = 0;
        R_cf_user_apx_sp_final(ijk)  = mean(R_cf_user_apxsp(:,1));
    else
    R_cf_user_apx_emp_final(ijk) = mean(R_cf_user_apxemp(:,1));
    R_cf_user_apx_sp_final(ijk)  = mean(R_cf_user_apxsp(:,1));
    end
        
    R_cf_SE_apx_emp_final(ijk) = mean(R_cf_SE_apxemp);
    R_cf_SE_apx_sp_final(ijk)  = mean(R_cf_SE_apxsp);    
    
    save(save_file_name)
end

etta =1;% deltaf*M/10^6;   %Spectral efficieny

omega_emp =1-(MN+Nguard)./(2*MN);
omega_sp = 1-N/(2*N);

figure(1)
hold on
if EVA  ==1
    
     if uncor==0
        plot(K_uv, etta*omega_emp*R_cf_user_apx_emp_final,'-ko', 'LineWidth',1.5) %uncorrelated curves
        hold on
        plot(K_uv, etta*omega_sp*R_cf_user_apx_sp_final,'-ro', 'LineWidth',1.5) %uncorrelated curves
    else
        plot(K_uv, etta*omega_emp*R_cf_user_apx_emp_final,'-k*', 'LineWidth',1.5)%correlated curves
        hold on
        plot(K_uv, etta*omega_sp*R_cf_user_apx_sp_final,'-r*', 'LineWidth',1.5)%correlated curves
    end  
    
elseif EVB  ==1
    
    if uncor==0
        plot(K_uv, etta*omega_emp*R_cf_user_apx_emp_final,'-bo', 'LineWidth',1.5) %uncorrelated curves
        hold on
        plot(K_uv, etta*omega_sp*R_cf_user_apx_sp_final,'-go', 'LineWidth',1.5) %uncorrelated curves
    else
        plot(K_uv, etta*omega_emp*R_cf_user_apx_emp_final,'-b*', 'LineWidth',1.5)%correlated curves
        hold on
        plot(K_uv, etta*omega_sp*R_cf_user_apx_sp_final,'-g*', 'LineWidth',1.5)%correlated curves
    end  
    
elseif Derrick==1
    
    if uncor==0
        plot(K_uv, etta*omega_emp*R_cf_user_apx_emp_final,'-ko', 'LineWidth',1.5) %uncorrelated curves
        hold on
        plot(K_uv, etta*omega_sp*R_cf_user_apx_sp_final,'-ro', 'LineWidth',1.5) %uncorrelated curves
    else
        plot(K_uv, etta*omega_emp*R_cf_user_apx_emp_final,'-k*', 'LineWidth',1.5)%correlated curves
        hold on
        plot(K_uv, etta*omega_sp*R_cf_user_apx_sp_final,'-r*', 'LineWidth',1.5)%correlated curves
    end  
    
end

%'MarkerSize',6
box on
grid on
xlim([K_uv(1), K_uv(end)])
h =gca;
set(h, 'FontSize',12);
hx=xlabel('User transmit power (W)');
set(hx, 'FontSize',12,'FontName','Arial');
hy=ylabel('Average Uplink SE (Mbits/s)');
set(hy, 'FontSize',12,'FontName','Arial');



