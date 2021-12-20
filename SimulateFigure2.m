clc

%%This Matlab script generates Figure 2 in the paper:
%
%Mohammadali Mohammadi, Hien Quoc Ngo and Michail Matthaiou, "Cell-Free Massive MIMO Meets OTFS
%Modulation," submitted in IEEE Transactions on Communications
%
%Download article: https://www.researchgate.net/publication/357158978_Cell-Free_Massive_MIMO_Meets_OTFS_Modulation
%
%This is version 1.0 (Last edited: 2021-12-18)
%
%License: This code is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our
%paper as described above.
%==========================================================================

global  M_a K_u M N deltaf elmax kmax Lpqmax Pp_sp Pp_em Pd

%% Inital parameters
EVA   = 0;  %Set EVA =1 to model extended vehicular A model
EVB   = 0;  %Set EVB =1 to model extended vehicular B model
TWC   = 1;
%==========================================================================
%==========================================================================
%%%----------------------- OTFS parameters --------------------------------
%==========================================================================
%==========================================================================
% We consider an OTFS system with $\Delta f=15$ kHz, carrier frequency $f_c=4$ GHz
% and vehicular speed of  $300$ kmph. We consider two 3GPP vehicular models, namely
% (i) extended vehicular A (EVA) with $\ell_{pq}=9$ and  $\tau_{max}=2.51~\mu sec$ and
% (ii) extended vehicular B (EVB) with  $\ell_{pq}=6$ and  $\tau_{max}=10~\mu sec$.
% The block duration is assumed  as $N=128$ which corresponds to $T_f=8.85~msec$,
% and the number of sub-carriers is $M=512$.

clight  = 3*10^8;         % speed of light in m/s
deltaf  = 15*10^3;        % Frequency spacing between adjacent sub-carriers
f       = 4000;           % Carrier frequency in MHz 
if EVA==1
    %channel parameters for vehicular speed of up to 300 kmph
    vspeed  = 300;            % vehicular speed in kmph
    N  = 128;                 % number of symbol N  = 128;
    M  = 512;                 % number of subcarriers
    MN = N*M;                 % number of symbols per frame
    
    taumax = 2.5*10^-6;        %maximum delay spreed according to EVA,
    elmax  = floor(taumax*deltaf*M);
    
    numax  = f *vspeed ./(3.6*clight/1000000);  %f is in kHz
    kmax   = floor (numax*N/deltaf);
    khat    = 3;   % reduced guard interval length
    
    Lpqmax  = 9;   % maximum number of path between the APs and uers
elseif EVB==1    
    vspeed  = 300;           % vehicular speed in kmph    
    N  = 120;                % number of symbol
    M  = 512;                % number of subcarriers
    MN = N*M;                % number of symbols per frame
    
    taumax = 20*10^-6;       %maximum delay spreed according to EVA,
    elmax  = floor(taumax*deltaf*M);
    
    numax  = f *vspeed ./(3.6*clight./1000000);  %f is in kHz
    kmax   = floor (numax*N/deltaf);
    khat   = 3;              % reduced guard interval length    
    Lpqmax = 6;              % maximum number of path between the APs and uers
    
elseif TWC==1
    
    %  TWC.2021.Performance Analysis of Coded OTFS Systems over High-Mobility Channels
    %  Integer delay and Doppler case and set the maximum delay index
    %  as lmax = 3 and the maximum Doppler index as kmax = 5, which is
    %  corresponding to a relative speed around 250 km/h with 4 GHz
    %  carrier frequency and 1.5 kHz sub-carrier spacing.
    
    %channel parameters for vehicular speed of up to 250 kmph
    vspeed  = 250;             % vehicular speed in kmph
    f       = 4000;            % Carrier frequency in MHz (1900)
    tumax   = 2.5*10^-6;
    
    N  = 4;                   % number of symbol
    M  = 9;                   % number of subcarriers
    MN = N*M;                  % number of symbols per frame
    
    numax  = f *vspeed ./(3.6*clight/1000000);  %f is in MHz
    kmax   = floor (numax*N/deltaf);
    
    elmax   = floor (tumax*M*deltaf);   % maximum delay index
    khat    = 1;                        % reduced guard interval length
    Lpqmax  = 3;                        % maximum number of path between the APs and uers
end

% For each channel realization, we randomly select the delay and Doppler
% indices according to the uniform distribution, such that we have
% -kmax<ki<kmax and 0<li<lmax

Nguard = (4*kmax+4*khat)*(2*elmax+1);
%==========================================================================
%==========================================================================
%%%-------------------  Cell free Network parameters ----------------------
%==========================================================================
%==========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Downplink
% Consider all bettas are set to one. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ma_v = [10:25:110, 250]; %Number of APs
K_u  = 10;                       %Number of r useterminals

power_f     = 15;      % Normalized transmit SNR in Downlink (dB)
power_u     = 10;      % Normalized transmit SNR in Uplink (dB)
power_p_em  = 15;      % Normalized transmit pilot SNR (emb) in (dB)
power_p_sp  = power_u; % Normalized transmit pilot SNR (sp) in (dB)
noise_p     = 1;       %noise power

Ed       = 10^(power_f/10);     %nomalized receive SNR
Pp_em    = 10^(power_p_em/10); 
Pp_sp    = 10^(power_p_sp/10);   
Pu       = 10^(power_u/10);

%Evaluating gammq for the scalling laws
c_1      = K_u - (4*kmax+4*khat+1)./N;
Lpq      = Lpqmax;
gamq_ep  = 0.9*Pp_em./ (1+Pp_em + Pu *Lpqmax *c_1 ./N);
gamq_sp  = 0.9*Pp_sp./ (1+Pp_sp + Pp_sp*Lpqmax*(K_u-1) + Pu *Lpqmax *K_u);


my_date   = date;
if EVA  ==1
    save_file_name = ['Sim_EVA_DLSE_v_MA_betta1_',my_date,'_N_',num2str(N),'_M_',num2str(M),'_Lpqmax_',num2str(Lpqmax),'_Ku_',num2str(K_u),'_Pd_',num2str(power_f),'.mat'];
elseif EVB==1
    save_file_name = ['Sim_EVB_DLSE_v_MA_betta1_',my_date,'_Speed_',num2str(vspeed),'kmph_N_',num2str(N),'_M_',num2str(M),'_Lpqmax_',num2str(Lpqmax),'_Ku_',num2str(K_u),'_Pd_',num2str(power_f),'.mat'];
elseif TWC==1
    save_file_name = ['Sim_TWC_DLSE_v_MA_betta1_',my_date,'_Speed_',num2str(vspeed),'kmph_N_',num2str(N),'_M_',num2str(M),'_Lpqmax_',num2str(Lpqmax),'_Ku_',num2str(K_u),'_Pd_',num2str(power_f),'.mat'];
end

j = sqrt(-1);

R_cf_user_emp_final   = zeros(length(Ma_v),K_u);
R_cf_user_sp_final    = zeros(length(Ma_v),K_u);
R_cf_SE_emp_final     = zeros(1,length(Ma_v));
R_cf_SE_sp_final      = zeros(1,length(Ma_v));
SINR_sim_em           = zeros(1,length(Ma_v));
SINR_sim_sp           = zeros(1,length(Ma_v));
R_cf_Slaw_emp         = zeros(1,length(Ma_v));
R_cf_Slaw_sp          = zeros(1,length(Ma_v));
R_cf_UP_sp            = zeros(1,length(Ma_v));
R_cf_UP_ep            = zeros(1,length(Ma_v));

ItR =250;

for ijk = 1:length(Ma_v)
    
    M_a     = Ma_v(ijk)
    Pd      = (Ed)./M_a^2;     
    
    % Matrix contains nummber of paths between AP m and user k in (m,k)-th entry
    Lpq_vec = randi([1,Lpqmax],M_a,K_u);
    BETAA   = ones(M_a,K_u);
    
    %Initializing the parameters
    R_cf_user_emp   = zeros(ItR,K_u);
    R_cf_user_sp    = zeros(ItR,K_u);
    R_cf_SE         = zeros(ItR,1);
    R_cf_SE_emp     = zeros(ItR,1);
    R_cf_SE_sp      = zeros(ItR,1);
    
    Iq1_emk = zeros(ItR,K_u);
    Iq1_spk = zeros(ItR,K_u);
    Iq2_emk = zeros(ItR,1);
    Iq2_spk = zeros(ItR,1);
    DS1_emk = zeros(ItR,K_u);
    DS1_spk = zeros(ItR,K_u);
    
    X_ep = zeros(M*N,M*N);
    X_sp = zeros(M*N,M*N);
    Y_ep = zeros(M*N,M*N);
    Y_sp = zeros(M*N,M*N);
    
    for itr=1:ItR
        
        
        %% Create Gamma matrix parameters (variances of the channel estimates and errors)
        Gammaa_emp = zeros(M_a,K_u);  % Embedded-pilot CHE
        Gammaa_sp  = zeros(M_a,K_u);  % Superimposed-pilot CHE
        c_emp      = zeros(M_a,K_u);  % Embedded-pilot CHE
        c_sp       = zeros(M_a,K_u);  % Embedded-pilot CHE
        er_emp     = zeros(M_a,K_u);  % estimation error variance for the channel estimate realizations
        er_sp      = zeros(M_a,K_u);  % estimation error variance for the channel estimate realizations
        
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
                
                Gammaa_emp(m,k) = Pp_em*BETAA(m,k)^2/( Pp_em*BETAA(m,k) + Pu*( mau(m,k)-nu(m,k))+1);
                c_emp(m,k)      = sqrt(Pp_em)*BETAA(m,k)/( Pp_em*BETAA(m,k) + Pu*( mau(m,k)-nu(m,k))+1);
                er_emp(m,k)     = Pu*( mau(m,k)-nu(m,k))+1;
                
                Gammaa_sp(m,k)   = Pp_sp*BETAA(m,k)^2/(Pp_sp*BETAA(m,k)+ N*Pu*mau(m,k) + N*Pp_sp*mau(m,k)- Pp_sp*Lpq_vec(m,k).*BETAA(m,k)+ 1 );
                c_sp(m,k)        = sqrt(Pp_sp)*BETAA(m,k)/(Pp_sp*BETAA(m,k)+ N*Pu*mau(m,k) + N*Pp_sp*mau(m,k)- Pp_sp*Lpq_vec(m,k).*BETAA(m,k)+ 1 );
                er_sp(m,k)       = N*Pu*mau(m,k) + N*Pp_sp*mau(m,k)- Pp_sp*Lpq_vec(m,k).*BETAA(m,k)+ 1;                
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
        
        %% OTFS channel generation between APs and UEs %%%%
        %=============================================================================
        %%%----------------- OTFS channel generation and Simulation-------------------
        %=============================================================================
        [DS1_emk(itr,:),DS1_spk(itr,:),Iq1_emk(itr,:),Iq1_spk(itr,:),Iq2_emk(itr),...
            Iq2_spk(itr),D_qq_ep,D_qq_sp,sum_ep,sum_sp] =...
            DL_rate_sim_PCSI(BETAA,Lpq_vec,c_emp,c_sp,er_emp,er_sp,etaa_emp,etaa_sp);   
        
        %% Upper  bound parameters        
        X_ep = X_ep +D_qq_ep;
        X_sp = X_sp +D_qq_sp;
        
        Y_ep = Y_ep +sum_ep;
        Y_sp = Y_sp +sum_sp;
        
        %% Compute Rate
        SINR_emp = zeros(1,K_u);
        SINR_sp  = zeros(1,K_u);
        
        R_cf     = zeros(1,K_u);
        R_cf_emp = zeros(1,K_u);
        R_cf_sp  = zeros(1,K_u);
        
        % Note: assume that gamma_{pq,i} =\gamma_{pq} and \bata_{pq,i}=\bta_{pq}:
        Gammaa_Lemp = Lpq_vec.* Gammaa_emp;
        Gammaa_Lsp  = Lpq_vec.* Gammaa_sp;
        BETAA_L     = Lpq_vec.* BETAA;
        
        for k=1:K_u
            
            num_emp =0;
            Iq2_emp = 0;
            var_emp =0;
            Iq1_emp =0;
            
            num_sp =0;
            Iq2_sp = 0;
            var_sp =0;
            Iq1_sp =0;
            
            for m=1:M_a                
                % embeded-pilot CHE
                num_emp       = num_emp + (etaa_emp(m)^(1/2))*Lpq_vec(m,k)*Gammaa_emp(m,k);
                var_emp       = var_emp + etaa_emp(m)*BETAA(m,k)* Gammaa_emp(m,k)*Lpq_vec(m,k);
                Iq1_emp       = Iq1_emp + etaa_emp(m)*BETAA(m,k)* Gammaa_emp(m,k)*Lpq_vec(m,k)*(Lpq_vec(m,k)-1);
                Iq2_emp       = Iq2_emp + etaa_emp(m)*BETAA_L(m,k)* (sum(Gammaa_Lemp(m,:))-Gammaa_Lemp(m,k));
                
                % superimposed pilot CHE
                num_sp       = num_sp + (etaa_sp(m)^(1/2))*Lpq_vec(m,k)*Gammaa_sp(m,k);
                var_sp       = var_sp + etaa_sp(m)*BETAA(m,k)* Gammaa_sp(m,k)*Lpq_vec(m,k);
                Iq1_sp       = Iq1_sp + etaa_sp(m)*BETAA(m,k)* Gammaa_sp(m,k)*Lpq_vec(m,k)*(Lpq_vec(m,k)-1);
                Iq2_sp       = Iq2_sp + etaa_sp(m)*BETAA_L(m,k)* (sum(Gammaa_Lsp(m,:))-Gammaa_Lsp(m,k));
                
            end
            
            SINR_emp(k)   = Pd*num_emp^2/(1 + Pd*(var_emp + Iq1_emp + Iq2_emp));
            SINR_sp(k)    = Pd*num_sp^2/(1 + Pd*(var_sp + Iq1_sp + Iq2_sp));
            
            % DL Rate of each user based on Theorem 1:
            R_cf_emp(k)   = log2(1+ SINR_emp(k));
            R_cf_sp(k)    = log2(1+ SINR_sp(k));
            
        end
        
        R_cf_user_emp(itr,:)   = R_cf_emp;
        R_cf_user_sp(itr,:)    = R_cf_sp;
        
        R_cf_SE_emp(itr) = sum(R_cf_emp);
        R_cf_SE_sp(itr)  = sum(R_cf_sp);        
    end
    
    % Upper bound for embedded and superimposed pilot-based channel
    % estimations methods
    
    X_bar_ep  = X_ep./ItR;
    X_bar_sp  = X_sp./ItR;
    
    Y_bar_ep  = Y_ep./ItR;
    Y_bar_sp  = Y_sp./ItR;
    
    Si_ep  = eye(M*N) + Pd*Y_bar_ep-Pd*(X_bar_ep*X_bar_ep');
    Si_sp  = eye(M*N) + Pd*Y_bar_sp-Pd*(X_bar_sp*X_bar_sp');
    
    R_cf_UP_ep(ijk) = (1/MN)*log2(real(det(eye(M*N) + Pd*X_bar_ep'*(Si_ep^(-1))*X_bar_ep)));
    R_cf_UP_sp(ijk) = (1/MN)*log2(real(det(eye(M*N) + Pd*X_bar_sp'*(Si_sp^(-1))*X_bar_sp)));    
    
    R_cf_user_emp_final(ijk,:) = mean(R_cf_user_emp);
    R_cf_user_sp_final(ijk,:)  = mean(R_cf_user_sp);
    
    R_cf_SE_emp_final(ijk) = mean(R_cf_SE_emp);
    R_cf_SE_sp_final(ijk)  = mean(R_cf_SE_sp);
    
    %  Simulation Results for embedded-pilot based channel estimation
    v1bar_em     = mean(DS1_emk(:,1));
    nom_em       = Pd*abs(v1bar_em)^2;
    denom1_em    = Pd*mean(abs(DS1_emk(:,1) - v1bar_em).^2);
    denom2_em    = Pd*mean(Iq1_emk(:,1));
    denom3_em    = Pd*mean(Iq2_emk);    
    SINR_sim_em(ijk)  = nom_em./ (denom1_em+denom2_em+denom3_em+1);
    
    %  Simulation Results for superimposed pilot-based channel estimation
    v1bar_sp     = mean(DS1_spk(:,1));
    nom_sp       = Pd*abs(v1bar_sp)^2;
    denom1_sp    = Pd*mean(abs(DS1_spk(:,1) - v1bar_sp).^2);
    denom2_sp    = Pd*mean(Iq1_spk(:,1));
    denom3_sp    = Pd*mean(Iq2_spk);
    SINR_sim_sp(ijk)  = nom_sp./ (denom1_sp+denom2_sp+denom3_sp+1);
    
    % Scalling Law curves for both channel estimation method
    R_cf_Slaw_emp(ijk) = log2(1 + (Ed*Lpq*gamq_ep)./K_u);
    R_cf_Slaw_sp(ijk)  = log2(1 + (Ed*Lpq*gamq_sp)./K_u);
    
    save(save_file_name)
end

R_cf_user_sim_emp_final = log2(1+SINR_sim_em);
R_cf_user_sim_sp_final  = log2(1+SINR_sim_sp);

etta =1;% deltaf*M/10^6;   %Spectral efficieny

figure(1)
hold on
if EVA  ==1
    plot(Ma_v, etta*R_cf_user_emp_final(:,1),'-k*', 'LineWidth',1.5)%correlated curves
    hold on
    plot(Ma_v, etta*R_cf_user_sp_final(:,1),'-ko', 'LineWidth',1.5)%correlated curves
    hold on
    plot(Ma_v, etta*R_cf_user_sim_emp_final,'-.k*', 'LineWidth',1.5) %uncorrelated curves
    hold on
    plot(Ma_v, etta*R_cf_user_sim_sp_final,'-.ko', 'LineWidth',1.5) %uncorrelated curves
    hold on
    plot(Ma_v, etta*R_cf_Slaw_emp,'-.r', 'LineWidth',1.5) %Scalling low curves
    hold on
    plot(Ma_v, etta*R_cf_Slaw_sp,'-.r', 'LineWidth',1.5) %Scalling low curves
    hold on
    plot(Ma_v, etta*R_cf_UP_ep,'-.b', 'LineWidth',1.5) %Upperbound curves
    hold on
    plot(Ma_v, etta*R_cf_UP_sp,'-.b', 'LineWidth',1.5) %Upperbound curves
    
elseif EVB  ==1
    
    plot(Ma_v, etta*R_cf_user_emp_final(:,1),'-k*', 'LineWidth',1.5)%correlated curves
    hold on
    plot(Ma_v, etta*R_cf_user_sp_final(:,1),'-ko', 'LineWidth',1.5)%correlated curves
    hold on
    plot(Ma_v, etta*R_cf_user_sim_emp_final,'-.k*', 'LineWidth',1.5) %uncorrelated curves
    hold on
    plot(Ma_v, etta*R_cf_user_sim_sp_final,'-.ko', 'LineWidth',1.5) %uncorrelated curves
    hold on
    plot(Ma_v, etta*R_cf_Slaw_emp,'-.r', 'LineWidth',1.5) %Scalling low curves
    hold on
    plot(Ma_v, etta*R_cf_Slaw_sp,'-.r', 'LineWidth',1.5) %Scalling low curves
    hold on
    plot(Ma_v, etta*R_cf_UP_ep,'-.b', 'LineWidth',1.5) %Upperbound curves
    hold on
    plot(Ma_v, etta*R_cf_UP_sp,'-.b', 'LineWidth',1.5) %Upperbound curves
    
elseif TWC==1
    
    plot(Ma_v, etta*R_cf_user_emp_final(:,1),'-k*', 'LineWidth',1.5)%correlated curves
    hold on
    plot(Ma_v, etta*R_cf_user_sp_final(:,1),'-ko', 'LineWidth',1.5)%correlated curves
    hold on
    plot(Ma_v, etta*R_cf_user_sim_emp_final,'-.k*', 'LineWidth',1.5) %uncorrelated curves
    hold on
    plot(Ma_v, etta*R_cf_user_sim_sp_final,'-.ko', 'LineWidth',1.5) %uncorrelated curves
    hold on
    plot(Ma_v, etta*R_cf_Slaw_emp,'-.r', 'LineWidth',1.5) %Scalling low curves
    hold on
    plot(Ma_v, etta*R_cf_Slaw_sp,'-.r', 'LineWidth',1.5) %Scalling low curves
    hold on
    plot(Ma_v, etta*R_cf_UP_ep,'-.b', 'LineWidth',1.5) %Upperbound curves
    hold on
    plot(Ma_v, etta*R_cf_UP_sp,'-.b', 'LineWidth',1.5) %Upperbound curves
    
end

%'MarkerSize',6
box on
grid on
h =gca;
set(h, 'FontSize',12);
hx=xlabel('Number of APs (M_a)');
set(hx, 'FontSize',12,'FontName','Arial');
hy=ylabel('Average Downlink SE (Mbis/s/Hz)');
set(hy, 'FontSize',12,'FontName','Arial');
