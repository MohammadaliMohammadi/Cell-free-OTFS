


function [DS_emp,DS_sp,Iq1_emk,Iq1_spk,Iq2_emk,Iq2_spk,D_qq_ep,D_qq_sp,Y_ep,Y_sp] ...
    = DL_rate_sim_PCSI(BETAA,Lpq_vec,c_emp,c_sp,er_emp,er_sp,etaa_emp,etaa_sp)

global M_a K_u M N  elmax kmax Pp_sp Pp_em 

PI0 = eye(M*N);
PI  = circshift(PI0,1);

z     = exp(sqrt(-1)*2*pi/(M*N));
Delta = diag(z.^(0:1:M*N-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% DD channel (MN x MN) %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hpq_mat     = zeros(M*N,M*N,M_a); %cell(M_a,K_u);

DS_emp   = zeros(1,K_u);
DS_sp    = zeros(1,K_u);
Iq1_em1k = zeros(1,K_u);
Iq1_sp1k = zeros(1,K_u);

Iq2_em1k =0;
Iq2_sp1k =0;

D_qq_ep   = zeros(M*N,M*N);
D_qq_sp   = zeros(M*N,M*N);
Y_ep      = zeros(M*N,M*N);
Y_sp      = zeros(M*N,M*N);

for k=1:K_u
    
    D_qqr_ep  = zeros(M*N,M*N);
    D_qqr_sp  = zeros(M*N,M*N);

    for m=1:M_a
        
        Hpq     = zeros(M*N,M*N);
        hHpq_ep = zeros(M*N,M*N);
        hHpq_sp = zeros(M*N,M*N);
        Lpq     = Lpq_vec(m,k);
        
        ltui     = randi([0,elmax],1,Lpq);      %[3,4,1];   % delay indices in [1:M]
        ktui     = randi([-kmax,kmax],1,Lpq);   %[3,0,1];   % Doppler indices in   [1:N]
        kaptui   = -0.5+rand(1,Lpq);            %[0.2,0,0]; % fractional Doppler shifts in [-0.5, 05];
        
          
        %channel coefficeints
        hpq = sqrt(BETAA(m,k))*( sqrt(1/2)*randn(1,Lpq)+1i*sqrt(1/2)*randn(1,Lpq));
        
        %Channel estimates EM
        hpq_ep = c_emp(m,k)*sqrt(Pp_em)*hpq + ...
            c_emp(m,k)*(sqrt(er_emp(m,k)))*((sqrt(1/2)*randn(1,Lpq)+1i*sqrt(1/2)*randn(1,Lpq)));
        
        %Channel estimates SP
        hpq_sp = c_sp(m,k)*sqrt(Pp_sp)*hpq+ ...
            c_sp(m,k)*sqrt(er_sp(m,k))*((sqrt(1/2)*randn(1,Lpq)+1i*sqrt(1/2)*randn(1,Lpq)));
        
        for ii=1:Lpq
            Dmat = (PI^(ltui(ii)))*(Delta^(ktui(ii)+kaptui(ii)));
            Tii  = kron(1/sqrt(N)*dftmtx(N),eye(M))* Dmat* kron((1/sqrt(N)*dftmtx(N))',eye(M));
            
            Hpq     = Hpq + hpq(ii).*Tii;
            hHpq_ep = hHpq_ep + hpq_ep(ii).*Tii;
            hHpq_sp = hHpq_sp + hpq_sp(ii).*Tii;
        end
        
        hHpq_ep =hHpq_ep';
        hHpq_sp =hHpq_sp';
        
        DS_emp(k)       = DS_emp(k) + (etaa_emp(m))^(1/2)*Hpq(1,:)* hHpq_ep(:,1); 
        DS_sp(k)        = DS_sp(k) + (etaa_sp(m))^(1/2)*Hpq(1,:)* hHpq_sp(:,1);          
       
        Iq1_em1k(k)      = Iq1_em1k(k) + (etaa_emp(m))^(1/2)*sum(Hpq(1,:)* hHpq_ep(:,2:end));
        Iq1_sp1k(k)      = Iq1_sp1k(k) + (etaa_sp(m))^(1/2)*sum(Hpq(1,:)* hHpq_sp(:,2:end));
        
        if k==1            
            % Here we keep the Hp1 for next step to calculate the Iq2
            Hpq_mat(:,:,m)     = Hpq;
            D_qq_ep    = D_qq_ep + (etaa_emp(m))^(1/2)* Hpq*hHpq_ep;
            D_qq_sp    = D_qq_sp + (etaa_sp(m))^(1/2) * Hpq*hHpq_sp;           

        else
            % Inter-user interference is calculated only for the first user
            Iq2_em1k = Iq2_em1k + (etaa_emp(m))^(1/2)*sum(Hpq_mat(1,:,m)* hHpq_ep(:,1:end));
            Iq2_sp1k = Iq2_sp1k + (etaa_sp(m))^(1/2)*sum(Hpq_mat(1,:,m)* hHpq_sp(:,1:end));
            
            D_qqr_ep    = D_qqr_ep + (etaa_emp(m))^(1/2)* Hpq_mat(:,:,m)*hHpq_ep;
            D_qqr_sp    = D_qqr_sp + (etaa_sp(m))^(1/2) * Hpq_mat(:,:,m)*hHpq_sp;
        end
    end
    
    Y_ep  = Y_ep + D_qqr_ep *D_qqr_ep';
    Y_sp  = Y_sp + D_qqr_sp *D_qqr_sp';    
    
end
Iq1_emk = abs(Iq1_em1k).^2;
Iq1_spk = abs(Iq1_sp1k).^2;

Iq2_emk = abs(Iq2_em1k).^2;
Iq2_spk = abs(Iq2_sp1k).^2;

Y_ep  = (Y_ep + D_qq_ep*D_qq_ep'); 
Y_sp  = (Y_sp + D_qq_sp*D_qq_sp'); 

end