function [AD2,AD1,I,E1,E2,K1,KD2] = Cheb_matrices_symbolic(N,H,a,t,alpha,chi,beta)

f = pi;
cutoff = 0.0000001;

%%% Gauss-Chebyshev points %%%%%%%%%
K = 2000;
zeta = cos((2*[1:K]-1)*pi/(2*(K)));
alph = pi/(K);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% identity <W_p, V_k> %%%%%
I = zeros(N+1,N+1);

for k = 0:N
    for p = 0:N
        
        if k == p+2
            I(k+1,p+1) = -c(p)/4;
            
        elseif k == p
            I(k+1,p+1) = ( 2 + 2*c(p) - c(p-1) )/4;
            
        elseif k == p-2
            I(k+1,p+1) = -( 2 + c(p-2) )/4;
            
        elseif k == p-4
            I(k+1,p+1) = 1/4;
            
        else
            I(k+1,p+1) = 0;
        end
        
    end
end

I = (f/2)*I;


%%%%% matrix E1 <exp(-(1/4)*H*(zeta+1)/a^gamma)*W_p, V_k> %%%%%

fid = fopen('matrixE1.mat','r');
if fid == -1
    E1   = zeros(N+1,N+1);
    for k = 0:N
        for p = 0:N
            z = H*(zeta+1)/(2*a);
            val = alph*sum(exp(-alpha*z/2).*w(p,zeta).*v(k,zeta));
            E1(k+1,p+1) = val;
        end
    end
    Hsav = H;
    Nsav = N;
    asav = a;
    alphasav = alpha;
    tsav = t;
    E1sav = E1;
    save matrixE1 E1sav Hsav Nsav asav alphasav tsav
else
    fclose(fid);
    load matrixE1 E1sav Hsav Nsav asav alphasav tsav
    if ( (H == Hsav)  & (Nsav >= N) & (a == asav) & (alphasav == alpha) & (t == tsav))
        E1   = E1sav(1:N+1,1:N+1);
    else
        E1   = zeros(N+1,N+1);
        for k = 0:N
            for p = 0:N
                z = H*(zeta+1)/(2*a);
                val = alph*sum(exp(-alpha*z/2).*w(p,zeta).*v(k,zeta));
                E1(k+1,p+1) = val;
            end
        end
        Hsav = H;
        Nsav = N;
        asav = a;
        alphasav = alpha;
        tsav = t;
        E1sav = E1;
        save matrixE1 E1sav Hsav Nsav asav alphasav tsav
    end
end

%%%%% matrix E2 <exp((1/2)*H*((1/2)*alpha-1)*(zeta+1)/a^gamma)*W_p, V_k> %%%%%

fid = fopen('matrixE2.mat','r');
if fid == -1
    E2   = zeros(N+1,N+1);
    for k = 0:N
        for p = 0:N
            z = H*(zeta+1)/(2*a);
            val = alph*sum(exp(alpha*z/2).*primprof(z,t,chi).*w(p,zeta).*v(k,zeta));
            E2(k+1,p+1) = val;
        end
    end
    Hsav = H;
    Nsav = N;
    asav = a;
    alphasav = alpha;
    tsav = t;
    E2sav = E2;
    save matrixE2 E2sav Hsav Nsav asav alphasav tsav
else
    fclose(fid);
    load matrixE2 E2sav Hsav Nsav asav alphasav tsav
    if ( (H == Hsav)  & (Nsav >= N) & (a == asav) & (alphasav == alpha) & (t == tsav))
        E2   = E2sav(1:N+1,1:N+1);
    else
        E2   = zeros(N+1,N+1);
        for k = 0:N
            for p = 0:N
                z = H*(zeta+1)/(2*a);
                val = alph*sum(exp(alpha*z/2).*primprof(z,t,chi).*w(p,zeta).*v(k,zeta));
                E2(k+1,p+1) = val;
            end
        end
        Hsav = H;
        Nsav = N;
        asav = a;
        alphasav = alpha;
        tsav = t;
        E2sav = E2;
        save matrixE2 E2sav Hsav Nsav asav alphasav tsav
    end
end


fid = fopen('matrixK1.mat','r');
if fid == -1
  K1   = zeros(N+1,N+1);
  for k = 0:N
    for p = 0:N
        z = H*(zeta+1)/(2*a);
        val = alph*sum((1./permeability(z, beta)).*w(p,zeta).*v(k,zeta));
        K1(k+1,p+1) = val;
    end
  end
  Hsav = H;
  Nsav = N;
  asav = a;
  alphasav = alpha;
  tsav = t;
  betasav = beta;
  K1sav = K1;
  save matrixK1 K1sav Hsav Nsav asav alphasav tsav betasav
else
  fclose(fid);
  load matrixK1 K1sav Hsav Nsav asav alphasav tsav betasav
  if ( (H == Hsav)  & (Nsav >= N) & (a == asav) & (alphasav == alpha) & (t == tsav) & (beta == betasav))    
    K1 = K1sav(1:N+1,1:N+1);
  else    
    K1   = zeros(N+1,N+1);
    for k = 0:N
        for p = 0:N
          z = H*(zeta+1)/(2*a);
          val = alph*sum((1./permeability(z, beta)).*w(p,zeta).*v(k,zeta));
          K1(k+1,p+1) = val;
        end
    end
    Hsav = H;
    Nsav = N;
    asav = a;
    alphasav = alpha;
    tsav = t;
    betasav = beta;
    K1sav = K1;
    save matrixK1 K1sav Hsav Nsav asav alphasav tsav betasav
  end
end

% fid = fopen('matrixK2.mat','r');
% if fid == -1
%
%   disp('Creating file...')
%   K2   = zeros(N+1,N+1);
%   for k = 0:N
%     for p = 0:N
%         z = H*(zeta+1)/(2*a);
%         val = alph*sum((-permeability_prime(z, beta)./permeability(z, beta)).*Dw(p,zeta).*v(k,zeta));
%         if abs(val) < cutoff
%             K2(k+1,p+1) = 0;
%         else
%             K2(k+1,p+1) = (2/H)*val;
%         end
%     end
%   end
%   Hsav = H;
%   Nsav = N;
%   asav = a;
%   alphasav = alpha;
%   tsav = t;
%   betasav = beta;
%   K2sav = K2;
%   save matrixK2 K2sav Hsav Nsav asav alphasav tsav betasav
% else
%   fclose(fid);
%   load matrixK2 K2sav Hsav Nsav asav alphasav tsav betasav
%   if ( (H == Hsav)  & (Nsav >= N) & (a == asav) & (alphasav == alpha) & (t == tsav) & (beta == betasav))
%     disp('Loading existing file...');
%     K2 = K2sav(1:N+1,1:N+1);
%   else
%     %disp('Overwriting file...')
%     K2   = zeros(N+1,N+1);
%     for k = 0:N
%         for p = 0:N
%           z = H*(zeta+1)/(2*a);
%           val = alph*sum((-permeability_prime(z, beta)./permeability(z, beta)).*Dw(p,zeta).*v(k,zeta));
%           if abs(val) < cutoff
%               K2(k+1,p+1) = 0;
%           else
%               K2(k+1,p+1) = (2/H)*val;
%           end
%         end
%     end
%     Hsav = H;
%     Nsav = N;
%     asav = a;
%     alphasav = alpha;
%     tsav = t;
%     betasav = beta;
%     K2sav = K2;
%     save matrixK2 K2sav Hsav Nsav asav alphasav tsav betasav
%   end
% end

fid = fopen('matrixKD2.mat','r');
if fid == -1
    KD2   = zeros(N+1,N+1);
    for k = 0:N
        for p = 0:N
            z = H*(zeta+1)/(2*a);
            val = alph*sum(-(1./permeability(z, beta)).*Dw(p,zeta).*Dv(k,zeta));
            val = val + alph*sum(-(1./permeability(z, beta)).*(zeta./(1-zeta.^2)).*Dw(p,zeta).*v(k,zeta));
            if abs(val) < cutoff
                KD2(k+1,p+1) = 0;
            else
                KD2(k+1,p+1) = (4/H^2)*val;
            end
        end
    end
    Hsav = H;
    Nsav = N;
    asav = a;
    alphasav = alpha;
    tsav = t;
    betasav = beta;
    KD2sav = KD2;
    save matrixKD2 KD2sav Hsav Nsav asav alphasav tsav betasav
else
    fclose(fid);
    load matrixKD2 KD2sav Hsav Nsav asav alphasav tsav betasav
    if ( (H == Hsav)  & (Nsav >= N) & (a == asav) & (alphasav == alpha) & (t == tsav) & (beta == betasav))
        KD2 = KD2sav(1:N+1,1:N+1);
    else
        KD2   = zeros(N+1,N+1);
        for k = 0:N
            for p = 0:N
                z = H*(zeta+1)/(2*a);
                val = alph*sum(-(1./permeability(z, beta)).*Dw(p,zeta).*Dv(k,zeta));
                val = val + alph*sum(-(1./permeability(z, beta)).*(zeta./(1-zeta.^2)).*Dw(p,zeta).*v(k,zeta));
                if abs(val) < cutoff
                    KD2(k+1,p+1) = 0;
                else
                    KD2(k+1,p+1) = (4/H^2)*val;
                end
            end
        end
        Hsav = H;
        Nsav = N;
        asav = a;
        alphasav = alpha;
        tsav = t;
        betasav = beta;
        KD2sav = KD2;
        save matrixKD2 KD2sav Hsav Nsav asav alphasav tsav betasav
    end
end

%%%%% matrix <(4/H^2)*D2W_p, V_k> %%%%%

AD2 = zeros(N+1,N+1);
for k = 0:N
    for p = 0:N
        
        if k == p
            AD2(k+1,p+1) = -p*( p + 3 ) - 2*c(p);
            
        elseif k == p-2
            AD2(k+1,p+1) = p*( p - 3 ) + 2*c(p);
            
        else
            AD2(k+1,p+1) = 0;
        end
        
    end
end

AD2 = (4/H^2)*(f/2)*AD2;


%%%%% matrix <(2/H)*DW_p, V_k> %%%%%
AD1 = zeros(N+1,N+1);
for k = 0:N
    for p = 0:N
        
        if k == p+1
            AD1(k+1,p+1) = ( -2*c(p) - p ) / 2;
            
        elseif k == p-1
            AD1(k+1,p+1) = ( 5*p + 2 - c(p-1)*( p + 2 ) - c(p-2)*( d(p-3) + 1 )*p ) / 2;
            
        elseif k == p-3
            AD1(k+1,p+1) = ( 2 - p ) / 2;
            
        else
            AD1(k+1,p+1) = 0;
        end
        
    end
end

AD1 = (2/H)*(f/2)*AD1;
