function [ dev_ineq, dev, J_ineq, J   ] = deviation_steady_state(x) 
 % [dev_ineq, dev, J_ineq, J] % use this for Jacobian with fmincon
 % [dev]  % use this if no Jacobian in use
 % [dev, J] % use if only solve model w/o objective, i.e. constraints only
 
% This function computes the deviation 'dev' of the equilibrium equations given 
% a vector x of endogenous variables, and provides the Jacobian 'J'.  
 
global L sig zeta nu kappa k rho Gamma Beta theta
global N S
global totvarnumb tradebalanceratios
global T_fixed g_fixed GDPlevels % added on 3/8/20
global nu_tilde % added on 18/10/21
global alpha % added on 13/11/21
global WAGEfrom WAGEto Yfrom Yto Pfrom Pto  RDlaborfrom RDlaborto THOLDfrom THOLDto 
global ETAfrom ETAto Tfrom Tto TAUfrom TAUto DELTAfrom DELTAto Ffrom Fto PsiCDfrom PsiCDto;  
global display_constr_value_prob; % added Jan 20, 2022
global solvector % added Jan 2022
global DELTA % added Jan 2022 for counterfactual exercises
global ETA_per ETA T_per T TAU
global TBratio
global f_origin f_dest g_0_exo
global obj_is_scalar

% ----------------------------------------------------------------------- %
%   could add inequality constraints here (fmincon syntax requires these) %
% ----------------------------------------------------------------------- %
dev_ineq = []; J_ineq = []; % don't need inequality constraints
     
% ----------------------------------------------------------------------- %
%   extract variables, and reshape to get correct format                  %
% ----------------------------------------------------------------------- %

x=real(x); 
P_input = [1; x(Pfrom:Pto)];
% P_input = ones(N,1); % to shut off for J check
WAGE = x(WAGEfrom:WAGEto);
% WAGE = ones(N,1); % to shut off for J check
Y = x(Yfrom:Yto);
% Y = ones(N,1); % to shut off for J check
Y = Y./P_input; % Simon suggests this works better
RDlabor_per = zeros(1,S+1,N);
RDlabor_per(:,2:end,:) = reshape(x(RDlaborfrom:RDlaborto), 1, S, N);  
% RDlabor_per = 0.*RDlabor_per + 1; % ad hoc to check J
RDlabor = permute(RDlabor_per, [3 2 1]);
THOLD=ones(N,S+1,N);
THOLD(:,2:end,:) = reshape(x(THOLDfrom:THOLDto), N, S, N); 
% THOLD = 0.*THOLD+2; % ad hoc to check J

% ----------------------------------------------------------------------- %
%   From here: compute variables                                          %
% ----------------------------------------------------------------------- %
gPsi_s = sum(ETA_per.*k./(k-1).*RDlabor_per.^(1-kappa),3)-zeta;
gPsi_s(1) = g_0_exo; % replace s=0 growth rate by exogenous value
g = sum(Beta./(sig-1).*gPsi_s,2)/sum(Beta.*alpha,2);

r = rho + g./Gamma;  

Phi = permute(T, [3 2 1]).*(TAU.*permute(WAGE,[3 2 1]).^alpha.*permute(P_input,[3 2 1]).^(1-alpha)).^-theta;
% ----------------------------------------------------------------------- %
%   aggregate qualities                                                   %
% ----------------------------------------------------------------------- %
PsiM_P_ND = (ETA_per.*RDlabor_per.^(1-kappa).*k./(k-1).* ...
    THOLD.^(1-k))./(gPsi_s+DELTA+nu+zeta);  

PsiM_NP = (ETA_per.*RDlabor_per.^(1-kappa).*k./(k-1).* ...
    (1-THOLD.^(1-k))+DELTA.*PsiM_P_ND)./(gPsi_s+nu+zeta); 

PsiM_PL = nu.*PsiM_P_ND ./ (gPsi_s+DELTA+nu_tilde+zeta);

PsiM_PD = nu_tilde.*PsiM_PL ./ (gPsi_s+DELTA+zeta);

PsiCL = (nu.*PsiM_NP+DELTA.*PsiM_PL) ./ (gPsi_s+nu_tilde+zeta); 

PsiM = PsiM_P_ND+PsiM_NP+PsiM_PL+PsiM_PD;     

PsiCD = 1-sum(PsiM+PsiCL,3);

% ----------------------------------------------------------------------- %
%   other implied variables                                               %
% ----------------------------------------------------------------------- %

numerator_PMoverP =     (sig./(sig-1)).^(1-sig) .* sum(PsiM.*Phi.^((sig-1)./theta),3);
numerator_PCLoverP = sum(PsiCL.*Phi.^((sig-1)./theta),3);
numerator_PCDoverP = PsiCD.*sum(Phi,3).^((sig-1)./theta);
denominator_allPoverP = numerator_PMoverP + numerator_PCLoverP + numerator_PCDoverP;

PMoverP = (numerator_PMoverP./denominator_allPoverP).^(1./(1-sig));

PCLoverP = (numerator_PCLoverP./denominator_allPoverP).^(1./(1-sig)); % denominator is the same as above.

PCDoverP = (numerator_PCDoverP./denominator_allPoverP).^(1./(1-sig)); % denominator is the same as above.
% go on here..still seems to be mistake in pcd over p. DIG DEEPER..check partials!!!!!!.
P = prod((gamma((theta+1-sig)./theta).*denominator_allPoverP).^(Beta./(1-sig)),2);

XM = PsiM.*Phi.^((sig-1)./theta) ./ sum(PsiM.*Phi.^((sig-1)./theta),3) .* PMoverP.^(1-sig) .* Beta.*P.*Y;
XM(isnan(XM)) = 0; % because of sector 0 issues with Inf and 0  

XCL = PsiCL.*Phi.^((sig-1)./theta) ./ sum(PsiCL.*Phi.^((sig-1)./theta),3) .* PCLoverP.^(1-sig) .* Beta.*P.*Y;
XCL(isnan(XCL)) = 0; % because of sector 0 issues with Inf and 0  

XCD = Phi./sum(Phi,3) .* PCDoverP.^(1-sig) .* Beta.*P.*Y;

% XCL = 0.*XCL;  % AD HOC, REMOVE
% XCL(:,2:end,:) = 0.04;% AD HOC, REMOVE
% XCD = 0.*XCD + 0.04;% AD HOC, REMOVE

% evening of sept 23. now i know that issue has to be with xcl and xcd..coz if these are set exogenously, then J is correct
% now the tricky part: I checked all individual derivatives for XCL, still didnt find the issue. But what I know is e.g. from changing RDlabor that
% "own sector" (for sectors 1 and 2, i.e. patenting sectors) gives problem, while the indirect one is fine. So, most likely seems a mistake on the
% "Main diagonal"... howevver, what is weird is that all individual derivatives seem ok.

Profnorm = XM ./ (sig.*PsiM.*permute(WAGE,[3 2 1]));
Profnorm(isnan(Profnorm)) = 0; % because of sector 0 issues with Inf and 0  



% ----------------------------------------------------------------------- %   
% ----------------------------------------------------------------------- %
%     stady state equations                                               %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %  
 
dev_sales_income = 1./permute(WAGE.*L, [3 2 1]) .* (sum(alpha.*sum(XCL+XCD+(1-1./sig).*XM,1),2) ...
    + permute(WAGE, [3 2 1]).*sum(RDlabor_per+sum(f_dest.*ETA.*RDlabor.^(1-kappa).*permute(THOLD, [3 2 1]).^-k + f_origin.*ETA_per.*RDlabor_per.^(1-kappa).*THOLD.^-k,1),2)) ...
    - ones(1,1,N);

TB_numerator = (permute(sum(sum(XCL+XCD+XM,1),2), [3 2 1]) + sum(sum(WAGE.*f_dest.*ETA_per(:,2:end,:).*RDlabor_per(:,2:end,:).^(1-kappa).*THOLD(:,2:end,:).^-k,3),2) ...
    - permute(sum(sum(WAGE.*f_dest.*ETA_per(:,2:end,:).*RDlabor_per(:,2:end,:).^(1-kappa).*THOLD(:,2:end,:).^-k,1),2), [3 2 1]) ...
    - (P.*Y)); % this will be used again for the derivatives
dev_trade_balance = TB_numerator ./ (sum(sum(sum(XCL+XCD+XM,3),2),1))  -  TBratio;

dev_RD = 1./RDlabor_per.^kappa .* ETA_per./(k-1).*sum(Profnorm.* ...
    k./(r+zeta+nu-g+gPsi_s) + (f_origin+f_dest.*WAGE./permute(WAGE,[3 2 1])) .* THOLD.^-k, 1) - ones(1,S+1,N);
dev_RD = dev_RD(:,2:end,:);

dev_PatThold = ( (f_origin+f_dest.*WAGE./permute(WAGE,[3 2 1])) .* (r+zeta(2:end)+DELTA(:,2:end)-g+gPsi_s(2:end)) ...
    .* (r+zeta(2:end)+nu(2:end)-g+gPsi_s(2:end)+DELTA(:,2:end)) ./ (Profnorm(:,2:end,:).*nu(2:end)) ) ./ THOLD(:,2:end,:) - ones(N,S,N);

dev_PriceRatio = P(2:end)./P(1) - P_input(2:end);

dev_PriceNorm = P(1) - 1;

%% compute deviation
dev =  [1.*dev_sales_income(:); ...
       1.*dev_trade_balance(:); ...
       1.*dev_RD(:); ...
       1.*dev_PatThold(:); ...
       1.*dev_PriceRatio; ...
       1.*dev_PriceNorm];


display_probability = 0.1;
if rand<display_probability 
    norm(dev)  % display current norm of deviation
end


if nargout>1

%%
% ----------------------------------------------------------------------- %   
% ----------------------------------------------------------------------- %
%     derivatives                                                         %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %                   

% ----------------------------------------------------------------------- %
%   g_Psi_s - w.r.t. RDlabor, (plus ETA for calibration)                  %
% ----------------------------------------------------------------------- %
Dg_Psi_s=zeros(S, totvarnumb);  
Dg_Psi_s(:, RDlaborfrom:RDlaborto) = reshape(ETA_per(:,2:end,:).*k./(k-1) ...
    .*(1-kappa).*RDlabor_per(:,2:end,:).^(-kappa).*eye(S), S,[]);  
% Dg_Psi_s(:, ETAfrom:ETAto) = reshape(k./(k-1).*RDlabor_per(:,2:end,:).^(1-kappa).*eye(S), S,[]);  % for later (calibration)


% ----------------------------------------------------------------------- %
%   g            %
% ----------------------------------------------------------------------- %
Dg__g_Psi_s = Beta(:,2:end)./(sig(:,2:end)-1)./(sum(Beta.*alpha,2));
Dg = sum(Dg__g_Psi_s.*permute(Dg_Psi_s, [2 1]),2);


% ----------------------------------------------------------------------- %
%   Phi_nsi derivatives - w.r.t. P_input, WAGE, (TAU, T)                  %
% ----------------------------------------------------------------------- %
DPhi_nsi=zeros(N,S+1,N, totvarnumb);  

DPhi_nsi(:,:,:, Pfrom:Pto) = -(1-alpha).*theta.*Phi./ ...
    permute(P_input(2:end), [2 3 4 1]).*permute([zeros(1,N-1); eye(N-1)], [3 4 1 2]); 

%  abcd = -(1-alpha).*theta.*Phi./ ...
%     permute(P_input(1:end), [2 3 4 1]).*permute([eye(N)], [3 4 1 2])
% abcd = abcd(:,:,:,2:end); %  alternative, but the above seems indeed to be the same. can eventually delete this once J correct.

DPhi_nsi(:,:,:, WAGEfrom:WAGEto) = -alpha.*theta.*Phi./ ...
    permute(WAGE, [2 3 4 1]).*permute(eye(N), [3 4 1 2]); 



% ----------------------------------------------------------------------- %
% Psi_P_ND derivatives - w.r.t. RDlabor (order: 1,s,i!), THOLD, g Psi, (ETA, DELTA) %
% ----------------------------------------------------------------------- %
DPsiM_P_ND=zeros(N,S,N, totvarnumb); 

DPsiM_P_ND(:,:,:, RDlaborfrom:RDlaborto) = reshape((1-kappa).* ...
    ETA_per(:,2:end,:).*permute(RDlabor_per(:,2:end,:),[4 5 6 1 2 3]).^-kappa.*k./(k-1).* ...
    THOLD(:,2:end,:).^(1-k)./(gPsi_s(2:end)+DELTA(:,2:end)+nu(2:end)+zeta(2:end)).* permute(eye(N) ...
    , [3 4 1 5 6 2]).*permute(eye(S), [3 1 4 5 2]), N,S,N,[]);  

DPsiM_P_ND(:,:,:, THOLDfrom:THOLDto) = reshape(ETA_per(:,2:end,:).* ...
    RDlabor_per(:,2:end,:).^(1-kappa).*-k.*permute(THOLD(:,2:end,:), [4 5 6 1 2 3]).^(-k)./ ...
    (gPsi_s(2:end)+DELTA(:,2:end)+nu(2:end)+zeta(2:end)).*permute(eye(N),[3 4 1 5 6 2]).* ...
    permute(eye(S),[3 1 4 5 2 ]).*permute(eye(N),[1 3 4 2]), N,S,N,[]);

DPsiM_P_ND__g_Psi_s=zeros(N,S,N, S);
DPsiM_P_ND__g_Psi_s = -PsiM_P_ND(:,2:end,:)./(gPsi_s(2:end)+DELTA(:,2:end)+nu(2:end)+zeta(2:end)) ...
    .*permute(eye(S),[3 1 4 2]);



DPsiM_P_ND = DPsiM_P_ND + sum(permute(DPsiM_P_ND__g_Psi_s, [1 2 3 5 4]).*permute(Dg_Psi_s, [3 4 5 2 1]), 5);



% ----------------------------------------------------------------------- %
%   PsiM_NP derivatives - w.r.t. RDlabor (1si order!), THOLD, g Psi, (ETA, DELTA) %
% ----------------------------------------------------------------------- %
DPsiM_NP=zeros(N,S,N, totvarnumb);


DPsiM_NP(:,:,:, RDlaborfrom:RDlaborto) = reshape( ((1-kappa).*ETA_per(:,2:end,:).* ...
    RDlabor_per(:,2:end,:).^-kappa.*k./(k-1).*(1-THOLD(:,2:end,:).^(1-k))) ./ ...
    (gPsi_s(2:end)+nu(2:end)+zeta(2:end)).*permute(eye(N), [3 4 1 5 6 2]).* ...
    permute(eye(S), [3 1 4 5 2]), N,S,N,[]);  % i think this is correct


DPsiM_NP(:,:,:, THOLDfrom:THOLDto) = reshape(ETA_per(:,2:end,:).* ...
    RDlabor_per(:,2:end,:).^(1-kappa).*k.*permute(THOLD(:,2:end,:), [4 5 6 1 2 3]).^(-k) ...
    ./(gPsi_s(2:end)+nu(2:end)+zeta(2:end)).*permute(eye(N), [3 4 1 5 6 2]).* ...
    permute(eye(S),[3 1 4 5 2]).*permute(eye(N),[1 3 4 2]),N,S,N,[]);

DPsiM_NP__g_Psi_s=zeros(N,S,N, S);
DPsiM_NP__g_Psi_s = -PsiM_NP(:,2:end,:)./(gPsi_s(2:end)+nu(2:end)+zeta(2:end)) ...
    .*permute(eye(S),[3 1 4 2]);

  %  .*permute(eye(S),[3 1 4 2]); % i think this is not needed
  

DPsiM_NP__PsiM_P_ND=zeros(N,S,N, N,S,N);
DPsiM_NP__PsiM_P_ND = DELTA(:,2:end)./(gPsi_s(:,2:end)+nu(:,2:end)+zeta(:,2:end)).*permute(eye(N), [3 4 1 5 6 2]).* ...
    permute(eye(S), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);

DPsiM_NP = DPsiM_NP + sum(permute(DPsiM_NP__g_Psi_s, [1 2 3 5 4]).*permute(Dg_Psi_s, [3 4 5 2 1]), 5) ...
    + sum(sum(sum(permute(DPsiM_NP__PsiM_P_ND, [1 2 3 7 4 5 6]).*permute(DPsiM_P_ND, [5 6 7 4 1 2 3]), 7), 6), 5);



% ----------------------------------------------------------------------- %
% NEW: PsiM_P_L derivatives - w.r.t. RDlabor (order: 1,s,i!), THOLD, g Psi, (ETA, DELTA) %
% ----------------------------------------------------------------------- %
DPsiM_P_L=zeros(N,S,N, totvarnumb); 

DPsiM_P_L__g_Psi_s=zeros(N,S,N, S);
DPsiM_P_L__g_Psi_s = -nu(2:end).*PsiM_P_ND(:,2:end,:)./(gPsi_s(2:end)+DELTA(:,2:end)+nu_tilde(2:end)+zeta(2:end)).^2 ...
    .*permute(eye(S),[3 1 4 2]);

DPsiM_P_L__PsiM_P_ND=zeros(N,S,N, N,S,N);
DPsiM_P_L__PsiM_P_ND = nu(:,2:end)./(gPsi_s(:,2:end)+DELTA(:,2:end)+nu_tilde(2:end)+zeta(:,2:end)).*permute(eye(N), [3 4 1 5 6 2]).* ...
    permute(eye(S), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);

DPsiM_P_L = DPsiM_P_L + sum(permute(DPsiM_P_L__g_Psi_s, [1 2 3 5 4]).*permute(Dg_Psi_s, [3 4 5 2 1]), 5) ...
    + sum(sum(sum(permute(DPsiM_P_L__PsiM_P_ND, [1 2 3 7 4 5 6]).*permute(DPsiM_P_ND, [5 6 7 4 1 2 3]), 7), 6), 5);




% ----------------------------------------------------------------------- %
%   PsiM_P_D derivatives - w.r.t. RDlabor (1,s,i order), THOLD, g Psi, (ETA)                %
% ----------------------------------------------------------------------- %
DPsiM_PD=zeros(N,S,N, totvarnumb);

DPsiM_P_D__g_Psi_s = -nu_tilde(2:end).*PsiM_PL(:,2:end,:)./ ...
    (gPsi_s(2:end)+DELTA(:,2:end)+zeta(2:end)).^2 ...
    .*permute(eye(S),[3 1 4 2]);

DPsiM_P_D__PsiM_P_L=zeros(N,S,N, N,S,N);
DPsiM_P_D__PsiM_P_L = nu_tilde(:,2:end)./(gPsi_s(:,2:end)+DELTA(:,2:end)+zeta(:,2:end)).*permute(eye(N), [3 4 1 5 6 2]).* ...
    permute(eye(S), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);

DPsiM_PD = DPsiM_PD + sum(permute(DPsiM_P_D__g_Psi_s, [1 2 3 5 4]).*permute(Dg_Psi_s, [3 4 5 2 1]), 5) ...
    + sum(sum(sum(permute(DPsiM_P_D__PsiM_P_L, [1 2 3 7 4 5 6]).*permute(DPsiM_P_L, [5 6 7 4 1 2 3]), 7), 6), 5);



% ----------------------------------------------------------------------- %
%  PsiC_L derivatives - w.r.t. RDlabor (order: 1,s,i!), THOLD, g psi, (ETA, DELTA) %
% ----------------------------------------------------------------------- %
DPsiC_L=zeros(N,S,N, totvarnumb); 

DPsiC_L__g_Psi_s = -1./(gPsi_s(2:end)+nu_tilde(2:end)+zeta(2:end)).^2 ...
    .* (nu(2:end).*PsiM_NP(:,(2:end),:) + DELTA(:,2:end).*PsiM_PL(:,(2:end),:)) ...
    .* permute(eye(S), [3 1 4 2]);

DPsiC_L__PsiM_NP=zeros(N,S,N, N,S,N);
DPsiC_L__PsiM_NP = nu(:,2:end)./(gPsi_s(:,2:end)+nu_tilde(:,2:end)+zeta(:,2:end)).*permute(eye(N), [3 4 1 5 6 2]).* ...
    permute(eye(S), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);

DPsiC_L__PsiM_P_L=zeros(N,S,N, N,S,N);
DPsiC_L__PsiM_P_L = DELTA(:,2:end)./(gPsi_s(:,2:end)+DELTA(:,2:end)+zeta(:,2:end)).*permute(eye(N), [3 4 1 5 6 2]).* ...
    permute(eye(S), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);

DPsiC_L = DPsiC_L + sum(permute(DPsiC_L__g_Psi_s, [1 2 3 5 4]).*permute(Dg_Psi_s, [3 4 5 2 1]), 5) ...
        + sum(sum(sum(permute(DPsiC_L__PsiM_NP, [1 2 3 7 4 5 6]).*permute(DPsiM_NP, [5 6 7 4 1 2 3]), 7), 6), 5) ...
            + sum(sum(sum(permute(DPsiC_L__PsiM_P_L, [1 2 3 7 4 5 6]).*permute(DPsiM_P_L, [5 6 7 4 1 2 3]), 7), 6), 5);




% ----------------------------------------------------------------------- %
%   PsiM (aggregate) derivatives  - w.r.t. all variables                  %
% ----------------------------------------------------------------------- %
DPsiM = zeros(N,S,N, totvarnumb);

DPsiM = DPsiM_P_ND + DPsiM_NP + DPsiM_P_L + DPsiM_PD; 




% ----------------------------------------------------------------------- %
%   PM derivatives  - w.r.t. (symbol:"__") Phi_nsi, PsiM, PsiCL           %
% ----------------------------------------------------------------------- %
Ddenominator__Phi_nsi = ((sig./(sig-1)).^(1-sig).*permute(PsiM, [4 5 6 1 2 3]).*(sig-1)./theta ... % this expression will be used a few times
    .*permute(Phi.^((sig-1)./theta-1),[4 5 6 1 2 3]) ...
    +permute(PsiCL, [4 5 6 1 2 3]).*(sig-1)./theta.*permute(Phi.^((sig-1)./theta-1), [4 5 6 1 2 3]) ...
    +permute(PsiCD, [4 5 6 1 2 3]).*(sig-1)./theta.*(sum(Phi,3)).^((sig-1)./theta-1)) ...
    .*permute(eye(S+1), [3 1 4 5 2]) .* permute(eye(N), [1 3 4 2]); % last line added on sept 21, 
% I think its necessary (or not? if always impose dimension-consistency later..but then, still cannot hurt?). mistakes are currently here in these derivatives.

DPMoverP__Phi_nsi = zeros(N,S+1,1, N,S+1,N); 
DPMoverP__Phi_nsi = 1./(1-sig).*PMoverP.^sig ...
    .*( (sig/(sig-1)).^(1-sig).*permute(PsiM, [4 5 6 1 2 3]).*(sig-1)./theta ...
    .*permute(Phi.^((sig-1)./theta-1),[4 5 6 1 2 3])./denominator_allPoverP ... 
    - numerator_PMoverP.*Ddenominator__Phi_nsi./denominator_allPoverP.^2) ...
    .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);
DPMoverP__Phi_nsi(isnan(DPMoverP__Phi_nsi)) = 0; % because of sector 0 issues with Inf and 0  


DPMoverP__PsiM = zeros(N,S+1,1, N,S+1,N); 
DPMoverP__PsiM = 1./(1-sig).*PMoverP.^sig ...
    .*( (sig/(sig-1)).^(1-sig).*permute(Phi.^((sig-1)./theta),[4 5 6 1 2 3])./denominator_allPoverP ...
    - (sig/(sig-1)).^(1-sig).*permute(sum(PsiM ...
    .*Phi.^((sig-1)./theta), 3), [4 5 6 1 2 3]).*((sig/(sig-1)).^(1-sig).*permute(Phi.^((sig-1)./theta),[4 5 6 1 2 3]) ...
    - permute(sum(Phi, 3).^((sig-1)./theta), [4 5 6 1 2 3]) ) ./ denominator_allPoverP.^2 ) ...
    .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);
DPMoverP__PsiM(isnan(DPMoverP__PsiM)) = 0; % because of sector 0 issues with Inf and 0  


DPMoverP__PsiCL = zeros(N,S+1,1, N,S+1,N); 
DPMoverP__PsiCL = 1./(1-sig).*PMoverP.^sig ...
    .* -(sig/(sig-1)).^(1-sig).*permute(sum(PsiM.*Phi.^((sig-1)./theta),3), [4 5 6 1 2 3]) ...
    .*(permute(Phi.^((sig-1)./theta)-sum(Phi,3).^((sig-1)./theta), [4 5 6 1 2 3])) ...
    ./ denominator_allPoverP.^2  ...
    .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);
DPMoverP__PsiCL(isnan(DPMoverP__PsiCL)) = 0; % because of sector 0 issues with Inf and 0  


DPMoverP = sum(sum(sum(permute(DPMoverP__Phi_nsi, [1 2 3 7 4 5 6]).*permute(DPhi_nsi, [5 6 7 4 1 2 3]), 7), 6), 5);
DPMoverP(:,2:end,:,:) = DPMoverP(:,2:end,:,:) ...
    + sum(sum(sum(permute(DPMoverP__PsiM(:,2:end,:,:,2:end,:), [1 2 3 7 4 5 6]).*permute(DPsiM,  [5 6 7 4 1 2 3]), 7), 6), 5) ...
    + sum(sum(sum(permute(DPMoverP__PsiCL(:,2:end,:,:,2:end,:),[1 2 3 7 4 5 6]).*permute(DPsiC_L,[5 6 7 4 1 2 3]), 7), 6), 5);


% ----------------------------------------------------------------------- %
%   PCL derivatives  - w.r.t. (symbol:"__") Phi_nsi, PsiM, PsiCL           %
% ----------------------------------------------------------------------- %
DPCLoverP__Phi_nsi = zeros(N,S+1,1, N,S+1,N); 
DPCLoverP__Phi_nsi = 1./(1-sig).*PCLoverP.^sig ...
    .*( permute(PsiCL, [4 5 6 1 2 3]).*(sig-1)./theta ...
    .*permute(Phi.^((sig-1)./theta-1),[4 5 6 1 2 3])./denominator_allPoverP ... 
    - numerator_PCLoverP.*Ddenominator__Phi_nsi./denominator_allPoverP.^2) ...
    .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);
DPCLoverP__Phi_nsi(isnan(DPCLoverP__Phi_nsi)) = 0; % because of sector 0 issues with Inf and 0  


DPCLoverP__PsiCL = zeros(N,S+1,1, N,S+1,N); 
DPCLoverP__PsiCL = 1./(1-sig).*PCLoverP.^sig ...
    .*( permute(Phi.^((sig-1)./theta),[4 5 6 1 2 3])./denominator_allPoverP ...
    - permute(sum(PsiCL ...
    .*Phi.^((sig-1)./theta), 3), [4 5 6 1 2 3]).*(permute(Phi.^((sig-1)./theta),[4 5 6 1 2 3]) ...
    - permute(sum(Phi, 3).^((sig-1)./theta), [4 5 6 1 2 3]) ) ./ denominator_allPoverP.^2 ) ...
    .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);
DPCLoverP__PsiCL(isnan(DPCLoverP__PsiCL)) = 0; % because of sector 0 issues with Inf and 0  


DPCLoverP__PsiM = zeros(N,S+1,1, N,S+1,N); 
DPCLoverP__PsiM = 1./(1-sig).*PCLoverP.^sig ...
    .* -permute(sum(PsiCL.*Phi.^((sig-1)./theta),3), [4 5 6 1 2 3]) ...
    .*((sig./(sig-1)).^(1-sig).*permute(Phi.^((sig-1)./theta), [4 5 6 1 2 3])-permute(sum(Phi,3).^((sig-1)./theta), [4 5 6 1 2 3])) ... % sept 24: finally found a small mistake here!!
    ./ denominator_allPoverP.^2  ...
    .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);
DPCLoverP__PsiM(isnan(DPCLoverP__PsiM)) = 0; % because of sector 0 issues with Inf and 0  


DPCLoverP = sum(sum(sum(permute(DPCLoverP__Phi_nsi, [1 2 3 7 4 5 6]).*permute(DPhi_nsi, [5 6 7 4 1 2 3]), 7), 6), 5);
DPCLoverP(:,2:end,:,:) = DPCLoverP(:,2:end,:,:) ...
    + sum(sum(sum(permute(DPCLoverP__PsiM(:,2:end,:,:,2:end,:), [1 2 3 7 4 5 6]).*permute(DPsiM,  [5 6 7 4 1 2 3]), 7), 6), 5) ...
    + sum(sum(sum(permute(DPCLoverP__PsiCL(:,2:end,:,:,2:end,:),[1 2 3 7 4 5 6]).*permute(DPsiC_L,[5 6 7 4 1 2 3]), 7), 6), 5);



% ----------------------------------------------------------------------- %
%   PCD derivatives  - w.r.t. (symbol:"__") Phi_nsi, PsiM, PsiCL           %
% ----------------------------------------------------------------------- %
DPCDoverP__Phi_nsi = zeros(N,S+1,1, N,S+1,N); 
DPCDoverP__Phi_nsi = 1./(1-sig).*PCDoverP.^sig ...
    .*( permute(PsiCD.*(sig-1)./theta.*sum(Phi,3).^((sig-1)./theta-1),[4 5 6 1 2 3])./denominator_allPoverP ... 
    - numerator_PCDoverP.*Ddenominator__Phi_nsi./denominator_allPoverP.^2) ...
    .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);
DPCDoverP__Phi_nsi(isnan(DPCDoverP__Phi_nsi)) = 0; % because of sector 0 issues with Inf and 0  


DPCDoverP__PsiCL = zeros(N,S+1,1, N,S+1,N); 
DPCDoverP__PsiCL = 1./(1-sig).*PCDoverP.^sig ...
    .*( -permute(sum(Phi,3).^((sig-1)./theta),[4 5 6 1 2 3])./denominator_allPoverP ... % here was mistake with sum. now ok it seems
    - (permute(PsiCD .* sum(Phi,3).^((sig-1)./theta), [4 5 6 1 2 3])...
    .*(permute(Phi.^((sig-1)./theta)-sum(Phi,3).^((sig-1)./theta),[4 5 6 1 2 3]))) ...
    ./ denominator_allPoverP.^2 ) ...
    .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);
DPCDoverP__PsiCL(isnan(DPCDoverP__PsiCL)) = 0; % because of sector 0 issues with Inf and 0  


DPCDoverP__PsiM = zeros(N,S+1,1, N,S+1,N); 
DPCDoverP__PsiM = 1./(1-sig).*PCDoverP.^sig ... % exploit symmetry (almost!! thers is small difference in derivative of denominator)  
    .*( -permute(sum(Phi,3).^((sig-1)./theta),[4 5 6 1 2 3])./denominator_allPoverP ...   % here was mistake with sum. now ok it seems
    - (permute(PsiCD .* sum(Phi,3).^((sig-1)./theta), [4 5 6 1 2 3])...  
    .*(permute((sig./(sig-1)).^(1-sig).*Phi.^((sig-1)./theta)-sum(Phi,3).^((sig-1)./theta),[4 5 6 1 2 3]))) ... % here is the only small difference: (sig./(sig-1)).^(1-sig) 
    ./ denominator_allPoverP.^2 ) ...
    .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);
DPCDoverP__PsiM(isnan(DPCDoverP__PsiM)) = 0; % because of sector 0 issues with Inf and 0  


DPCDoverP = sum(sum(sum(permute(DPCDoverP__Phi_nsi, [1 2 3 7 4 5 6]).*permute(DPhi_nsi, [5 6 7 4 1 2 3]), 7), 6), 5);
DPCDoverP(:,2:end,:,:) = DPCDoverP(:,2:end,:,:) ...
    + sum(sum(sum(permute(DPCDoverP__PsiM(:,2:end,:,:,2:end,:), [1 2 3 7 4 5 6]).*permute(DPsiM,  [5 6 7 4 1 2 3]), 7), 6), 5) ...
    + sum(sum(sum(permute(DPCDoverP__PsiCL(:,2:end,:,:,2:end,:),[1 2 3 7 4 5 6]).*permute(DPsiC_L,[5 6 7 4 1 2 3]), 7), 6), 5);


DPCDoverP = 1./(1-sig).*numerator_PCDoverP.^(sig./(1-sig)) ...
    .* (-(1-sig).*numerator_PMoverP.^(-sig./(1-sig)).*DPMoverP-(1-sig).*numerator_PCLoverP.^(-sig./(1-sig)).*DPCLoverP); % different approach -- sept 24 ff

% ----------------------------------------------------------------------- %
%   P_n derivatives  - w.r.t. (symbol:"__") Phi_nsi, PsiM, PsiCL          %
% ----------------------------------------------------------------------- %
    BracketTerm = (sig./(sig-1)).^(1-sig).*sum(PsiM.*Phi.^((sig-1)./theta),3)+sum(PsiCL.*Phi.^((sig-1)./theta),3)...
        + PsiCD.*sum(Phi,3).^((sig-1)./theta);


    DP__Phi = P.*permute(Beta./(1-sig), [1 3 4 5 2]).*1./permute(BracketTerm, [3 4 5 1 2])...
        .* permute( ((sig./(sig-1)).^(1-sig).*PsiM.*(sig-1)./theta.*Phi.^((sig-1)./theta-1) ...
        + PsiCL.*(sig-1)./theta.*Phi.^((sig-1)./theta-1) ...
        + PsiCD.*(sig-1)./theta.*sum(Phi,3).^((sig-1)./theta-1)), [4 5 6 1 2 3]) ...
        .* permute(eye(N), [1 3 4 2]);


    DP__PsiM = P.*permute(Beta./(1-sig), [1 3 4 5 2]).*1./permute(BracketTerm, [3 4 5 1 2])...
        .* permute( (sig./(sig-1)).^(1-sig).*Phi.^((sig-1)./theta) ...
        - sum(Phi,3).^((sig-1)./theta), [4 5 6 1 2 3]) ...
        .* permute(eye(N), [1 3 4 2]);


    DP__PsiCL = P.*permute(Beta./(1-sig), [1 3 4 5 2]).*1./permute(BracketTerm, [3 4 5 1 2])...
        .* permute( Phi.^((sig-1)./theta) ...
        - sum(Phi,3).^((sig-1)./theta), [4 5 6 1 2 3]) ...
        .* permute(eye(N), [1 3 4 2]);


    DP = sum(sum(sum(...
        permute(DP__Phi, [1 2 3 7 4 5 6]).*permute(DPhi_nsi, [5 6 7 4 1 2 3]) ...
        ,7), 6), 5);
    DP = DP + sum(sum(sum(...
        permute(DP__PsiM(:,:,:,:,2:end,:), [1 2 3 7 4 5 6]).*permute(DPsiM, [5 6 7 4 1 2 3]) ...
        + permute(DP__PsiCL(:,:,:,:,2:end,:), [1 2 3 7 4 5 6]).*permute(DPsiC_L, [5 6 7 4 1 2 3]) ...
        ,7), 6), 5);


    % ----------------------------------------------------------------------- %
    %   XM (monop) trade flows derivatives                                    %
    % ----------------------------------------------------------------------- %
    DXM__Phi_nsi = (sig-1)./theta.*(XM./permute(Phi, [4 5 6 1 2 3]).*permute(eye(N), [3 4 1 5 6 2]) ...
        -XM./sum(PsiM.*Phi.^((sig-1)./theta),3) ...
        .* permute(PsiM .* Phi.^((sig-1)./theta-1), [4 5 6 1 2 3])) ...
        .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);
    DXM__Phi_nsi(isnan(DXM__Phi_nsi)) = 0; % because of sector 0 issues with Inf and 0  


    DXM__PsiM = (XM./permute(PsiM, [4 5 6 1 2 3]).*permute(eye(N), [3 4 1 5 6 2]) ...
        - XM./sum(PsiM.*Phi.^((sig-1)./theta), 3).*permute(Phi.^((sig-1)./theta), [4 5 6 1 2 3])) ...
        .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);
    DXM__PsiM(isnan(DXM__PsiM)) = 0; % because of sector 0 issues with Inf and 0  


    DXM__PMoverP = (1-sig).*XM./permute(PMoverP, [3 4 5 1 2]) ...
        .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);


    DXM__P = XM./permute(P, [2 3 4 1]).*permute(eye(N), [1 3 4 2]);


    DXM = zeros(N,S+1,N,totvarnumb);
    DXM(:,:,:,Yfrom:Yto) = XM./permute(Y, [2 3 4 1]).*permute(eye(N), [1 3 4 2]);

    DXM = DXM + sum(sum(sum(...
        permute(DXM__Phi_nsi, [1 2 3 7 4 5 6]).*permute(DPhi_nsi, [5 6 7 4 1 2 3]) ...
        , 7), 6), 5) ...
        + sum(sum(sum(...
        permute(DXM__PsiM(:,:,:,:,2:end,:), [1 2 3 7 4 5 6]).*permute(DPsiM, [5 6 7 4 1 2 3]) ...
        , 7), 6), 5) ...
        + sum(sum(permute(DXM__PMoverP, [1 2 3 6 4 5]).*permute(DPMoverP, [3 5 6 4 1 2]) ...
        , 6), 5) ...
        + sum(permute(DXM__P, [1 2 3 5 4]).*permute(DP, [2 3 5 4 1]), 5);



    % ----------------------------------------------------------------------- %
    %   XCL (locally diffused, comp) trade flows derivatives                  %
    % ----------------------------------------------------------------------- %
    DXCL__Phi_nsi = (sig-1)./theta.*(XCL./permute(Phi, [4 5 6 1 2 3]).*permute(eye(N), [3 4 1 5 6 2]) ...
        -XCL./sum(PsiCL.*Phi.^((sig-1)./theta),3) ...
        .* permute(PsiCL .* Phi.^((sig-1)./theta-1), [4 5 6 1 2 3])) ...
        .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]); % old approach. but below seems to be the same..weird, Im kind of convinced mistake is in DXCL... (DXM is fine).
    DXCL__Phi_nsi = (sig-1)./theta.*XCL./permute(Phi, [4 5 6 1 2 3]).*permute(eye(N), [3 4 1 5 6 2]) .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]) ...
        -XCL .* permute(1./sum(PsiCL.*Phi.^((sig-1)./theta),3).*PsiCL .*(sig-1)./theta.* Phi.^((sig-1)./theta-1), [4 5 6 1 2 3]) ...
        .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);  % alternative. but I think its just equivalent to above.. (so can pick either one)
    DXCL__Phi_nsi = (sig-1)./theta.*(XCL./permute(Phi, [4 5 6 1 2 3]).*permute(eye(N), [3 4 1 5 6 2]) ...
        -XCL./sum(PsiCL.*Phi.^((sig-1)./theta),3) ...
        .* permute(PsiCL .* Phi.^((sig-1)./theta-1), [4 5 6 1 2 3])) ...
        .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]); %finally, on sept 24, copy pasted from DXM coz it should be virtually the same... but, also that one is just identical to above
    DXCL__Phi_nsi(isnan(DXCL__Phi_nsi)) = 0; % because of sector 0 issues with Inf and 0  


    DXCL__PsiCL = (XCL./permute(PsiCL, [4 5 6 1 2 3]).*permute(eye(N), [3 4 1 5 6 2]) ...
        - XCL./sum(PsiCL.*Phi.^((sig-1)./theta), 3).*permute(Phi.^((sig-1)./theta), [4 5 6 1 2 3])) ...
        .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);
    DXCL__PsiCL = (XCL./permute(PsiCL, [4 5 6 1 2 3]).*permute(eye(N), [3 4 1 5 6 2]) ...
        - XCL./sum(PsiCL.*Phi.^((sig-1)./theta), 3).*permute(Phi.^((sig-1)./theta), [4 5 6 1 2 3])) ...
        .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]); %finally, on sept 24, copy pasted from DXM coz it should be virtually the same...
    DXCL__PsiCL(isnan(DXCL__PsiCL)) = 0; % because of sector 0 issues with Inf and 0  


    DXCL__PCLoverP = (1-sig).*XCL./permute(PCLoverP, [3 4 5 1 2]) ...
        .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);
    DXCL__PCLoverP = (1-sig).*XCL./permute(PCLoverP, [3 4 5 1 2]) ...
        .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);%finally, on sept 24, copy pasted from DXM coz it should be virtually the same...

    DXCL__P = XCL./permute(P, [2 3 4 1]).*permute(eye(N), [1 3 4 2]);
    DXCL__P = XCL./permute(P, [2 3 4 1]).*permute(eye(N), [1 3 4 2]);


    DXCL = zeros(N,S+1,N,totvarnumb);
    DXCL(:,:,:,Yfrom:Yto) = XCL./permute(Y, [2 3 4 1]).*permute(eye(N), [1 3 4 2]);

    DXCL = DXCL + sum(sum(sum(...
        permute(DXCL__Phi_nsi, [1 2 3 7 4 5 6]).*permute(DPhi_nsi, [5 6 7 4 1 2 3]) ...
        , 7), 6), 5) ...
        + sum(sum(sum(...
        permute(DXCL__PsiCL(:,:,:,:,2:end,:), [1 2 3 7 4 5 6]).*permute(DPsiC_L, [5 6 7 4 1 2 3]) ...
        , 7), 6), 5) ...
        + sum(sum(permute(DXCL__PCLoverP, [1 2 3 6 4 5]).*permute(DPCLoverP, [3 5 6 4 1 2]) ...
        , 6), 5) ...
        + sum(permute(DXCL__P, [1 2 3 5 4]).*permute(DP, [2 3 5 4 1]), 5);

%     DXCL = 0.*DXCL; % AD HOC, REMOVE!!!!!!!! XXXXXXX
    
    % ----------------------------------------------------------------------- %
    %   XCD (intl' diffused, comp) trade flows derivatives                    %
    % ----------------------------------------------------------------------- %
    DXCD__Phi_nsi = (XCD./permute(Phi, [4 5 6 1 2 3]).*permute(eye(N), [3 4 1 5 6 2]) ...
        -XCD./permute(sum(Phi,3), [4 5 6 1 2 3])) ...
        .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);
    DXCD__Phi_nsi(isnan(DXCD__Phi_nsi)) = 0; % because of sector 0 issues with Inf and 0  


    DXCD__PCDoverP = (1-sig).*XCD./permute(PCDoverP, [3 4 5 1 2]) ...
        .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]);


    DXCD__P = XCD./permute(P, [2 3 4 1]).*permute(eye(N), [1 3 4 2]);


    DXCD = zeros(N,S+1,N,totvarnumb);
    DXCD(:,:,:,Yfrom:Yto) = XCD./permute(Y, [2 3 4 1]).*permute(eye(N), [1 3 4 2]);

    DXCD = DXCD + sum(sum(sum(...
        permute(DXCD__Phi_nsi, [1 2 3 7 4 5 6]).*permute(DPhi_nsi, [5 6 7 4 1 2 3]) ...
        , 7), 6), 5) ...
        + sum(sum(permute(DXCD__PCDoverP, [1 2 3 6 4 5]).*permute(DPCDoverP, [3 5 6 4 1 2]) ...
        , 6), 5) ...
        + sum(permute(DXCD__P, [1 2 3 5 4]).*permute(DP, [2 3 5 4 1]), 5);
    
%         DXCD = 0.*DXCD; % AD HOC, REMOVE!!!!!!!! XXXXXXX

    % ----------------------------------------------------------------------- %
    %   Profnorm (normalized profits) derivatives                             %
    % ----------------------------------------------------------------------- %
    DProfnorm__XM = Profnorm./permute(XM, [4 5 6 1 2 3]) ...
        .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]).*permute(eye(N), [3 4 1 5 6 2]);
    DProfnorm__XM(isnan(DProfnorm__XM)) = 0; % because of sector 0 issues with Inf and 0  


    DProfnorm__XPsiM = -Profnorm./permute(PsiM, [4 5 6 1 2 3]) ...
        .*permute(eye(S+1), [3 1 4 5 2]).*permute(eye(N), [1 3 4 2]).*permute(eye(N), [3 4 1 5 6 2]);
    DProfnorm__XPsiM(isnan(DProfnorm__XPsiM)) = 0; % because of sector 0 issues with Inf and 0  
    
    
    DProfnorm = zeros(N,S+1,N,totvarnumb);    
    DProfnorm(:,:,:,WAGEfrom:WAGEto) = -Profnorm./permute(WAGE, [2 3 4 1]) ...
        .*permute(eye(N), [3 4 1 2]);    
    
    DProfnorm = DProfnorm + sum(sum(sum(...
        permute(DProfnorm__XM, [1 2 3 7 4 5 6]).*permute(DXM, [5 6 7 4 1 2 3]) ...
        , 7), 6), 5) ...
        + sum(sum(sum(...
        permute(DProfnorm__XPsiM(:,:,:,:,2:end,:), [1 2 3 7 4 5 6]).*permute(DPsiM, [5 6 7 4 1 2 3]) ...
        , 7), 6), 5);
    
    

%% Equilibrium equations derivatives                                                

% ----------------------------------------------------------------------- %
%   1. PRODUCTION=SALES, (1,1,N),  w.r.t. XM, XC, WAGE           %
% ----------------------------------------------------------------------- %
Dsales_income__XM = zeros(1,1,N, N,S,N);
Dsales_income__XM = 1./permute(WAGE.*L, [3 2 1]) .* ones(1,1,N).*permute(ones(N,1,N).*(1-1./sig),[4 5 6 1 2 3])...
    .*permute(eye(N), [3 4 1 5 6 2]).*permute(alpha, [1 3 4 5 2]); 


Dsales_income__XCL = zeros(1,1,N, N,S,N);
Dsales_income__XCL = 1./permute(WAGE.*L, [3 2 1]) .* ones(1,1,N).*permute(ones(N,1,N),[4 5 6 1 2 3])...
    .*permute(eye(N), [3 4 1 5 6 2]).*permute(alpha, [1 3 4 5 2]); 


Dsales_income__XCD = Dsales_income__XCL;


Dsales_income = zeros(1,1,N, totvarnumb);
Dsales_income(:,:,:, WAGEfrom:WAGEto) = reshape(- 1./permute(WAGE.^2.*L, [3 2 1]) .* permute(eye(N), [3 4 1 2]) .* (sum(alpha.*sum(XCL+XCD+(1-1./sig).*XM,1),2) ...
    + permute(WAGE, [3 2 1]).*sum(RDlabor_per+sum(f_dest.*ETA.*RDlabor.^(1-kappa).*permute(THOLD, [3 2 1]).^-k + f_origin.*ETA_per.*RDlabor_per.^(1-kappa).*THOLD.^-k,1),2)) ...
    + 1./permute(WAGE.*L, [3 2 1]) .* permute(eye(N), [3 4 1 2]).* ( ...
    permute(sum(sum(f_dest.*ETA_per.*RDlabor_per.^(1-kappa).*THOLD.^-k,3),2), [3 2 1]) ...
    + sum(sum(f_origin.*ETA_per.*RDlabor_per.^(1-kappa).*THOLD.^-k,1)+RDlabor_per,2)) ...
    , 1,1,N,[]); 
% permute((-L+sum(sum(f_dest.*ETA_per.*RDlabor_per.^(1-kappa).*THOLD.^-k+permute(f_origin.*ETA_per.*RDlabor_per.^(1-kappa).*THOLD.^-k, [3 2 1]), 3) + RDlabor,2)) .* permute(eye(N), [1 3 4 2]) ...
%         , [2  3  1  4]) % this may be the more logical formulation, and it coincides with the above (thus above seems fine).

Dsales_income(:,:,:, RDlaborfrom:RDlaborto) = 1./permute(WAGE.*L, [3 2 1]) .* (reshape(permute(WAGE, [3 2 1]) ...
    .*(1+sum((permute(f_dest.*ETA_per(:,2:end,:).*(1-kappa).*RDlabor_per(:,2:end,:).^-kappa.*THOLD(:,2:end,:).^-k, [4 5 6 1 2 3]) ...
    + permute(sum(f_origin.*ETA_per(:,2:end,:).*(1-kappa).*RDlabor_per(:,2:end,:).^-kappa.*THOLD(:,2:end,:).^-k,1).*permute(eye(N), [1 3 2]), [4 5 6 1 2 3])) ...
    .*permute(eye(N), [3 4 1 2]),4)), 1,1,N,[]));

Dsales_income(:,:,:, RDlaborfrom:RDlaborto) = 1./permute(WAGE.*L, [3 2 1]) .* (reshape(permute(WAGE.*permute(eye(N), [1 3 4 5 6 2])+ WAGE.*permute(eye(N), [1 3 4 5 6 2]) ...
     .*permute( sum( f_origin.*ETA_per(:,2:end,:).*(1-kappa).*RDlabor_per(:,2:end,:).^-kappa.*THOLD(:,2:end,:).^-k,1) , [4 5 6 1 2 3]) ...
        + sum(permute(eye(N), [1 3 4 2]).*WAGE.*permute(f_dest.*ETA_per(:,2:end,:).*(1-kappa).*RDlabor_per(:,2:end,:).^-kappa.*THOLD(:,2:end,:).^-k, [4 5 6 1 2 3]), 4) ...
         , [2 3 1 4 5 6]), 1,1,N,[])); % IS THIS BETTER THAN ABOVE? GO ON HERE!!! (ITS DIFFERENT FOR SURE) BUT CORRECT??.
    
    
Dsales_income(:,:,:, THOLDfrom:THOLDto) = reshape( 1./permute(WAGE.*L, [3 2 1]) .* (...
    permute(permute(WAGE, [3 2 1]) ...
    .*f_dest.*ETA_per(:,2:end,:).*RDlabor_per(:,2:end,:).^(1-kappa).*(-k).*THOLD(:,2:end,:).^(-k-1), [4 5 6 1 2 3]).*permute(eye(N), [3 4 1 2]) ...
    + permute(permute(WAGE, [3 2 1]) ...
    .*(f_origin).*ETA_per(:,2:end,:).*RDlabor_per(:,2:end,:).^(1-kappa).*(-k).*THOLD(:,2:end,:).^(-k-1), [4 5 6 1 2 3]).*permute(eye(N), [3 4 1 5 6 2])) ...
    , 1,1,N,[]); % old notation..think this was wrong. see below for new approach..can delete this eventually.

Dsales_income(:,:,:, THOLDfrom:THOLDto) = reshape( 1./permute(WAGE.*L, [3 2 1]) .* (...
    permute(WAGE, [3 2 1]) .* permute( ...
    f_dest.*ETA_per(:,2:end,:).*RDlabor_per(:,2:end,:).^(1-kappa).*(-k).*THOLD(:,2:end,:).^(-k-1), [4 5 6 1 2 3]).*permute(eye(N), [3 4 1 2]) ...
    + permute(WAGE, [3 2 1]) .* permute( ...
    f_origin.*ETA_per(:,2:end,:).*RDlabor_per(:,2:end,:).^(1-kappa).*(-k).*THOLD(:,2:end,:).^(-k-1), [4 5 6 1 2 3]).*permute(eye(N), [3 4 1 5 6 2])) ...
    , 1,1,N,[]); % new approach. CONT HEREXXXXXXXX

%    reshape(permute( WAGE.*permute(f_dest.*ETA_per(:,2:end,:).*RDlabor_per(:,2:end,:).^(1-kappa).*(-k).*THOLD(:,2:end,:).^(-k-1), [4 5 6 1 2 3]).*permute(eye(N), [1 3 4 2]) ...
%     + WAGE.*permute(f_origin.*ETA_per(:,2:end,:).*RDlabor_per(:,2:end,:).^(1-kappa).*(-k).*THOLD(:,2:end,:).^(-k-1), [4 5 6 1 2 3]).*permute(eye(N), [1 3 4 5 6 2]) ... 
%     , [2 3 1 4 5 6]), 1,1,N,[]) % same as above -- so can ignore once all is correct.


Dsales_income = Dsales_income + sum(sum(sum( ...
    + permute(Dsales_income__XM, [1 2 3 7 4 5 6]).*permute(DXM, [5 6 7 4 1 2 3]) ...
    + permute(Dsales_income__XCL, [1 2 3 7 4 5 6]).*permute(DXCL, [5 6 7 4 1 2 3]) ...
    + permute(Dsales_income__XCD, [1 2 3 7 4 5 6]).*permute(DXCD, [5 6 7 4 1 2 3]), 7), 6), 5);




% ----------------------------------------------------------------------- %
%   2. TB, i.e. TRADE (IM)BALANCE                                         %
% ----------------------------------------------------------------------- %
DTB__XM = ones(N,1,1,N,S+1,N) ...
    .*(permute(eye(N), [1 3 4 5 6 2])./(sum(sum(sum(XCL+XCD+XM,3),2),1)) - TB_numerator./(sum(sum(sum(XCL+XCD+XM,3),2),1)).^2);
% DTB__XM = (1./(sum(sum(sum(XCL+XCD+XM,3),2),1))).*ones(1,1,1,N,S+1,N) ...
%     .*(permute(eye(N), [1 3 4 5 6 2]).*1./(sum(sum(sum(XCL+XCD+XM,3),2),1)) - TB_numerator./(sum(sum(sum(XCL+XCD+XM,3),2),1)).^2); % old notation - but seems wrong. can eventually delete here.
DTB__XCL = DTB__XM;
DTB__XCD = DTB__XM;


DTB__P = -Y./(sum(sum(sum(XCL+XCD+XM,3),2),1)).*permute(eye(N),[1 3 4 2]);

DTB = zeros(N,1,1,totvarnumb);
DTB(:,:,:,Yfrom:Yto) = -P./(sum(sum(sum(XCL+XCD+XM,3),2),1)).*permute(eye(N),[1 3 4 2]);


DTB(:,:,:,WAGEfrom:WAGEto) = ( permute(sum(sum(f_dest.*ETA_per.*RDlabor_per.^(1-kappa).*THOLD.^-k, 3),2), [4 5 6 1 2 3]) ...
    .*permute(eye(N), [1 3 4 2]) ...
    - permute(sum(f_dest.*ETA_per.*RDlabor_per.^(1-kappa).*THOLD.^-k,2), [3 2 4 1]) ) ...
    ./ (sum(sum(sum(XCL+XCD+XM,3),2),1));


DTB(:,:,:,RDlaborfrom:RDlaborto) = (1./(sum(sum(sum(XCL+XCD+XM,3),2),1))) ...
    .*reshape(WAGE.*f_dest.*permute(sum(ETA_per(:,2:end,:).*(1-kappa) ...
    .*RDlabor_per(:,2:end,:).^-kappa.*THOLD(:,2:end,:).^-k ...
    .*permute(eye(N),[1 3 2]),1), [4 5 6 1 2 3]).*permute(eye(N), [1 3 4 5 6 2]), N,1,1,[]); %old stuff..seems wrong, but lets see. can delete eventually

DTB(:,:,:,RDlaborfrom:RDlaborto) = reshape(-permute(sum(WAGE.*f_dest.*ETA_per(:,2:end,:).*(1-kappa).*RDlabor_per(:,2:end,:).^-kappa.*THOLD(:,2:end,:).^-k, 1) ...
    , [4 5 6 1 2 3]) ./ (sum(sum(sum(XCL+XCD+XM,3),2),1)) .* permute(eye(N), [1 3 4 5 6 2]) ...
    + permute(WAGE.*f_dest.*ETA_per(:,2:end,:).*(1-kappa).*RDlabor_per(:,2:end,:).^-kappa.*THOLD(:,2:end,:).^-k, [1 4 5 6 2 3]) ./ (sum(sum(sum(XCL+XCD+XM,3),2),1)) ...
    , N,1,1,[]); % new approach -- but so far dont know if this is correct either. PROBABLY CONTINUE HERE XXXXXX.


DTB(:,:,:,THOLDfrom:THOLDto) = (1./(sum(sum(sum(XCL+XCD+XM,3),2),1))) ...
    .*reshape(permute(-WAGE.*f_dest.*ETA_per(:,2:end,:).*RDlabor_per(:,2:end,:).^(1-kappa).*(-k).*THOLD(:,2:end,:).^(-k-1) ...
    , [4 5 6 1 2 3]).*permute(eye(N), [1 3 4 5 6 2]) ...
    + permute(WAGE.*f_dest.*ETA_per(:,2:end,:).*RDlabor_per(:,2:end,:).^(1-kappa).*(-k).*THOLD(:,2:end,:).^(-k-1)...
    , [4 5 6 1 2 3]).*permute(eye(N), [1 3 4 2]), N,1,1,[]);


DTB = DTB + sum(sum(sum( ...
      permute(DTB__XM,  [1 2 3 7 4 5 6]).*permute(DXM,  [5 6 7 4 1 2 3]) ...
    + permute(DTB__XCL, [1 2 3 7 4 5 6]).*permute(DXCL, [5 6 7 4 1 2 3]) ...
    + permute(DTB__XCD, [1 2 3 7 4 5 6]).*permute(DXCD, [5 6 7 4 1 2 3]), 7), 6), 5) ...
    + sum(permute(DTB__P, [1 2 3 5 4]) .* permute(DP, [2 3 5 4 1]),5);


% ----------------------------------------------------------------------- %
%   3. RD equations                                                       %
% ----------------------------------------------------------------------- %      
DRD__Profnorm = 1./RDlabor_per(:,2:end,:).^kappa .* (permute(eye(N), [3 4 1 5 6 2]).*permute(eye(S), [3 1 4 5 2]) ...
    .*ETA_per(:,2:end,:)./(k-1).*k./(r+zeta(:,2:end,:)+nu(:,2:end,:)-g+gPsi_s(:,2:end,:)).*ones(1,1,1,N));


DRD__g = 1./RDlabor_per(:,2:end,:).^kappa .* ETA_per(:,2:end,:)./(k-1).*sum(Profnorm(:,2:end,:).*(-k)./(r+zeta(:,2:end,:)+nu(:,2:end,:)-g+gPsi_s(:,2:end,:)).^2.*(1./Gamma-1),1);


DRD__gPsi = 1./RDlabor_per(:,2:end,:).^kappa .* permute(eye(S), [3 1 4 5 2]).*ETA_per(:,2:end,:)./(k-1).*sum(Profnorm(:,2:end,:).*(-k)./(r+zeta(:,2:end,:)+nu(:,2:end,:)-g+gPsi_s(:,2:end,:)).^2,1);

DRD = zeros(1,S,N, totvarnumb);
DRD(:,:,:, THOLDfrom:THOLDto) = reshape(1./RDlabor_per(:,2:end,:).^kappa .* ETA_per(:,2:end,:)./(k-1) ...
    .*permute(-k.*THOLD(:,2:end,:).^(-k-1).*(f_origin+f_dest.*WAGE./permute(WAGE,[3 2 1])), [4 5 6 1 2 3]) ...
    .*permute(eye(N), [3 4 1 5 6 2]).*permute(eye(S), [3 1 4 5 2]), 1,S,N, []);


DRD(:,:,:, WAGEfrom:WAGEto) = 1./RDlabor_per(:,2:end,:).^kappa .* (ETA_per(:,2:end,:)./(k-1).*(f_dest.*1./permute(WAGE,[3 2 1]).*permute(THOLD(:,2:end,:).^-k, [4 2 3 1])) ...
    + ETA_per(:,2:end,:)./(k-1).*sum(f_dest.*-WAGE./permute(WAGE,[3 2 1]).^2.*THOLD(:,2:end,:).^-k,1).*permute(eye(N), [3 4 1 2]));


DRD(:,:,:, RDlaborfrom:RDlaborto) = reshape(-kappa .* RDlabor_per(:,2:end,:).^(-kappa-1) .* ETA_per(:,2:end,:)./(k-1).*sum(Profnorm(:,2:end,:).* ...
    k./(r+zeta(:,2:end)+nu(:,2:end,:)-g+gPsi_s(:,2:end,:)) + (f_origin+f_dest.*WAGE./permute(WAGE,[3 2 1])) .* THOLD(:,2:end,:).^-k, 1) ...
    .*permute(eye(N), [3 4 1 5 6 2]).*permute(eye(S), [3 1 4 5 2]), 1,S,N, []);


DRD = DRD + sum(sum(sum(permute(DRD__Profnorm, [1 2 3 7 4 5 6]).*permute(DProfnorm(:,2:end,:,:), [5 6 7 4 1 2 3]), 7), 6), 5) ...
    + DRD__g.*permute(Dg, [2 3 4 1]) + sum(DRD__gPsi.*permute(Dg_Psi_s, [3 4 5 2 1]), 5);


% ----------------------------------------------------------------------- %
%   4. patent threshold equations                                         %
% ----------------------------------------------------------------------- %
DPatThold__Profnorm = zeros(N,S,N, N,S,N);
DPatThold__Profnorm = 1./THOLD(:,2:end,:) .* -(f_origin+f_dest.*WAGE./permute(WAGE,[3 2 1]))...
    .*(r+zeta(:,2:end,:)-g+gPsi_s(:,2:end,:)+DELTA(:,2:end,:)).*(r+zeta(:,2:end,:)+nu(:,2:end,:)-g+gPsi_s(:,2:end,:)+DELTA(:,2:end,:)) ...
    ./ (nu(:,2:end,:).*Profnorm(:,2:end,:).^2) ...
    .*permute(eye(N), [1 3 4 2]).*permute(eye(N), [3 4 1 5 6 2]).*permute(eye(S), [3 1 4 5 2]);
    

DPatThold = zeros(N,S,N, totvarnumb);
DPatThold(:,:,:, THOLDfrom:THOLDto) = reshape(-1./THOLD(:,2:end,:).^2 .* (f_origin+f_dest.*WAGE./permute(WAGE,[3 2 1]))...
    .*(r+zeta(:,2:end,:)-g+gPsi_s(:,2:end,:)+DELTA(:,2:end,:)).*(r+zeta(:,2:end,:)+nu(:,2:end,:)-g+gPsi_s(:,2:end,:)+DELTA(:,2:end,:)) ...
    ./ (nu(:,2:end,:).*Profnorm(:,2:end,:)) ...
    .*permute(eye(N), [1 3 4 2]).*permute(eye(N), [3 4 1 5 6 2]).*permute(eye(S), [3 1 4 5 2]), N,S,N, []);


DPatThold__g = zeros(N,S,N,1);
DPatThold__g = 1./THOLD(:,2:end,:) .* (f_origin+f_dest.*WAGE./permute(WAGE, [3 2 1]))./(Profnorm(:,2:end,:).*nu(:,2:end,:)) ...
    .*((1./Gamma-1).*(r+zeta(:,2:end,:)+nu(:,2:end,:)-g+gPsi_s(:,2:end,:)+DELTA(:,2:end,:)) ...
    +(1./Gamma-1).*(r+zeta(:,2:end,:)-g+gPsi_s(:,2:end,:)+DELTA(:,2:end,:)));


DPatThold__gPsi = 1./THOLD(:,2:end,:) .* (f_origin+f_dest.*WAGE./permute(WAGE, [3 2 1]))./(Profnorm(:,2:end,:).*nu(:,2:end,:)) ...
    .*((r+zeta(:,2:end,:)-g+gPsi_s(:,2:end,:)+DELTA(:,2:end,:))+(r+zeta(:,2:end,:)+nu(:,2:end,:)-g+gPsi_s(:,2:end,:)+DELTA(:,2:end,:))) ...
    .*permute(eye(S), [3 1 4 5 2]);


DPatThold(:,:,:, WAGEfrom:WAGEto) = 1./THOLD(:,2:end,:) .* (f_dest./permute(WAGE,[3 2 1]) ...
    .*(r+zeta(:,2:end,:)-g+gPsi_s(:,2:end,:)+DELTA(:,2:end,:)).*(r+zeta(:,2:end,:)+nu(:,2:end,:)-g+gPsi_s(:,2:end,:)+DELTA(:,2:end,:)) ...
    ./(Profnorm(:,2:end,:).*nu(:,2:end,:)) .* permute(eye(N), [1 3 4 2]) ...
    - WAGE./permute(WAGE,[3 2 1]).^2.*f_dest ...
    .*(r+zeta(:,2:end,:)-g+gPsi_s(:,2:end,:)+DELTA(:,2:end,:)).*(r+zeta(:,2:end,:)+nu(:,2:end,:)-g+gPsi_s(:,2:end,:)+DELTA(:,2:end,:)) ...
    ./(Profnorm(:,2:end,:).*nu(:,2:end,:)) .* permute(eye(N), [3 4 1 2]));


DPatThold = DPatThold + sum(sum(sum(permute(DPatThold__Profnorm, [1 2 3 7 4 5 6]).*permute(DProfnorm(:,2:end,:,:), [5 6 7 4 1 2 3]), 7), 6), 5) ...
    + DPatThold__g.*permute(Dg, [2 3 4 1]) ...
    + sum(DPatThold__gPsi.*permute(Dg_Psi_s, [3 4 5 2 1]), 5);





% ----------------------------------------------------------------------- %
%   5. price ratios                                       %
% ----------------------------------------------------------------------- %
DPriceRatio = zeros(N-1,1,1, totvarnumb); 
DPriceRatio(:,:,:, Pfrom:Pto) = -permute(eye(N-1), [1 3 4 2]);

DNumOverPhi = (sig./(sig-1)).^(1-sig).*PsiM.*(sig-1)./theta.*Phi.^((sig-1)./theta-1) + PsiCL.*(sig-1)./theta.*Phi.^((sig-1)./theta-1) ...
    + PsiCD.*(sig-1)./theta.*sum(Phi,3).^((sig-1)./theta-1);
DPriceRatio__Phi = permute([0; P(2:end)]./P(1) .* Beta./(1-sig).*1./denominator_allPoverP, [3 4 5 1 2]) ...
    .* permute(DNumOverPhi, [4 5 6 1 2 3]).*permute(eye(N), [1 3 4 2]) ...
    + permute([0; P(2:end)]./P(1) .* Beta./(1-sig).*-1./denominator_allPoverP(1,:), [1 3 4 5 2]) ... % I think this part is correct, but check again if derivative mistake!
    .* permute(DNumOverPhi(1,:,:), [4 5 6 1 2 3]).*permute([1; zeros(N-1,1)], [2 3 4 1]);


DPriceRatio__PsiM = permute([0; P(2:end)]./P(1) .* Beta./(1-sig).*1./denominator_allPoverP, [3 4 5 1 2]) ...
    .* permute((sig./(sig-1)).^(1-sig).*Phi.^((sig-1)./theta)-sum(Phi,3).^((sig-1)./theta), [4 5 6 1 2 3]).*permute(eye(N), [1 3 4 2]) ...
    + permute([0; P(2:end)]./P(1) .* Beta./(1-sig).*-1./denominator_allPoverP(1,:), [1 3 4 5 2]) ...% I think this part is correct, but check again if derivative mistake!
    .* permute((sig./(sig-1)).^(1-sig).*Phi(1,:,:).^((sig-1)./theta)-sum(Phi(1,:,:),3).^((sig-1)./theta), [4 5 6 1 2 3]).*permute([1; zeros(N-1,1)], [2 3 4 1]);


DPriceRatio__PsiCL = permute([0; P(2:end)]./P(1) .* Beta./(1-sig).*1./denominator_allPoverP, [3 4 5 1 2]) ...
    .* permute(Phi.^((sig-1)./theta)-sum(Phi,3).^((sig-1)./theta), [4 5 6 1 2 3]).*permute(eye(N), [1 3 4 2]) ...
    + permute([0; P(2:end)]./P(1) .* Beta./(1-sig).*-1./denominator_allPoverP(1,:), [1 3 4 5 2]) ...% I think this part is correct, but check again if derivative mistake!
    .* permute(Phi(1,:,:).^((sig-1)./theta)-sum(Phi(1,:,:),3).^((sig-1)./theta), [4 5 6 1 2 3]).*permute([1; zeros(N-1,1)], [2 3 4 1]);


DPriceRatio = DPriceRatio + sum(sum(sum(permute(DPriceRatio__Phi(2:end,:,:,:,:,:), [1 2 3 7 4 5 6]).*permute(DPhi_nsi, [5 6 7 4 1 2 3]), 7), 6), 5) ...
    + sum(sum(sum(permute(DPriceRatio__PsiM(2:end,:,:,:,2:end,:), [1 2 3 7 4 5 6]).*permute(DPsiM, [5 6 7 4 1 2 3]), 7), 6), 5) ...
    + sum(sum(sum(permute(DPriceRatio__PsiCL(2:end,:,:,:,2:end,:), [1 2 3 7 4 5 6]).*permute(DPsiC_L, [5 6 7 4 1 2 3]), 7), 6), 5);
    


% ----------------------------------------------------------------------- %
%   6. price normalization P1=1                                       %
% ----------------------------------------------------------------------- %
DPriceNorm__Phi1 = P(1).*Beta./(1-sig).*1./(numerator_PMoverP(1,:,:)+numerator_PCLoverP(1,:,:)+numerator_PCDoverP(1,:,:)) ...
    .*DNumOverPhi(1,:,:);


DPriceNorm__PsiM = P(1).*Beta./(1-sig).*1./(numerator_PMoverP(1,:,:)+numerator_PCLoverP(1,:,:)+numerator_PCDoverP(1,:,:)) ...
    .*((sig./(sig-1)).^(1-sig).*Phi(1,:,:).^((sig-1)./theta)-sum(Phi(1,:,:),3).^((sig-1)./theta));


DPriceNorm__PsiCL = P(1).*Beta./(1-sig).*1./(numerator_PMoverP(1,:,:)+numerator_PCLoverP(1,:,:)+numerator_PCDoverP(1,:,:)) ...
    .*(Phi(1,:,:).^((sig-1)./theta)-sum(Phi(1,:,:),3).^((sig-1)./theta));


DPriceNorm = zeros(totvarnumb,1);
DPriceNorm = sum(sum(DPriceNorm__Phi1.*DPhi_nsi(1,:,:,:),3),2) ...
    + sum(sum(DPriceNorm__PsiM(:,2:end,:).*DPsiM(1,:,:,:),3),2) ...
    + sum(sum(DPriceNorm__PsiCL(:,2:end,:).*DPsiC_L(1,:,:,:),3),2);



%% Jacobian. Should be a matrix of size (totvarnumb)*(totvarnumb+1),  
% i.e., there is one more equation than endogenous variables because of the P1=1 normalization.

J = [1.*reshape(Dsales_income,N,[]).' ...
     1.*reshape(DTB,N,[]).' ... % eventually, can move this to objective function.
     1.*reshape(DRD,S*N,[]).' ...
     1.*reshape(DPatThold,N^2*S,[]).' ...
     1.*reshape(DPriceRatio,N-1,[]).' ...
     1.*reshape(DPriceNorm,1,[]).' ... 
     ];
 
%  J(RDlaborfrom:RDlaborto,:) = 0; % to shut off for J check



% XXXX CONT HERE. SO now it seems that up to Profnorm, derivatives could be correct (at least for element 7..which was the issue before; but maybe check others too). 
% however, the four first equil conditions are still wrong..in particular, for SalesIncome, its still element 7 that is wrong..so dig deeper from here.

tem = J(Yfrom:Yto,:);
J(Yfrom:Yto,:) =  J(Yfrom:Yto,:)./permute(P_input.^1, [1 2]);
J(Pfrom:Pto,:) = J(Pfrom:Pto,:) + permute(-Y(2:end)./P_input(2:end).^1, [1 2]) .*  tem(2:end,:);


%% FORMAT OF OBJECTIVE AND GRADIENT -- DEPENDING ON SOLVER
if obj_is_scalar==1
    J = J.';
    J = (sum(2*J.*dev,1))';
    dev = sum(dev.^2);   
end
 
%   J = J';

end

end








