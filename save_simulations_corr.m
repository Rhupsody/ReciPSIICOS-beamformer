% 1. FORWARD MODEL
% for dense and sparse matrices

G3 = load('R0015_G.mat'); % forward model matrix 20002 sources
G3_red = load('RR015_G_5k.mat'); % reduced forward model matrix 5001 sources
chans = load('R0015_chans.mat'); % channels

R = G3.GridLoc; % source location, dense matrix
R_red = G3_red.GridLoc; % source location reduced matrix
ChUsed = find(strcmp({chans.Channel.Type}, 'MEG GRAD')); % set to use gradiometers only
Nch = length(ChUsed);

[G2d, G2d0, Nsites,V] = G3toG2_aux(G3, ChUsed);
[G2d_red, G2d0_red, Nsites_red] = G3toG2(G3_red, ChUsed);


% 2. REDUCING SENSOR SPACE for sparse matrix
GainSVDTh = 0.001; % 0.05 results into 47 eigensensors and makes it run faster but produces less contrasting subcorr scans
                   % for a more reliable preformance use 0.01 to get all the sensor on board but be ready to wait;
[ug sg vg] = spm_svd(G2d_red*G2d_red',GainSVDTh);
UP = ug'; % direction of dimention reduction
G2dU_red = UP*G2d_red;
G2d0U_red = UP*G2d0_red;

[ug0 sg0 vg0] = spm_svd(G2d0*G2d0',GainSVDTh);

G2d0U = ug0'*G2d0;

G2d0U2 = G2d0U(:,1:2:end);
G2d02 = G2d0(:,2:2:end);
V2 = V(:,2:2:end);


ind_left = find(R(:,2)>0.01);
ind_right = find(R(:,2)<-0.01);





% 3. PROJECTION MATRICES
RankG = size(G2d0U_red,1);

% now find span of the target space
bLoad = true;
Nsrc = size(G2d0U_red,2)/2;
Swp = zeros(4);
Swp(2,3) = 1;
Swp(3,2) = 1;
Swp(1,1) = 1; 
Swp(4,4) = 1;
RankG = size(G2d0U_red,1);
NSites = fix(size(G2d0U_red,2)/2);


% 4. PSIICOS projection for sparse matrix
%[Upwr, ds, Apwr] = ProjectorOnlyAwayFromPowerComplete(G2d0U_red, 1500, 0);
load("Upwr.mat") % using already calculated projector
Rnk = 500;
PrPwr = Upwr(:,1:500)*Upwr(:,1:500)';






% 5. SIMULATIONS
corr = 0:0.11:1; % corr of sources
snr = 1.5;
c = lines(7); % colors
Fs = 500; % sampling frequency
srcF = 1; % sources frequency
Ntr = 100; % number of simulated trials
T = Fs; % number of time points in one trial
t = 1:T;
Nmc = 500;

Zp = zeros(2, length(corr), Nmc, Nsites_red);
Zbf = zeros(2, length(corr), Nmc, Nsites_red);
picked_src = zeros(Nmc, 2);

% load picked_src
% load noise
load('NoiseAveragedNormalized1.mat')
far_sites = find(R(:,2)>0.02); % sources from left hemisphere >2 cm far from midline (6542/10001)
for corr_i = 1:length(corr)
    range = 1:T;
	for tr = 1:Ntr
		S(1 + corr_i,range) = sin(2*pi*srcF*t/Fs + (1 - corr(corr_i))* pi/2 + randn*0);
		range = range + T;
	end
end

dd_sym = zeros(1, Nsites);
for mc = 1:Nmc

    % using pre-generated brain noise
    % generate signal for dense forward model matrix

    rand_idx = randperm(length(far_sites));
    picked_src(mc,1) = far_sites(rand_idx(1));  % pick the first source from the left hemisphere >2 cm far from midline

    clear dd_sym
    for i = 1:Nsites % find the symmetrical source
        dd_sym(i) = norm([R(picked_src(mc,1),1),-R(picked_src(mc,1),2), R(picked_src(mc,1),3)]-R(i,:));
    end
    [val, picked_src(mc,2)] = min(dd_sym);
    picked_src_oriented = picked_src(mc,:)*2;

    for corr_i = 1:length(corr)
        range = 1:T; % activation functions
		for tr = 1:Ntr
			S(1,range) = sin(2*pi*srcF*t/Fs + randn*0.1);
			range = range + T;
		end
        
		X = G2d0(:,picked_src_oriented(1))*S(1,:) + G2d0(:,picked_src_oriented(2))*S(1 + corr_i,:);
		X_av = mean(reshape(X, [Nch, T, Ntr]), 3);

		X_av_0 = X_av/norm(X_av);
		Data  = snr*X_av_0 + Noise_av_0; % add noise to the data
		Ca = UP*Data*Data'*UP'; % compute covariance in the virtual sensor space

		% LCMV BF
		Zbf(1, corr_i, mc, :) = lcmv(G2dU_red, Ca);

		% ReciPSIICOS beamformer
		Cap = reshape(PrPwr*Ca(:), size(Ca));
		[e a] = eig(Cap);
		Cap = e*abs(a)*e';
		iCap = tihinv(Cap, 0.01);

		range2d = 1:2;
		for i=1:Nsites_red
			g = G2dU_red(:,range2d);
			m = inv(g'*iCap*g);
			[u ss v] = svd(m);
			Zp(1, corr_i, mc, i) = ss(1,1);
			range2d = range2d+2;
		end
    end
    mc
end

Z_total = zeros(4, 2, length(d), Nmc, Nsites_red);
Z_total(1, :, :, :, :) = Zp;
disp "Zp done"
Z_total(2, :, :, :, :) = Zpw;
disp "Zpw done"
Z_total(3, :, :, :, :) = Zbf;
disp "Zbf done"
Z_total(4, :, :, :, :) = Zmne;
save ZtotalCorr Z_total
save pickedSrcCorr picked_src
