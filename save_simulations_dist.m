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
[ug, sg, vg] = spm_svd(G2d_red*G2d_red',GainSVDTh);
UP = ug'; % direction of dimention reduction
G2dU_red = UP*G2d_red;
G2d0U_red = UP*G2d0_red;

[ug0, sg0, vg0] = spm_svd(G2d0*G2d0',GainSVDTh);

G2d0U = ug0'*G2d0;

G2d0U2 = G2d0U(:,1:2:end);
G2d02 = G2d0(:,2:2:end);
V2 = V(:,2:2:end);


ind_left = find(R(:,2)>0.01);
ind_right = find(R(:,2)<-0.01);





% 3. PROJECTION MATRICES
RankG = size(G2d0U_red,1);

% now find span of the target space
NSites = fix(size(G2d0U_red,2)/2);
Wks =  load('C_re.mat','C_re');
% this is correlation subspace correlation matrix
Ccorr = Wks.C_re;
   
% 4. PSIICOS projection for sparse matrix
%[Upwr, ds, Apwr] = ProjectorOnlyAwayFromPowerComplete(G2d0U_red, 1500, 0);
Rnk = 500;
load("UpwrApwr.mat") % using already calculated projector

Cpwr = Apwr*Apwr';
Wpwr = sqrtm(inv(Cpwr+0.001*trace(Cpwr)/(RankG^2)*eye(size(Cpwr))));
WCcorrW = Wpwr*double(Ccorr)*Wpwr';
ProjRnkMax = 1280;
[u, ~] = eigs(WCcorrW,ProjRnkMax);
UcorrW_rnk = u(:,1:ProjRnkMax);
[u, ~] = eigs(double(Ccorr),ProjRnkMax);
Ucorr_rnk = u(:,1:ProjRnkMax);

PrFromCorr_W = inv(Wpwr+0.05*trace(Wpwr)/size(Wpwr,1))*(eye(size(UcorrW_rnk,1))-UcorrW_rnk(:,1:500)*UcorrW_rnk(:,1:500)')*Wpwr;
PrPwr = Upwr(:,1:Rnk)*Upwr(:,1:Rnk)';







% 5. SIMULATIONS
d = 0.005:0.005:0.09;
snr = 4; % snr level in the data
c = lines(7); % colors
Fs = 500; % sampling frequency
srcF = 20; % sources frequency
Ntr = 100; % number of simulated trials
T = Fs; % number of time points in one trial
t = 1:T;
Nmc = 500;

Zp = zeros(2, length(d), Nmc, Nsites_red);
Zpw = zeros(2, length(d), Nmc, Nsites_red);
Zmne = zeros(2, length(d), Nmc, Nsites_red);
Zbf = zeros(2, length(d), Nmc, Nsites_red);
picked_src = zeros(Nmc, length(d), 2);

%Pre-computing MNE kernel
Cs = eye(Nsites_red*2);
Cn = eye(size(G2dU_red,1));
GGt = G2dU_red*Cs*G2dU_red';
lambda = 0.1;
Wmne = G2dU_red'/(G2dU_red*Cs*G2dU_red'+lambda*trace(GGt)/size(GGt,1)*Cn);

% load noise
load('NoiseAveragedNormalized1.mat')

% Example of noise generation
%range = 1:T;
%for tr = 1:Ntr
%    Noise(:,range) = GenerateBrainNoise(G2d0,T,200,500,Fs);
%    range = range+T;
%end
%Noise_av = mean(reshape(Noise,[Nch, T, Ntr]), 3); % average by trial
%Noise_av_0 = Noise_av/norm(Noise_av); % normalized average noise

% For each distance compute set of all sources that are dist apart from midline(Y axis), with an error of 1mm
dd = zeros(length(d), size(R, 1)); 
for i = 1:length(d)
	dSet = find((R(:, 2) > d(i)/2 - 0.001) & (R(:, 2) < d(i)/2 + 0.001));
    dSet = cat(1, dSet, zeros(size(R, 1) - size(dSet, 1), 1)); %fill with zeros
    dd(i, :) = dSet;
end
range = 1:T; % activation functions
for tr = 1:Ntr
	S(2,range) = sin(2*pi*srcF*t/Fs + randn*0);
	S(3,range) = cos(2*pi*srcF*t/Fs + randn*0);
	range = range + T;
end

dd_sym = zeros(1, Nsites);
for mc = 1:Nmc
    % generate signal for dense forward model matrix
    
    for dist_i = 1:length(d)
        % for random source in set find symmetrical source
        ddNonZero = find(dd(dist_i, :) ~= 0);
        rand_idx = randperm(length(ddNonZero));  %Pick random source from the set of possible for reqiured distance
        picked_src(mc, dist_i, 1) = dd(dist_i, rand_idx(1));  % pick the first source from the left hemisphere >2 cm far from midline
        clear dd_sym
        for i = 1:Nsites % find the symmetrical source
            dd_sym(i) = norm([R(picked_src(mc, dist_i,1),1),-R(picked_src(mc, dist_i,1),2), R(picked_src(mc, dist_i,1),3)]-R(i,:));
        end
        [val, picked_src(mc, dist_i,2)] = min(dd_sym);
        picked_src_oriented = picked_src(mc, dist_i,:)*2;
        
		range = 1:T; % activation functions
		for tr = 1:Ntr
			S(1,range) = sin(2*pi*srcF*t/Fs + randn*0.1);
			range = range + T;
		end
        for synch = 1:2
            % Data generation
            X = G2d0(:,picked_src_oriented(1))*S(1,:) + G2d0(:,picked_src_oriented(2))*S(synch+1,:);
            X_av = mean(reshape(X, [Nch, T, Ntr]), 3);
            X_av_0 = X_av/norm(X_av);
            % using pre-generated brain noise   Noise_av_0
            Data  = snr*X_av_0 + Noise_av_0; % add noise to the data
            Ca = UP*Data*Data'*UP'; % compute covariance in the virtual sensor space

            % LCMV BF
            Zbf(synch, dist_i, mc, :) = lcmv(G2dU_red, Ca);
            
            % Minimum norm estimate
            Zmne(synch, dist_i, mc, :) = mne(Wmne, Ca);

            % ReciPSIICOS beamformer
            Cap = reshape(PrPwr*Ca(:), size(Ca));
            [e, a] = eig(Cap);
            Cap = e*abs(a)*e';
            iCap = tihinv(Cap, 0.01);

            range2d = 1:2;
            for i=1:Nsites_red
                g = G2dU_red(:,range2d);
                m = inv(g'*iCap*g);
                [u, ss, v] = svd(m);
                Zp(synch, dist_i, mc, i) = ss(1,1);
                range2d = range2d+2;
            end
            
            % Whitened ReciPSIICOS beamformer
            Cap =reshape(PrFromCorr_W*Ca(:), size(Ca));
            [e, a] = eig(Cap);
            Cap = e*abs(a)*e';
            iCap = tihinv(Cap, 0.01);

            range2d = 1:2;
            for i=1:Nsites_red
                g = G2dU_red(:,range2d);
                m = inv(g'*iCap*g);
                [u, ss, v] = svd(m);
                Zpw(synch, dist_i, mc, i) = ss(1,1);
                range2d = range2d+2;
            end
        end
    end
    mc %#ok<NOPTS>
end

Z_total = zeros(4, 2, length(d), Nmc, Nsites_red);
Z_total(1, :, :, :, :) = Zp;
disp "Zp done"
Z_total(2, :, :, :, :) = Zpw;
disp "Zpw done"
Z_total(3, :, :, :, :) = Zbf;
disp "Zbf done"
Z_total(4, :, :, :, :) = Zmne;
disp "Zmne done"
save Z_totalDist Z_total
disp "Z saved"
save pickedSrcDist picked_src
disp "Sources saved"