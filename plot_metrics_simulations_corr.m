%close all
%clear all
% 1. FORWARD MODEL

G3 = load('R0015_G.mat'); % forward model matrix 20002 sources
G3_red = load('RR015_G_5k.mat'); % reduced forward model matrix 5001 sources
chans = load('R0015_chans.mat'); % channels

R = G3.GridLoc; % source location, dense matrix
R_red = G3_red.GridLoc; % source location reduced matrix
ChUsed = find(strcmp({chans.Channel.Type}, 'MEG GRAD')); % set to use gradiometers only

[~, Nsites] = size(G3.Gain(ChUsed,1:3:end));
[Nch, Nsites_red] = size(G3_red.Gain(ChUsed,1:3:end));

% 2. SIMULATIONS 

% Z_total = (method, corrN, mc, src)
% methods: ReciPSIICOS, WReciPSIICOS, LCMV, MNE
load('Z_totalCorr.mat');
load('pickedSrcCorr.mat')

src_left_red = find(R_red(:,2)>0);
src_right_red = find(R_red(:,2)<0);

% assess the distance between the original and the sparse grids
clear val;
for i = 1:size(picked_src, 1)
    d  = R(picked_src(i, 1),:) - R_red;
    [val(i, 1), ~] = min(sum(d.*d,2));
    d  = R(picked_src(i, 2),:) - R_red;
    [val(i, 2), ~] = min(sum(d.*d,2));
end
%hist(sqrt(val),100)

corrSet = 0:0.11:1; % corr of sources
Nmc = size(picked_src, 1);;
Range_frac = [0.65; 0.25];
methodSet = 1:4;

%clear r var
for corrN = 1:length(corrSet) 
    for mc = 1:Nmc
        ind_generated = picked_src(mc, :);
		for method = methodSet
			for frac = 1:size(Range_frac, 2)
				Z = squeeze(Z_total(method, corrN, mc, :))';
				[max_val, max_ind] = max(Z);
				Z(Z < Range_frac(synch, frac)*max_val) = 0;
	  
				clear dist_l dist_r
				Z_left = Z(src_left_red); % values from left hemisphere
				Z_right = Z(src_right_red); % values from right hemisphere
				[max_val_l, max_l] = max(Z_left); % maximum in the left hemisphere
				max_ind_l = src_left_red(max_l);
				[max_val_r, max_r] = max(Z_right);
				max_ind_r = src_right_red(max_r);
	

				if (length(Z_left(Z_left>0)) == 0)|(length(Z_right(Z_right>0)) == 0)
					r(method, corrN, frac, mc) = NaN;
					var(method, corrN, frac, mc) = NaN;
					r_red(method, corrN, frac, mc) = NaN;
				   
					xyz_vec_trueL(method, corrN, frac, mc, :) = [NaN, NaN, NaN];
					xyz_vec_estL(method, corrN, frac, mc, :) =   [NaN, NaN, NaN];
					xyz_vec_trueR(method, corrN, frac, mc, :) = [NaN, NaN, NaN];
					xyz_vec_estR(method, corrN, frac, mc, :) =   [NaN, NaN, NaN];
				else
					
					d = R_red - R(ind_generated(1),:);
					[~, ind_generated_red(1)] = min(sum (d.*d,2));
					d = R_red - R(ind_generated(2),:);
					[~, ind_generated_red(2)] = min(sum (d.*d,2));;
					
					r_red(method, corrN, frac, mc) = (norm(R_red(max_ind_l,:)-R_red(ind_generated_red(1),:))+ ...
						norm(R_red(max_ind_r,:)-R_red(ind_generated_red(2),:)))/2;
					
					r(method, corrN, frac, mc) = (norm(R_red(max_ind_l,:)-R(ind_generated(1),:))+ ...
						norm(R_red(max_ind_r,:)-R(ind_generated(2),:)))/2;
			   
					xyz_vec_trueL(method, corrN, frac, mc,:) = R(ind_generated(1),:) ;
					xyz_vec_estL(method, corrN, frac, mc,:) =  R_red(max_ind_l,:);
				  
					xyz_vec_trueR(method, corrN, frac, mc,:) = R(ind_generated(2),:) ;
					xyz_vec_estR(method, corrN, frac, mc,:) =  R_red(max_ind_r,:);
					
					Z_left_norm = Z_left./sum(Z_left);
					Z_right_norm = Z_right./sum(Z_right);

					coord_active_l = R_red(src_left_red(Z_left_norm>0),:);
					for i = 1:size(coord_active_l,1)
						dist_l(i) = norm(coord_active_l(i,:)-R_red(max_ind_l,:));
					end
					
					coord_active_r = R_red(src_right_red(Z_right_norm>0),:);
					for i = 1:size(coord_active_r,1)
						dist_r(i) = norm(coord_active_r(i,:)-R_red(max_ind_r,:));
					end
					
					var(method, corrN, frac, mc) = (sum(Z_left_norm(Z_left_norm>0).*dist_l)+ ...
						sum(Z_right_norm(Z_right_norm>0).*dist_r))/2;
				end
			end
		end      
        mc
    end
    currentCorrelation = corrSet(corrN)
end

CorrPlot

