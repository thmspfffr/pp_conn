%% pp_conn_connectivity
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

% -------------------------
% VERSION 1
% -------------------------
v = 1;
lag = 1;
% include 28 subjects, as in pfeffer et al. (2018) plos biology
SUBJLIST = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
grid = 'cortex_lowres';
REG = 0.05;
% -------------------------

addpath ~/Documents/MATLAB/fieldtrip-20160919/
addpath ~/pconn/matlab/
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))

ft_defaults

outdir = '~/pp_conn/proc/conn/';
ord    = pconn_randomization;


if strcmp(grid,'aal_6mm')
  v_grid = 6;
elseif strcmp(grid, 'cortex_lowres')
  v_grid = 9;
end

freqoi=2.^(2:(1/4):6);
% 2-128 Hz as per Hipp et al. (2012) Nat Neurosci

%%
% -------------------------
for isubj = SUBJLIST
  
  % identify placebo condition (ord==1)
  im = find(ord(isubj,:)==1);
  
  for iblock = 1:2
    %
    fn = sprintf('pp_hh_conn_connectivity_s%d_b%d_v%d',isubj,iblock,v);
    if tp_parallel(fn,outdir,1,0)
      continue
    end
    %
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load cleaned meg data
      % data from 'pp_prepare_data.m'
      load(sprintf('~/pp/data/ham/pp_rest_s%d_b%d_v%d.mat',isubj,iblock,1))
    catch me
      continue
    end
    
    cfg=[];
    cfg.layout='CTF275.lay';
    lay = ft_prepare_layout(cfg);
    [~, outp.chanidx] = ismember(lay.label(1:275),label(startsWith(label,'M')));
    
    % bp-filter and resample pupil
    % ------
    k = 2; f_sample = 1000;
    fnq = f_sample/2;
    hil_hi = 0.005; hil_lo = 2;
    hil_Wn = [hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    pupil = filtfilt(bhil, ahil, pupil(:,end));
    pupil = resample(pupil,400,1000);
    % ------
    f_sample = 400;
    % align pupil and meg (at signal offset)
    % ------
    
    pupil = pupil(end:-1:1,:);
    dat = dat(:,end:-1:1);
    
    len = min([size(pupil,1) size(dat,2)]);
    if len/400 > 600
      len = 400*600;
    end
    
    dat = dat(:,1:len);
    pupil = pupil(1:len);
    
    dat = dat(:,end:-1:1);
    pupil = pupil(end:-1:1);
    % ------
    
    % pupil shift: 930 ms from hoeks & levelt (1992)
    if lag
      pup_shift = round(400*0.93);
      pupil = pupil(pup_shift:end); pupil(end+1:end+pup_shift-1)=nan;
    end
    
    pupil_df = diff(pupil);
    
    dat(:,isnan(pupil))=nan(size(dat,1),sum(isnan(pupil)));
    
    load(sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,im,iblock,v_grid));
    
    lf = sa.L_cortex_lowres;
    
    for ifreq = 1:length(freqoi)
      
      % COMPUTE SPATIAL FILTER ACROSS ALL SEGMENTS
      % this focuses the analysis on a time window
      [~,opt]=tp_mkwavelet(freqoi(ifreq),0.75,f_sample,0);
      
      % ------------
      % Compute cross spectrum (using wavelets)
      % ------------
      para          = [];
      para.freq     = freqoi(ifreq);
      para.fsample  = 400;
      para.overlap  = 0.5;
      [cs,dataf]   	= tp_compute_csd_wavelets(dat,para);
      % ------------
      % Compute spatial filter (LCMV)
      % ------------
      para          = [];
      para.reg      = REG;
      filt          = tp_beamformer(real(cs),lf,para);
      % ------------
      
      nwin = 60*400;
      nshift = nwin/3;
      nep=floor((size(dat,2)-nwin)/nshift+1);
      
      for iep = 1 : nep
        
        data_loc = dat(:,(iep-1)*nshift+1:(iep-1)*nshift+nwin);
        
        % compute local pupil size
        outp.loc_pup(iep,ifreq)=nanmean(pupil((iep-1)*nshift+1:(iep-1)*nshift+nwin));
        
        fprintf('Processing s%d m%d b%d f%d epoch%d  ...\n', isubj,im,iblock,ifreq,iep)
        % ------------
        % Compute cross spectrum (using wavelets)
        % ------------
        para          = [];
        para.freq     = freqoi(ifreq);
        para.fsample  = 400;
        para.overlap  = 0.5;
        [cs,dataf]   	= tp_compute_csd_wavelets(data_loc,para);
        % ------------
        % Project to source (using filter from above)
        % ------------
        src_dat = dataf'*filt;
        % ------------
        % Compute orthogonalized power enevelope correlations
        % ------------
            
        kk = 0;
        
        for jj=1:size(src_dat,1)
          
          datasf1=src_dat(jj,:);
          
          if any(isnan(datasf1(:,1)))
            warning('NaN detected')
            continue
          end
          
          kk=kk+1;
          
          for i1=1:size(filt,2)
            x1=datasf1(i1);
            x2=imag(datasf1*conj(x1)./abs(x1));
            y1=log10(abs(x1)^2);
            y2=log10(abs(x2).^2);
            if kk==1
              res1(i1,:)=y1*y2;
              res2(i1,:)=y1;
              res3(i1,:)=y2;
              res4(i1,:)=y1^2;
              res5(i1,:)=y2.^2;
            else
              res1(i1,:)=res1(i1,:)+y1*y2; % E[xy]
              res2(i1,:)=res2(i1,:)+y1;    % E[x]
              res3(i1,:)=res3(i1,:)+y2;    % E[y]
              res4(i1,:)=res4(i1,:)+y1^2;  % E[x^2]
              res5(i1,:)=res5(i1,:)+y2.^2; % E[y^2]
            end
          end
          
          %           datvar(j,:) = abs(datasf1).^2;
          
        end
        
        %         variance = var(datvar);
        % -------------------------------------
        res1=res1/kk;
        res2=res2/kk;
        res3=res3/kk;
        res4=res4/kk;
        res5=res5/kk;
        % -------------------------------------
        % resout is asymetrical (see hipp 2012), resout needs to be averaged
        resout=(res1-res2.*res3)./(sqrt((res4-res2.^2).*(res5-res3.^2))+eps);
        outp.powcorr(:,:,iep,ifreq)=tanh((atanh(resout)./2 + atanh(resout)'./2));
        
      end
    end
    
    save([outdir fn '.mat'],'outp')
    tp_parallel(fn,outdir,0)
    
    clear src_r all_nai outp
    
  end
end


error('!')

