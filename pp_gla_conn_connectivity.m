%% pp_gla_src_fooof
% correlate bandpass filtered (or via wavelets) pupil and MEG signals

clear
restoredefaultpath

% -------------------------
% VERSION 1: no pupil lag
% -------------------------
v = 1;
SUBJLIST  = 1:24;
freqoi    = 2.^(1:(1/4):7);
win_len = 800;
lag = 0;
overlap = 0.5;
% -------------------------
% VERSION 2: with pupil lag
% -------------------------
% v = 2;
% SUBJLIST  = 1:24;
% freqoi    = 2.^(1:(1/4):7);
% win_len = 800;
% lag = 1;
% overlap = 0.5;
% -------------------------

addpath('~/Documents/MATLAB/fieldtrip-20181231/')
addpath ~/pconn/matlab/
load(sprintf('~/pp/proc/pp_atlas_BNA.mat'))

ft_defaults

outdir = '~/pp_conn/proc/conn/';
ord    = pconn_randomization;

trans = pp_transfer_gla2hh;

%%
% -------------------------
for isubj = 1:24
  
  clear data dat pupil pup dataf src_r
  
  for iblock = 1:1
    %
    fn = sprintf('pp_gla_src_powerspectra_s%d_b%d_v%d',isubj,iblock,v);
    if tp_parallel(fn,outdir,1,0)
      continue
    end
    %
    fprintf('Processing subj%d block%d ...\n',isubj,iblock);
    
    try
      % load pupil data
      if isubj < 10
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/pupil/sub0%d_gla_pp.mat',isubj));
      else
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/pupil/sub%d_gla_pp.mat',isubj));
      end
      
      pupil = data.trial{1}';
      f_sample = data.fsample;
      
      % load meg data
      if isubj < 10
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/meg/sub0%d_gla_meg.mat',isubj));
      else
        load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/meg/sub%d_gla_meg.mat',isubj));
      end
      
      
      artifPnts=data.cfg.artfctdef.visual.artifact;
      
      cfg=[];
      cfg.layout='4D248.lay';
      lay = ft_prepare_layout(cfg);
      [~, outp.chanidx] = ismember(lay.label(1:248),data.label);
      
      chanidx = outp.chanidx;
      save(sprintf('~/pp/proc/src/chanidx_s%d.mat',isubj),'chanidx')

    catch me
      src_r = nan(246,25);
      save([outdir fn '.mat'],'src_r')
      continue
    end
    
    k = 2;
    fnq = f_sample/2;
    hil_hi = 0.005;
    hil_lo = 2;
    hil_Wn=[hil_hi/fnq hil_lo/fnq];
    [bhil, ahil] = butter(k, hil_Wn);
    
    pupil = filtfilt(bhil, ahil, pupil(:,4));
    
    if lag
        pup_shift = round(f_sample*0.93); % 930s from hoeks and levelt (1992?)
        pupil = pupil(pup_shift:end); pupil(end+1:end+pup_shift-1)=nan;
    end
    
    data.avg = data.trial{1}'; 

    if isubj < 10
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/leadfields/sub0%d_gla_lf_BNA5mm.mat',isubj))
    else
      load(sprintf('~/pp/data_gla/fw4bt/osfstorage/data/gla01/leadfields/sub%d_gla_lf_BNA5mm.mat',isubj))
    end
    
    
    for iart = 1 : size(artifPnts,1)
        data.avg(artifPnts(iart,1):artifPnts(iart,2),:)=NaN;
    end

    clear tp_csd
    
    for ifreq=1:numel(freqoi)
      
      fprintf('Freq: %d\n',ifreq)

      % -------------------------------
      % compute csd
      % -------------------------------
      para          = [];
      para.freq     = freqoi(ifreq);
      para.fsample  = f_sample;  
      para.overlap = 0.5;
      tp_csd(:,:,ifreq)=tp_compute_csd_wavelets(data.avg',para);
      
    end
    
    tp_csd = nanmean(tp_csd,3);
      
    % -------------------------------
    % beamforming
    % -------------------------------
    para      = [];
    para.iscs = 1;
    para.reg  = 0.05;
    filt   = tp_beamformer(real(tp_csd),lf,para);
    % -------------------------------
    
    opt.n_win = win_len; % 10s segment length, i.e., 0.1:0.1:100
    opt.n_shift = win_len*(1-overlap); % no overlap
    
    nseg=floor((size(data.avg,1)-opt.n_win)/opt.n_shift+1);
    clear pxx fxx pup pup_df
    if v<222
      ff = 2:1/(opt.n_win/400):128;
    else
      ff = 2:1/(800/400):128;  
    end
    
    pupil = pupil(1:size(data.avg,1));
    pup_nanidx = isnan(pupil);
    pupil_df = diff(pupil);
    
    data.avg(pup_nanidx,:)=nan;
     
    pxx = nan(size(ff,2),max(BNA.tissue_5mm(:)),nseg);
    
    for iseg = 1 : nseg
        fprintf('%d / %d\n',iseg,nseg)
        seg_dat = data.avg((iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win,:)*filt;
        
        seg_dat = seg_dat(:,trans);
        
        if any(isnan(seg_dat(:,1)))
%             pxx(:,:,iseg) = nan(size(fxx,1),size(tp_filt,2));
            pup(iseg) = nan;
            pup_df(iseg)=nan;
            continue        
        end
        if v<222
          [tmp_pxx,fxx]=pwelch(seg_dat,hanning(opt.n_win),0,ff,400,'power');
        else
          [tmp_pxx,fxx]=pwelch(seg_dat,hanning(800),0.5,ff,400,'power');  
        end
        for igrid = 1 : max(BNA.tissue_5mm(:))
          pxx(:,igrid,iseg) = mean(tmp_pxx(:,BNA.tissue_5mm == igrid),2);
        end
        
        pup(iseg)  = mean(pupil((iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win));
        if iseg ~= nseg
            pup_df(iseg) = mean(pupil_df((iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win));
        else
            pup_df(iseg) = mean(pupil_df((iseg-1)*opt.n_shift+1:(iseg-1)*opt.n_shift+opt.n_win-1));            
        end
    end
%     pxx=single(pxx);
    save([outdir fn '.mat'],'pxx','fxx','pup','pup_df')
    tp_parallel(fn,outdir,0)
    
    clear pxx fxx pup pup_df
    
  end
end

exit
error('!')
%%
