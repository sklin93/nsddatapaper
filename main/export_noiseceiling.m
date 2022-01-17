% notes:
% - only consider cases of 3 trials
% - use all the data for the main version (also do split-half for the split1 and split2 versions)
% - be careful about missing data (invalid voxels)
% - translate the vmetric to SNR
%
% history:
% - 2021/07/08 - add _split1 and _split2 versions.
% - 2020/07/14 - now use sq, mean, sqrt (std is unbiased for variance!!)
% - 2020/05/17 - re-run b3 (due to b3 update)

% load and define
a1 = load('~/NSD/data/nsddata/experiments/nsd/nsd_expdesign.mat');
nsess = [37 40 32 30 40 32 40 30];
hemis = {'lh' 'rh'};

% define
% Subject 1	[145 186 148]	[81 104 83]	227021	226601
% Subject 2	[146 190 150]	[82 106 84]	239633	239309
% Subject 3	[145 190 146]	[81 106 82]	240830	243023
% Subject 4	[152 177 143]	[85 99 80]	228495	227262
% Subject 5	[141 173 139]	[79 97 78]	197594	198908
% Subject 6	[152 202 148]	[85 113 83]	253634	259406
% Subject 7	[139 170 145]	[78 95 81]	198770	200392
% Subject 8	[143 184 139]	[80 103 78]	224364	224398
nslices1pt8 = [83 84 82 80 78 83 81 78];
nslices1pt0 = [148 150 146 143 139 148 145 139];

% do it
for subjix=1:1
  fprintf('*** subject %d\n',subjix);
  dir0 = sprintf('~/NSD/data/nsddata_betas/ppdata/subj%02d',subjix);

  % experimental design stuff
  ord = a1.masterordering(1:750*nsess(subjix));
  ordU = unique(ord);
  allixs = [];
  for qq=1:length(ordU)
    ix = find(ord==ordU(qq));
    if length(ix)==3
      allixs(:,end+1) = ix(:);
    end
  end

  % func1pt8mm
  bdirs = matchfiles([dir0 '/func1pt8mm/betas_*'])
  for pp=1:length(bdirs)
    fprintf('***** bdir %s\n',bdirs{pp}); tic;
    chunks = chunking(1:nslices1pt8(subjix),10);
    snr = single([]);
    snr1 = single([]);
    snr2 = single([]);
    for cc=1:length(chunks), fprintf('cc=%d,',cc);
    % for cc=4:5, fprintf('cc=%d,',cc);
      betas = single([]);
      for qq=nsess(subjix):-1:1, fprintf('qq=%d,',qq);     % backwards so we get the memory allocation
        % invalid voxels are all 0
        temp = h5read(sprintf('%s/betas_session%02d.hdf5',bdirs{pp},qq),'/betas',[1 1 chunks{cc}(1) 1],[Inf Inf range(chunks{cc})+1 Inf]);
        betas(:,:,:,:,qq) = single(temp);  % = cat(5,betas,single(temp));  % no need to /300
      end
      % betas: 81   104    10   750    37
      betas = calczscore(betas,4,[],[],0);  % per-sess voxelwise zscore, invalid voxels become NaN
      temp = reshape(betas(:,:,:,flatten(allixs)),size(betas,1),size(betas,2),size(betas,3),3,[]); % (x,y,10,3,num_unique_img_for_subjix)
      temp = std(temp,[],4).^2; % voxel variance (x,y,10,1,num_unique_img_for_subjix)
      vmetric = sqrt(nanmean(temp,5));  % we ignore NaNs that seep in, (x,y,10)
      vmetric1 = sqrt(nanmean(temp(:,:,:,:,1:2:end),5)); % temp(:,:,:,:,1:2:end) contains only odd img ids: (x,y,10,1,num_unique_img_for_subjix/2), vmetric1 (x,y,10)
      vmetric2 = sqrt(nanmean(temp(:,:,:,:,2:2:end),5)); % temp(:,:,:,:,1:2:end) contains only even img ids: (x,y,10,1,num_unique_img_for_subjix/2), vmetric2 (x,y,10)
      snr = cat(3,snr,translatevmetric(vmetric));
      snr1 = cat(3,snr1,translatevmetric(vmetric1));
      snr2 = cat(3,snr2,translatevmetric(vmetric2));
      %keyboard;
    end
    % keyboard;
    nsd_savenifti(snr,[1.8 1.8 1.8],sprintf('%s/ncsnr.nii.gz',bdirs{pp}),4/3);
    nsd_savenifti(snr1,[1.8 1.8 1.8],sprintf('%s/ncsnr_split1.nii.gz',bdirs{pp}),4/3);
    nsd_savenifti(snr2,[1.8 1.8 1.8],sprintf('%s/ncsnr_split2.nii.gz',bdirs{pp}),4/3);
    toc;
  end

end
