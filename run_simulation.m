clc;clear;
%% Parameter tuning
disp('Setting parameters...')

%1: VisTex, 2: UIUC
experiment=2;

%Numbers of Rotations for the Scattering transform, |R| in the Thesis
filt_opt.L=4;

%Do not change anything below
if experiment==1
    load VisTex.mat
    nTxt=640; %number of texture patches
    class_size=16; %number of patches in each class
    filt_opt.J=4; %lowest-frequency band-pass
    filt_opt.precision='double'; %single is default in ScatNet
    scat_opt.oversampling = 1; %critical sampling of scattering output
    scat_opt.M = 2; %maximum path length
    Wop=wavelet_factory_2d([128 128],filt_opt, scat_opt); %generate scattering wavelets
    dec_level=3; %decomposition level of multiresolution approach
elseif experiment==2
    load UIUC.mat
    nTxt=200;
    for k=1:nTxt
        textures(1:128,1:128,k)=imresize(textures(:,:,k), 0.5);
    end
    textures=textures(1:128,1:128,:);
    class_size=10;
    filt_opt.J=4; filt_opt.precision='double';
    scat_opt.oversampling = 1;
    scat_opt.M = 2;
    Wop=wavelet_factory_2d([128 128],filt_opt, scat_opt);
    dec_level=3;
end

%% Normalization

disp('Normalizing textures...')

for k=1:nTxt
    x=textures(:,:,k);
    x=x-mean(x(:)); %centering
    x=x/norm(x(:),2); %uniform variance
    textures(:,:,k)=x;
end

clear inp h experiment;

%% Extract Parameters for Subband Decomposition
disp('Extracting Features via Subband Decomposition...')

alpha_beta=ones(dec_level*3,2,nTxt); %GGD parameters for all textures
for k=1:nTxt
    alpha_beta(:,:,k)=fwt2alphabeta(textures(:,:,k),dec_level);
end

%% Extract Parameters for DT-CWT

disp('Extracting Features via DT_CWT...')

lambda_k_cwt=ones(dec_level*6,2,nTxt); %WD parameters for all textures
for k=1:nTxt
     [~,Yh] = dtwavexfm2(textures(:,:,k), dec_level,'antonini','qshift_a');
     for k_=1:dec_level
         for l=1:6
             x=Yh{k_}(:,:,l);
             lambda_k_cwt((k_-1)*6+l,:,k)=wblfit(abs(x(:)));
         end
     end
end

clear k k_ l x Yh;

%% Extract Parameters for Scattering Transform

disp('Extracting Features via Scattering Transform...')

[S]=scat(textures(:,:,1),Wop); %perform a scattering transform to get cell array dimensions
M1=length(S{2}.signal); M2=length(S{3}.signal); 
if scat_opt.M>2
    M3=length(S{4}.signal);

else
    M3=0;
end
nSide=size(S{2}.signal{1},1);

lambda_k=zeros(M1+M2+M3,2,nTxt); %Weibull parameters for all other layers

L=nSide^2;
for k=1:nTxt
    S=scat(textures(:,:,k),Wop);
    x=S{1}.signal{1}; x=x(:);
    for l=1:M1
        x=S{2}.signal{l};
        lambda_k(l,:,k)=wblfit(x(:));
    end
    for l=1:M2
        x=S{3}.signal{l};
        lambda_k(l+M1,:,k)=wblfit(x(:));
    end
    for l=1:M3
        x=S{4}.signal{l};
        lambda_k(l+M1+M2,:,k)=wblfit(x(:));
    end
end
clear S

%% Extract Parameters for Normalized Scattering Transform

disp('Extracting Features via Normalized Scattering Transform...')

[S]=scat(textures(:,:,1),Wop); %perform a scattering transform to get cell array dimensions
M1=length(S{2}.signal); M2=length(S{3}.signal); 
if scat_opt.M>2
    M3=length(S{4}.signal);

else
    M3=0;
end
nSide=size(S{2}.signal{1},1);
lambda_k_nwst=zeros(M1+M2+M3,2,nTxt); %Weibull parameters for all other layers
L=nSide^2;
for k=1:nTxt
    S=scat(textures(:,:,k),Wop);
    S=normalize_scat(S, mean(mean(abs(textures(:,:,k)))) );
    x=S{1}.signal{1}; x=x(:);
    for l=1:M1
        x=S{2}.signal{l};
        lambda_k_nwst(l,:,k)=wblfit(x(:));
    end
    for l=1:M2
        x=S{3}.signal{l};
        lambda_k_nwst(l+M1,:,k)=wblfit(x(:));
    end
    for l=1:M3
        x=S{4}.signal{l};
        lambda_k_nwst(l+M1+M2,:,k)=wblfit(x(:));
    end
end

clear L S Wop k l nSide scat_opt textures x scattering_normalization M* filt_opt

%% Evaluate Bhattacharryya
disp('Computing Distances...')
ds_wst=zeros(nTxt,nTxt); 
ds_nwst=ds_wst; 
ds_dtcwt=ds_wst;
ds_fwt=ds_wst;
for k=1:nTxt
    for l=1:nTxt
         ds_wst(k,l)=bhatt_weibull(lambda_k(:,:,k),lambda_k(:,:,l));
         ds_nwst(k,l)=bhatt_weibull(lambda_k_nwst(:,:,k),lambda_k_nwst(:,:,l));
         ds_dtcwt(k,l)=bhatt_weibull(lambda_k_cwt(:,:,k),lambda_k_cwt(:,:,l));
         ds_fwt(k,l)=cross_entropy_ggd(alpha_beta(:,:,k),alpha_beta(:,:,l));
    end
end

clear k l  lambda_k  lambda_k_cwt lambda_k_nwst alpha_beta
clear l

%% Calculate Success Rates

disp('Computing retrieval rates...')

%assure that the same patch is not considered for retrieval rate
ds_wst=ds_wst+eye(nTxt)*(max(ds_wst(:))-min(ds_wst(:))); 
ds_nwst=ds_nwst+eye(nTxt)*(max(ds_nwst(:))-min(ds_nwst(:))); 
ds_dtcwt=ds_dtcwt+eye(nTxt)*(max(ds_dtcwt(:))-min(ds_dtcwt(:)));
ds_fwt=ds_fwt+eye(nTxt)*(max(ds_fwt(:))-min(ds_fwt(:)));

sr_wst=zeros(nTxt,1);
sr_nwst=sr_wst; sr_dtcwt=sr_wst; sr_fwt=sr_wst;
for k=1:nTxt
    curClass=ceil(k/class_size);
    
    [~,minInd]=sort(ds_wst(k,:));
    minInd=minInd(1:class_size-1);
    classEqual=ceil(minInd/class_size)==curClass;
    sr_wst(k)=sum(classEqual)/(class_size-1);
    [~,minInd]=sort(ds_nwst(k,:));
    minInd=minInd(1:class_size-1);
    classEqual=ceil(minInd/class_size)==curClass;
    sr_nwst(k)=sum(classEqual)/(class_size-1);
    [~,minInd]=sort(ds_dtcwt(k,:));
    minInd=minInd(1:class_size-1);
    classEqual=ceil(minInd/class_size)==curClass;
    sr_dtcwt(k)=sum(classEqual)/(class_size-1);
    [~,minInd]=sort(ds_fwt(k,:));
    minInd=minInd(1:class_size-1);
    classEqual=ceil(minInd/class_size)==curClass;
    sr_fwt(k)=sum(classEqual)/(class_size-1);
end
clear ds_*
%evaluate rates
sr_dtcwt_class=zeros(class_size,nTxt/class_size); sr_wst_class=sr_dtcwt_class; sr_nwst_class=sr_dtcwt_class; sr_fwt_class=sr_wst_class;
sr_wst_avg=round(10000*mean(sr_wst))/100;
sr_nwst_avg=round(10000*mean(sr_nwst))/100;
sr_dtcwt_avg=round(10000*mean(sr_dtcwt))/100;
sr_fwt_avg=round(10000*mean(sr_fwt))/100;
sr_dtcwt_class(:)=sr_dtcwt; sr_dtcwt_class=100*mean(sr_dtcwt_class,1);
sr_wst_class(:)=sr_wst; sr_wst_class=100*mean(sr_wst_class,1);
sr_nwst_class(:)=sr_nwst; sr_nwst_class=100*mean(sr_nwst_class,1);
sr_fwt_class(:)=sr_fwt; sr_fwt_class=100*mean(sr_fwt_class,1);
clear class* curClass
clear dec_level k minInd nTxt