%% Primary MS analyses.

%requirements:
%run smriprep/fmriprep
%have applied parcellation .annot using .sh 2native script
%extracted corresponding .dat cortical features using the prep_msn script

%% Setup 
close all;
clear all;
clc;
%%
paths.home =fullfile('D:\aya_msn');
paths.scripts =fullfile(paths.home,'scripts');
paths.data =fullfile(paths.home,'data_clean');
paths.parcs = fullfile(paths.home,'parcs');
            Yeo500overlap = readmatrix(fullfile(paths.parcs,'Yeo_500_overlap.txt')); 
            vonEcon500overlap = readmatrix(fullfile(paths.parcs,'vonEcon_500_overlap.txt')); 
            
addpath(paths.home,paths.scripts,paths.data,paths.parcs)
addpath(fullfile(paths.scripts2,'NIfTI+toolbox'));
addpath(fullfile(paths.scripts2,'toolbox_matlab_nifti'));
addpath(fullfile(paths.scripts2,'BrainNetViewer'));
  GM_parc = load_nii(fullfile(paths.parcs,'DK_308','500.aparc_cortical_consecutive.nii.gz'));                    

%% prepare spin permutation
RAS = readmatrix(fullfile(paths.parcs,'centroids_500.txt')); %coordinates centroids 
permID = rotate_parcellation(RAS(1:152,:),RAS(153:end,:),10000); %v slow
save(fullfile(paths.data('permID_DK308.mat'),'permID'));

%%
%Importing the data:

%Regional values for  are provided for:
load(fullfile(paths.data,'PARC500_GC.dat')) % gauss curvature (GC)
load(fullfile(paths.data,'PARC500_MC.dat')) % mean curvature (MC)
load(fullfile(paths.data,'PARC500_CI.dat')) % internal curvature index (CI)
load(fullfile(paths.data,'PARC500_FI.dat')) %fold index (FI)
load(fullfile(paths.data,'PARC500_GM.dat')) %grey matter volume (GM)
load(fullfile(paths.data,'PARC500_SA.dat')) % surface area (SA)
load(fullfile(paths.data,'PARC500_CT.dat')) %cortical thickness (CT)

%load(fullfile(paths.data,'group.dat'))
group = readmatrix(fullfile(paths.data,'groupidx.xlsx')); %aya = 2 con = 1
load(fullfile(paths.data,'age.dat'))
load(fullfile(paths.data,'sex.dat'))
load(fullfile(paths.data,'subject.dat'))
load(fullfile(paths.data,'QC.dat')) %quality control
load(fullfile(paths.data,'QC_labels.mat'))

nregs=308; % number of regions ALTER
nsubs=length(group); % number of subjects- 
pats=find(group==2);
cons=find(group==1);

%% z-score the inputs:
PARC500_CT_zscore=zscore(transpose(PARC500_CT));
PARC500_SA_zscore=zscore(transpose(PARC500_SA));
PARC500_GM_zscore=zscore(transpose(PARC500_GM));
PARC500_MC_zscore=zscore(transpose(PARC500_MC));
PARC500_GC_zscore=zscore(transpose(PARC500_GC));
PARC500_CI_zscore=zscore(transpose(PARC500_CI));
PARC500_FI_zscore=zscore(transpose(PARC500_FI));

% Create a cell for each subject with all of the required inputs:
clear subj_features7
for subj=1:nsubs
    subj_features7{1,subj}(:,1)=PARC500_CT_zscore(:,subj);
    subj_features7{1,subj}(:,2)=PARC500_SA_zscore(:,subj);
    subj_features7{1,subj}(:,3)=PARC500_GM_zscore(:,subj);
    subj_features7{1,subj}(:,4)=PARC500_MC_zscore(:,subj);
    subj_features7{1,subj}(:,5)=PARC500_GC_zscore(:,subj);
    subj_features7{1,subj}(:,6)=PARC500_CI_zscore(:,subj);
    subj_features7{1,subj}(:,7)=PARC500_FI_zscore(:,subj);
end

% Calculate the MS matrices by correlating all inputs and set the diagonal to zero:
for subj=1:nsubs
    subj_MSN_7{1,subj}=corr(transpose(subj_features7{1,subj}));
    subj_MSN_7{1,subj}(logical(eye(size(subj_MSN_7{1,subj})))) = 0;
end

clear meanMS_regional
% Get regional MS values
for subj=1:nsubs
    meanMS_regional(subj,:)=sum(subj_MSN_7{1,subj})./(nregs-1);
end

%% Global differences in morphometric similarity:

x1=age;
x2=sex;
X = [ones(size(x1)) x1 x2 x1.*x2];

% Calculate regional residuals:
clear myresid_region
for region=1:nregs
    y = meanMS_regional(:,region);
    [b,bint,resid]=regress(y,X);
    YFIT = b(1) + b(2)*x1 + b(3)*x2 + b(4)*x1.*x2;
    myresid_region(region,:)=y-YFIT;
end

x = reshape(myresid_region(:,cons),[nregs*length(cons),1]);
y = reshape(myresid_region(:,pats),[nregs*length(pats),1]);

%calculate distribution differences.
[k_h,k_p] = kstest2(x,y);

% Plot histograms of regional residuals for control subjects and patients:
figure
h1 = histogram(x);
hold on
h2 = histogram(y);
h1.Normalization = 'probability';
h1.BinWidth = 0.01;
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
legend('Controls','Ayahuasca')
xlim([-0.18 0.18])
xlabel('MS- regional residuals')
ylabel('Relative frequency')

% Calculate mean MS:
for subj=1:nsubs
    meanMS(subj)=mean(meanMS_regional(subj,:));
end

tbl = table(age,sex,group,transpose(meanMS));
tbl.sex = categorical(tbl.sex);
tbl.group = categorical(tbl.group);
lm = fitlm(tbl,'Var4~age*sex+group');
p_mean=lm.Coefficients{4,4} % p-value for the effect of group on mean MS

% Plot box plot:

for subj=1:nsubs
    myresid_region_mean(subj)=mean(myresid_region(:,subj));
end

figure
boxplot(myresid_region_mean,group,'notch','on','Labels',{'Controls','Ayahuasca',})
ylim([-3*10^(-3) 6*10^(-3)])
ylabel('Mean residual')

%% %% Regional differences in morphometric similarity:

%%The following code calculates a t-statistic for regional differences in morphometric similarity:

dummy=meanMS_regional;

clear mytstat mypval
for region=1:nregs
  tbl = table(age,sex,group,dummy(:,region));
  tbl.sex = categorical(tbl.sex);
  tbl.group = categorical(tbl.group);
  lm = fitlm(tbl,'Var4~age*sex+group');
  mytstat(region)=lm.Coefficients{4,3};
  mypval(region)=lm.Coefficients{4,4};
end

mytstat=transpose(mytstat);
mypval=transpose(mypval);
meantstat=  mytstat;

pvalue_fdr = mafdr(mypval,'BHFDR',1); % FDR corrected p-values
sigregs=find(pvalue_fdr<0.05); % list of the statistically significant regions
sigregs_render = mytstat; 
sigregs_render(find(pvalue_fdr>=0.05)) = 0; %same but for later t-stat plotting and comparison if you wish 

dlmwrite(fullfile(paths.data,'mytstat_Maast.dat'),mytstat)
dlmwrite(fullfile(paths.data,'mypval_Maast.dat'),mypval)


%% replication of analyses using TIV
tiv =QC(:,strcmp(data_QC_labels, 'TIV'));

% now regional ms
dummy=meanMS_regional;
clear mytstat_tiv mypval_tiv
for region=1:nregs
  tbl = table(tiv,age,sex,group,dummy(:,region));
  tbl.sex = categorical(tbl.sex);
  tbl.group = categorical(tbl.group);
  lm = fitlm(tbl,'Var5~tiv+age*sex+group');
  mytstat_tiv(region)=lm.Coefficients{5,3};
  mypval_tiv(region)=lm.Coefficients{5,4};
end

mytstat_tiv=transpose(mytstat_tiv);
mypval_tiv=transpose(mypval_tiv);
meantstat_tiv=  mytstat_tiv;

pvalue_fdr_tiv = mafdr(mypval_tiv,'BHFDR',1); % FDR corrected p-values
sigregs_tiv=find(pvalue_fdr_tiv<0.05); % list of the statistically significant regions
sigregs_render_tiv = mytstat_tiv; 
sigregs_render_tiv(find(pvalue_fdr_tiv>=0.05)) = 0; %same but for plotting

%get binarised index of fdr for jaccard. 
sigregs_binary(:,1) = sigregs_render >0 | sigregs_render < 0;
sigregs_binary(:,2) = sigregs_render_tiv > 0 | sigregs_render_tiv < 0;
getJaccard(sigregs_binary(:,1),sigregs_binary(:,2))
perm_sphere_p( mytstat_tiv,mytstat,permID,'Pearson'); %compare
[h,p] = corr(mytstat_tiv,mytstat);


%% %% The code below calculates a grid figure for changes in regional MS

%First you need to extract meanMS_regional (as calculated above) for the control subjects only 

meanMS_con=mean(meanMS_regional(cons,:),1);

%Plot the correlation between the regional mean control MS and the mean t-statistic:
figure
scatter_kde(meanMS_con',meantstat, 'filled', 'MarkerSize', 40,'MarkerEdgeColor',[0 0 0])
refline(0,0)
hold on
plot([0 0], ylim)
plot([0 0], ylim,'-b')
xlabel('Mean control MS')
ylabel('Mean t-statistic')
l = lsline ;
set(l,'LineWidth', 3,'Color', 'black');
xticks([-0.1:0.05:0.1])
% Calculate the percentage of scatter points in each quadrant:

a=0;
b=0;
c=0;
d=0;

xvalues=meanMS_con;
yvalues=meantstat;
tspin_p = perm_sphere_p(meanMS_con',meantstat,permID,'Pearson');
tspin_r =corr(meanMS_con',meantstat);

for ind=1:nregs
    xval=xvalues(ind);
    yval=yvalues(ind);
    if ((xval<0)&&yval>0) % a is top left quadrant
        a=a+1;
        elseif ((xval>0)&&yval>0) % b is top right quadrant
        b=b+1;
        elseif ((xval<0)&&yval<0) % c is bottom left quadrant
        c=c+1;
        elseif ((xval>0)&&yval<0) % d is bottom right quadrant
        d=d+1;
    end
end

a=a/nregs % percentage of scatter points in the top left quadrant
b=b/nregs % percentage of scatter points in the top right quadrant
c=c/nregs % percentage of scatter points in the bottom left quadrant
d=d/nregs % percentage of scatter points in the bottom right quadrant


%% Now understand how changes fall within Yeo networks:

networks=Yeo500overlap;
%networks=vonEcon500overlap;
% To plot box plot for a specific network/class:
yeo_labels = {'VIS','SM','DA','VA','L','FPN','DMN'};
%network_labels = {'Primary motor','Association','Association Secondary sensory','Primary sensory','Limbic','Insular'};
% To calculate t-statistics and p-values (for the von Economo class and Yeo network Tables in the SI in the paper):
%yeo network based stats
clear myt_yeo myp_yeo
for class=1:7
    myregions=find(networks==class);
    for subj=1:nsubs
        classreg(subj)=sum(meanMS_regional(subj,myregions));
    end

    tbl = table(age,sex,group,transpose(classreg));
    tbl.sex = categorical(tbl.sex);
    tbl.group = categorical(tbl.group);
    lm = fitlm(tbl,'Var4~age*sex+group');
    myt_yeo(class)=lm.Coefficients{4,3};
    myp_yeo(class)=lm.Coefficients{4,4};
end

% vonEcon networks

% To plot box plot for a specific network/class:
networks=vonEcon500overlap;
econ_labels = {'Prim motor','Asso1','Asso2', 'Sec sens','Prim sens','Limbic','Insular'};
% To calculate t-statistics and p-values (for the von Economo class and Yeo network Tables in the SI in the paper):
%yeo network based stats
clear myt_vonecon myp_vonecon
for class=1:7
    myregions=find(networks==class);
    for subj=1:nsubs
        classreg(subj)=sum(meanMS_regional(subj,myregions));
    end

    tbl = table(age,sex,group,transpose(classreg));
    tbl.sex = categorical(tbl.sex);
    tbl.group = categorical(tbl.group);
    lm = fitlm(tbl,'Var4~age*sex+group');
    myt_vonecon(class)=lm.Coefficients{4,3};
    myp_vonecon(class)=lm.Coefficients{4,4};
end
%% Make spider plots of networks 
%plot absolute t values
figure;
subplot(1,2,1)
spider_plot(abs(myt_yeo),...
        'AxesLimits', [repmat(0.1,1,length(myt_yeo));repmat(6,1,length(myt_yeo))],...
        'AxesLabels', yeo_labels,...
      'AxesInterval', 2,...
      'FillOption', 'on',...
    'FillTransparency', 0.1,...
    'AxesLabelsEdge', 'none',...
    'AxesRadial', 'off')
hold on; 
subplot(1,2,2)
spider_plot(abs(myt_vonecon),...
        'AxesLimits', [repmat(0.5,1,length(myt_vonecon));repmat(5,1,length(myt_vonecon))],...
        'AxesLabels', econ_labels,...
      'AxesInterval', 2,...
      'FillOption', 'on',...
    'FillTransparency', 0.1,...
    'AxesLabelsEdge', 'none',...
    'AxesRadial', 'off')


 

























