%4/25/18
%GES
%a work flow to weight netmats by heritability per node
%filenames
f1='/projects/niblab/scripts/twin/hertibility/HD_correlation_p.mat'
f2='/projects/niblab/scripts/twin/hertibility/HD_average_corr_MZ.txt'
f3='/projects/niblab/scripts/twin/hertibility/HD_average_corr_DZ.txt'
f4='/projects/niblab/scripts/twin/hertibility/HD_additive_gene.txt'
f5='/projects/niblab/scripts/twin/hertibility/HD_unique_env.txt'
f6='/projects/niblab/scripts/twin/hertibility/HD_common_env.txt'
f7='/projects/niblab/scripts/twin/hertibility/HD_corr_1-p.txt'
f8='/projects/niblab/scripts/twin/hertibility/HD_tvals.txt'

%where dem maps at?
group_maps='/projects/niblab/data/HCP_PTN820/groupICA/groupICA_3T_HCP820_MSMAll_d25.ica';     
%where dem timeseries at?
ts_dir='/projects/niblab/data/HCP_PTN820/node_timeseries/3T_HCP820_MSMAll_d25_ts2/hdiscord';
%load timeseries
ts=nets_load(ts_dir,0.72,1,4);
%need to parse baddies?
ts.DD=[1:25];
%generate netmat correlation matrix run/sub x nodes
netmats1=  nets_netmats(ts,1,'corr');      
%where are the twins? This is a text file version of the mat file, use Vest2Text to convert
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loc=load('/projects/niblab/scripts/twin/hertibility/hdisc.txt');%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%remove first line
loc(:,1) = []; 
%how many columns are there?
cols = size(loc)
cols=cols(2)
%new blank matrix to fill
all=[]
%for every column in the design matrix
for pair=1:cols
	x=find(loc(:,pair)) %find the index of the runs for both pairs
	A=x(1:4) %split in half 
	B=x(5:8)
	twinA=zeros([1,376]); %make a blank array
	for value=A
		twinA(value)=1; %fill blank array with 1s where the twin ts is located
	end
	twin1=twinA*netmats1 %use matrix algebra to multiple the position by the netmat to get a sum of the timeseries per subject across runs
	twin1=reshape(twin1,[25,25]).'; %reshape into a square matrix
	%do the same for the second twin
	twinB=zeros([1,376]);
	for value=B
		twinB(value)=1
	end
	twin2=twinB*netmats1;
	twin2=reshape(twin2,[25,25]);
	
	[rho,pval]=corr(twin1,twin2); %correlation between the two twins

	all=cat(3,all,rho); %save dat shit
end

save(f1,'all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zygosity=load('/projects/niblab/scripts/twin/hertibility/wd_ev_zyg.txt'); %load in the zygosity by ev%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MZ=[] %blank matrices to be filled with respective zygosity
DZ=[]

for zyg=1:length(zygosity) %for each entry in zygosity
	if zygosity(zyg) == 0 %if it is 0, it mono, add it to the MZ matrix
		disp(zygosity(zyg))
		MZ=cat(3,MZ,all(:,:,zyg));
	elseif zygosity(zyg) == 1 %if it is di, add it to the DZ
		disp(zygosity(zyg))
		DZ=cat(3,DZ,all(:,:,zyg));
	end
end 

MZ_final=mean(MZ,3); %average the correlation per node
DZ_final=mean(DZ,3);
save(f2,'MZ_final','-ascii');
save(f3,'DZ_final','-ascii');

additive_gene=2*(MZ_final-DZ_final); %subtract the average MZ node from the average DZ node and multiply by 2 (Falconer's forum for H^2)
unique_env=1-MZ_final;
common_env=MZ_final-additive_gene;

save(f4,'additive_gene','-ascii');
save(f5,'unique_env','-ascii');
save(f6,'common_env','-ascii');


rows=size(netmats1) %how many rows are there
rows=rows(1)

corr_mat=[]

for row=1:rows %for each row, get row (which is sub by correlation expanded)
	sub=netmats1(row,:); %get the run and nodes
	add_gene=reshape(additive_gene,[625,1]); %decompose the add gene into a vector
	corrected_sub_mat=add_gene*sub; %multiple the matrices
	corrected_sub = diag(corrected_sub_mat) %get the diagonal
	corrected_sub=corrected_sub.'; %transpose
	corr_mat=cat(1,corr_mat,corrected_sub); %add the row to the new matrix
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[p_uncorrected,p_corrected,tstatistic]=nets_glm(corr_mat,'/projects/niblab/data/HCP_PTN820/node_timeseries/3T_HCP820_MSMAll_d25_ts2/hdiscord/PTN_hgl.mat','/projects/niblab/data/HCP_PTN820/node_timeseries/3T_HCP820_MSMAll_d25_ts2/hdiscord/PTN_hgl.con',0,10000, '/projects/niblab/data/HCP_PTN820/node_timeseries/3T_HCP820_MSMAll_d25_ts2/hdiscord/PTN_hgl.grp');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(f7,'p_corrected','-ascii');

save(f8,'tstatistic','-ascii');

