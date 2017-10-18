function make_Spcorr_JW(L,Jstr,Jdis,Pdist,Jseed)
%make_Spcorr_JW(L,Jstr,Jdis,Pdist,Jseed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make_Spcorr_JW
% function to make spin corr from components
% S.S = SzSz + 2*SxSx
% 
% Andrew Goldsborough - 21/02/2017
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('printing Spcorr\n');

%open files
fnameSz = strcat('./Szcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Pdist),'_',num2str(Jseed),'_Szcorr_JW_PBC.txt');
Spcorr = importdata(fnameSz);

fnameSx = strcat('./Sxcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Pdist),'_',num2str(Jseed),'_Sxcorr_JW_PBC.txt');
corr = importdata(fnameSx);
Spcorr(:,3) = Spcorr(:,3) + 2*corr(:,3);

%open files to write to
fnameSp = strcat('./Spcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Pdist),'_',num2str(Jseed),'_Spcorr_JW_PBC.txt');
fprintf(strcat(fnameSp,'\n'));
fidSpcorr = fopen(fnameSp, 'w');

%print to file
for i=1:size(Spcorr,1)
    fprintf(fidSpcorr,'%d %d %.15e\n',Spcorr(i,1),Spcorr(i,2),Spcorr(i,3));
end

%close file
fclose(fidSpcorr);