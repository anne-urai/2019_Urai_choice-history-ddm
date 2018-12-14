clear all;
close all;

%% For each subject fits using 50 random seeds are estimated and saved in 
%% the 'fits' directory. These fits are here combined into a single Matlab
%% structure per subject, saved in fits2. These fits2 structures are later
%% used by superScript_Step2 that is the refinmenet stage.
ss=1:33;
ss(ss==6)=[];
for s=1:length(ss)
  files=dir(['fits/out' num2str(ss(s)) '_*']);
  clear gf;
  clear parms;
  gf=[];
  parms=[];
  fileIndex = find(~[files.isdir]);
    for kkk=1:length(fileIndex)
       [s kkk] 
       try
        load(['fits/' files(fileIndex(kkk)).name]);
catch
	disp('Did not load warning..');
end
        gf=[gf out2];
        parms=[parms;out1];
        clear out1 out2
    end
    
    outC(ss(s))=length(gf);
     save(['fits2/out' num2str(ss(s)) '.mat'],'gf','parms');
end