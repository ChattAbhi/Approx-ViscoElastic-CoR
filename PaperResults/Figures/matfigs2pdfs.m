% A simple way to use this is to first open this file in the matlab editor,
% then navigate to the folder where .fig files are stored, and make that
% current folder. Then select the code below in the matlab editor ->
% right-click the selected code -> Click "Evaluate Selection".

clear all; close all; clc;
contents=dir; 
for i=1:length(contents)
    if contents(i).bytes>0 && ~contents(i).isdir
        fname=contents(i).name;
        spl_fname=split(contents(i).name,'.');
        nam = spl_fname{1};
        ext = spl_fname{2};
        if strcmp(ext,'fig')
            fig=openfig(fname);
            print(fig,strcat(nam,'.pdf'),'-dpdf');
            close(fig);
        end
    end
end