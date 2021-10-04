clear all; close all; clc; 

cd matfig
redo = true; updatelist = true;
list=[];
while redo
    if updatelist
        contents=dir;
        contents=contents(find([contents.bytes]>0));
        nl=length(list);
        for i=1:length(contents)
            list(nl+i).name=contents(i).name;
            list(nl+i).path=contents(i).folder;
            list(nl+i).isfold=contents(i).isdir;
        end
        updatelist=false;
    end
    
    
    rm_ind=[];
    for i=length(list):-1:1
        if ~list(i).isfold
            np=split(list(i).name,'.');
            pp=split(list(i).path,'matfig');
            if strcmp(np{2},'fig')
                open(list(i).name);
                fname=strcat(pp{1},'pdf',pp{2},'/',np{1},'.pdf');
                print(fname,'-dpdf');
                close;
                
            end
            rm_ind=[rm_ind,i];
            k=length(contents);
            while k>0
                if strcmp(contents(k).name,list(i).name) && strcmp(contents(k).folder,list(i).path)
                    contents(k)=[]; 
                    k=length(contents);
                else
                    k=k-1;
                end
            end
        end
    end
    list(rm_ind)=[];
    
    if isempty(contents)
        cd ..
    else
        cd 
    end
    
    
end
