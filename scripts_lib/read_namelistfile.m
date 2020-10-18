function A=read_namelistfile(filenamestr)
%% read the namelist file of fortran input
% -inputs:
% filenamestr   string, the full path of the namelist file
% -outputs:
% A             struct, the data contained in the namelist file

% This script function is updated from read_mfile(filenamestr), which is
% first edited by Shuying SUN on 2019/05/27, here is an updated version
% updated by Shuying SUN on 2019/05/29, to make it support comments and 
% other variable types.
% Contact: bangrsun@163.com, sunshuyingc@enn.cn
% ENN Group 1989-2019, all rights reserved.

fid=fopen(filenamestr,'r');
if(fid<=0)
    error(['Can not open file:', filenamestr]);
end
A={};
i=0;    % namelist index
j=0;    % variable index
linestr=fgetl(fid);
linenum=1;  % locate current line number
while ischar(linestr)
    trimedline=strtrim(linestr);
    %% remove comments, save effective string to effstr
    if(~isempty(strfind(trimedline,'!')))
    % this line contains comments
        tmpstrcell=strsplit(trimedline,'!');
        effstr=strtrim(tmpstrcell{1});
    else
    % this line does not contains comments
        effstr=trimedline;
    end
    %% handling the effstr
    if(~isempty(effstr))
        firstchar=effstr(1);
        % if the line is started with '&', it means this line is a start
        % line of a new namelist
        if(firstchar=='&')
            i=i+1;
            A{i}.namelist=effstr(2:end);
            A{i}.varname={};
            A{i}.vartype={};
            A{i}.varvalue={};
            j=0;
            sumstr='';
        elseif(firstchar=='/')
            newstrcell=strsplit(sumstr,',');
            for k=1:length(newstrcell)
                if(~isempty(strfind(newstrcell{k},'=')))
                    j=j+1;
                    assignstrcell=strsplit(newstrcell{k},'=');
                    A{i}.varname{j} =strtrim(assignstrcell{1});
                    A{i}.varvalue{j}=strtrim(assignstrcell{2});
                    A{i}.vartype{j} ='character';
                else
                    A{i}.varvalue{j}=strcat(A{i}.varvalue{j},',',strtrim(newstrcell{k}));
                end
            end
        else
            sumstr=strcat(sumstr,effstr);
        end
    end
    linestr=fgetl(fid);
    linenum=linenum+1;
end % while
fclose(fid);
%%
for i=1:length(A)
    for j=1:length(A{i}.varvalue)
        valuestr=A{i}.varvalue{j};
        if(strcmpi(valuestr,'F')||strcmpi(valuestr,'T'))
            A{i}.vartype{j}='logical';
        elseif((valuestr(1)=='+' || valuestr(1)=='-' || ...
                (abs(valuestr(1))>=48 && abs(valuestr(1))<=57)))
            A{i}.varvalue{j}=str2num(valuestr);
            if(~isempty(strfind(valuestr,'.')))
                A{i}.vartype{j}='double';
            else
                A{i}.vartype{j}='integer';
            end
        end
    end
end

end % function
