function mfile=read_mfile(path_mfile)
%% read the EFIT input data from mfile
% -inputs:
% path_mfile    string, the path of the mfile
% -outputs:
% mfile             struct, the data passed from mfile

% Edited by Shuying SUN in 2019/05/27
% Modified by Bin Chen in 2019/05/31

fid=fopen(path_mfile,'r');
if(fid<=0)
    error(['Can not open file:', path_mfile]);
end
A={};
i=0;    % namelist index
j=0;    % variable index
namestr='';
valuestr='';
linenum=0;  % locate current line number
while ~feof(fid)
    linestr=fgetl(fid);
    linenum=linenum+1;
    %% remove comments, save effective string to effline, and trim out the space
    if(~isempty(strfind(linestr,'!')))    % this line contains comments
        tmpstrcell=strsplit(linestr,'!');
        effline=tmpstrcell{1};
    else    % this line does not contains comments
        effline=linestr;
    end
    trimedline=strtrim(effline);
    %% if it's an invalid line (a line only contains comments)
    if(isempty(trimedline))
        continue;
    end
    firstchar=trimedline(1);
    % if the line is started with '&', it means this line is a start line
    % of a new namelist
    if(firstchar=='&' || firstchar=='$')
        i=i+1;
        A{i}.namelist=trimedline(2:end);
        A{i}.varname={};
        A{i}.varvalue={};
        A{i}.vartype={};
        j=0;
    % if the line is started with '/', it means this line is an end line of
    % the last namelist
    elseif(firstchar=='/')
        % handling the last variable in this namelist
        if(~isempty(valuestr))
            j=j+1;
            A{i}.varname{j}=namestr;
            if (~isempty(strfind(valuestr,'t')) || ~isempty(strfind(valuestr,'T')))
                A{i}.vartype{j}='logical';
                A{i}.varvalue{j}=1;
            elseif (~isempty(strfind(valuestr,'f')) || ~isempty(strfind(valuestr,'F')))
                A{i}.vartype{j}='logical';
                A{i}.varvalue{j}=0;
            elseif (~isempty(strfind(valuestr,'.')) || ~isempty(strfind(valuestr,'e')) || ~isempty(strfind(valuestr,'E')))
            % use this rather than 'contains' to make it compatible with octave
                A{i}.vartype{j}='double';
                A{i}.varvalue{j}=str2num(valuestr);
            else
                A{i}.vartype{j}='integer';
                A{i}.varvalue{j}=str2num(valuestr);
            end
            namestr='';
            valuestr='';
        end
    % if the line is started with English character, it means this line is
    % a start line of a new variable
    elseif( (abs(firstchar)>=65 && abs(firstchar) <=90) || ...
            (abs(firstchar)>=97 && abs(firstchar) <=122))
        % handling the previous variable in this namelist
        if(~isempty(valuestr))
            j=j+1;
            A{i}.varname{j}=namestr;
            if (~isempty(strfind(valuestr,'t')) || ~isempty(strfind(valuestr,'T')))
                A{i}.vartype{j}='logical';
                A{i}.varvalue{j}=1;
            elseif (~isempty(strfind(valuestr,'f')) || ~isempty(strfind(valuestr,'F')))
                A{i}.vartype{j}='logical';
                A{i}.varvalue{j}=0;
            elseif (~isempty(strfind(valuestr,'.')) || ~isempty(strfind(valuestr,'e')) || ~isempty(strfind(valuestr,'E')))
            % use this rather than 'contains' to make it compatible with octave
                A{i}.vartype{j}='double';
                A{i}.varvalue{j}=str2num(valuestr);
            else
                A{i}.vartype{j}='integer';
                A{i}.varvalue{j}=str2num(valuestr);
            end
            namestr='';
            valuestr='';
        end
        if(~isempty(strfind(linestr,'=')))    % this line contains '='
            newstr=strsplit(trimedline,'=');
        else
            continue;
        end
        namestr=strtrim(newstr{1});
        valuestr=strtrim(newstr{2});
    % if the line is started with number character, it means this line is a
    % continued line of the current variable
    elseif( firstchar=='+' || firstchar=='-' || firstchar=='.' || ...
            (abs(firstchar)>=48 && abs(firstchar)<=57))
        valuestr=strcat(valuestr,[' ',trimedline]);
    % if the line is started with '!', it means this line is a comment line
    elseif(firstchar=='!')
        % do nothing, no handling the comments
    % other situations, need to update
    else
        error(['Undefined handling method in read_mfile for mfile ',path_mfile,...
            10, 'at line ',num2str(linenum)]);
    end
end % while
fclose(fid);

for i=1:length(A)
    eval(['mfile.',A{i}.namelist,' = A{i};']);
end

end % function
