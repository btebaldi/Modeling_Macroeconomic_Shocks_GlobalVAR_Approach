function varvec=convertstr2vec(vcell)

%**************************************************************************
% PURPOSE: Converts a comma delimited string to a vector of strings. For
% example the string 'r, lr' is converted to a 2x1 string vector with 'r'
% as the first element and 'lr' as the second element
%**************************************************************************     
     
varvec = {};
acc = [];
for i=1:length(vcell)
    
    letter = {vcell(i)};
    if not(strcmp(letter,' '));
        if strcmp(letter,',')
            acc = {acc};
            varvec = [varvec; acc]; %#ok
            acc = [];
        else
            acc = [acc vcell(i)]; %#ok
        end
    end
    
    if i == length(vcell)
        acc = {acc};
        varvec = [varvec; acc]; %#ok
    end
end