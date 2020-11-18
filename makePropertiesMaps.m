function propertiesMaps = makePropertiesMaps(propMapsDefs)
mykeys = keys(propMapsDefs);

numkeys = numel(mykeys);
v = zeros(1,numkeys);
ind = 1;
for key = mykeys
    n = numel(propMapsDefs(key{1}));
    v(ind) = n;
    ind = ind + 1;
end
propertiesMaps = cell(v);
indArray = zeros(v);
for i=1:prod(v)
    indArray(i)=i;
end
i = zeros(1,numkeys);
str = '[';
for j = 1:n
    str = [str ' i(' num2str(j) ')'];
end
str = [str ']'];
for j = 1:prod(v)
    eval([str ' = ind2sub(v,' num2str(j) ')']);
    myValues = cell(1,numkeys);
    str2 = 'propertiesMaps{';
    for k = 1:numkeys
        key = mykeys{k};
        key_values = propMapsDefs(key);
        myValues{k} = key_values{i(k)};
        str2 = [str2 num2str(i(k)) ','];
    end
    str2 = [str2(1:end-1) '}'];
    props = containers.Map(mykeys,myValues);
    eval([str2 '= props']);
end
        
    

    
