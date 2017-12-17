function FAIRerror(err)

if nargin == 0
    FAIRerror('nice error');
    return;
end

fprintf(2,'\n\n ----- FAIRerror.m \n');

if ischar(err)
    fprintf(2,'%s',err);
elseif isa(err,'MException')
    fprintf(2,'%s - %s',err.identifier,err.message);
else
    fprintf(2,'No clue, sry.');
end

fprintf(2,'\n ----- \n\n');