mech = 'sanDiego';

fnm = strcat('csp.',mech,'.1200K.temp.bin');
fid = fopen(fnm);
out = fread(fid,[11,558],'double');
fclose(fid);
out  = out';
t    = out(:,1);
y    = out(:,2:10);
T    = out(:,11);

fnm = strcat('csp.',mech,'.1200K.nslow.bin');
fid = fopen(fnm);
out = fread(fid,[1,558],'int');
fclose(fid);
out   = out';
nslow = out(:);

root = '~/Documents/UIUC/thesis/litreview/notes/csp/data/';
fnm = strcat(root,mech,'.1200K.temp.dat');
fid = fopen(fnm,'w');
fprintf(fid,'%s %4s\n','t','temp');
fprintf(fid,'%6.4f %6.4f\n',[log10(t)'; T']);
edit(fnm)

root = '~/Documents/UIUC/thesis/litreview/notes/csp/data/';
fnm = strcat(root,mech,'.1200K.nslow.dat');
fid = fopen(fnm,'w');
fprintf(fid,'%s %5s\n','t','nslow');
fprintf(fid,'%6.4f %6.4f\n',[log10(t)'; nslow']);
edit(fnm)

