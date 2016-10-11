doplots = 0;
mech    = 'sanDiego';

fnm = strcat('csp.',mech,'.1200K.bin');
fid = fopen(fnm);
if(strcmp(mech,'sanDiego'))
    out = fread(fid,[11,558],'double');
elseif(strcmp(mech,'boivin'))
    out = fread(fid,[11,627],'double');
end
fclose(fid);

out  = out';
t = out(:,1);
y = out(:,2:10);
T = out(:,11);

figure(1)
plot(log10(t),T,'r','linewidth',2)
set(gca,'fontsize',18)
xlim([-5 -3.95])
