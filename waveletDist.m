function [WEDD]=waveletDist(x,y,level)

% wtOr=tempLoad.wt.(dataName).(leadName{2});
% wtRe=tempLoad.wt.(dataName).(leadName{1});

wtOr=wltTfm(x,level);
wtRe=wltTfm(y,level);

decOr=wtOr.DEC;
acOr=wtOr.ac;
dcOr=wtOr.dc;

% decRe=wtRe.DEC;
acRe=wtRe.ac;
dcRe=wtRe.dc;

energy=zeros(1,length(acOr)+length(dcOr));
wprd=zeros(1,length(acOr)+length(dcOr));
level=length(dcOr);
totEnergy=0;
WEDD=0;

for i=1:1:length(decOr)
    totEnergy=totEnergy+(decOr(i)).^2;  %Total energy of the original DEC
end

for i=1:1:length(acOr{1,1})
    energy(1,1)=energy(1,1)+(acOr{1,1}(i)).^2;  %energy of the approximation coefficient
end

for l=level:-1:1
    for i=1:1:length(dcOr{1,l});
        energy(1,level-l+2)=energy(1,level-l+2)+(dcOr{1,l}(i)).^2;   %energy of the detail coefficient
    end
end

weight=energy/totEnergy; %calculating weights....acL, dcL,...dc1

wprd(1)=sqrt(sum((acOr{1,1}-acRe{1,1}).^2)/sum(acOr{1,1}.^2)); %wprd of the approximation coefficient acL

for i=length(dcOr):-1:1
    wprd(length(dcOr)-i+2)=sqrt(sum((dcOr{1,i}-dcRe{1,i}).^2)/sum(dcOr{1,i}.^2));  %wprd of the detail coefficient dcL...dc1
end


% wprd(1)=sqrt(sum((abs(acOr{1,1})-abs(acRe{1,1})).^2)/sum(abs(acOr{1,1}).^2)); %wprd of the approximation coefficient acL
% 
% for i=length(dcOr):-1:1
%     wprd(length(dcOr)-i+2)=sqrt(sum((abs(dcOr{1,i})-abs(dcRe{1,i})).^2)/sum(abs(dcOr{1,i}).^2));  %wprd of the detail coefficient dcL...dc1
% end


% wprd(1)=sum(acOr{1,1}-acRe{1,1})/sum(acOr{1,1}); %wprd of the approximation coefficient acL
% 
% for i=length(dcOr):-1:1
%     wprd(length(dcOr)-i+2)=sum(dcOr{1,i}-dcRe{1,i})/sum(dcOr{1,i});  %wprd of the detail coefficient dcL...dc1
% end

for i=1:1:length(acOr)+length(dcOr)
    WEDD=WEDD+(weight(i)*wprd(i));
end

% WEDD*100