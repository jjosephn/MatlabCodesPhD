function [x]=invWltTfm(wltCoeff,len,level)

% wltCoeff=w1;

wname = 'bior4.4'; % biorthogonal 9/7 filter

[Rf,Df] = biorwavf(wname);

[~,~,Lo_R,Hi_R] = biorfilt(Df,Rf);

wltCoeff.L(1)=length(wltCoeff.ac{1});

for i=1:1:level
    wltCoeff.L(i+1)=length(wltCoeff.dc{level-i+1});
end

wltCoeff.L(level+2)=len;

wltCoeff.DEC(1:wltCoeff.L(1))=wltCoeff.ac{1};

k=wltCoeff.L(1)+1;
j=2;
for i=level:-1:1
    wltCoeff.DEC(k:k+wltCoeff.L(j)-1)=wltCoeff.dc{i};
    k=k+wltCoeff.L(j);
    j=j+1;
end

x = waverec(wltCoeff.DEC,wltCoeff.L,Lo_R,Hi_R);


% [Ea,Ed] = wenergy(wltCoeff.DEC,wltCoeff.L);

% wltCoeff.ac{1}=wltCoeff.DEC(1:wltCoeff.L(1));
% 
% k=wltCoeff.L(1)+1;
% j=2;
% for i=level:-1:1
%     wltCoeff.dc{i}=wltCoeff.DEC(k:k+wltCoeff.L(j)-1);
%     k=k+wltCoeff.L(j);
%     j=j+1;
% end

% err1 = norm(x1-x)
% plot(x);

