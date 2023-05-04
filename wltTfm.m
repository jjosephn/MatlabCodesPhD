function [wltCoeff]=wltTfm(x,level)

wname = 'bior4.4'; % biorthogonal 9/7 filter

[Rf,Df] = biorwavf(wname);

[Lo_D,Hi_D,~,~] = biorfilt(Df,Rf);

% wname = 'haar';
% 
% [Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(wname);



[wltCoeff.DEC,wltCoeff.L] = wavedec(x,level,Lo_D,Hi_D); %decomposition coefficient vector

% [Ea,Ed] = wenergy(wltCoeff.DEC,wltCoeff.L);

wltCoeff.ac{1}=wltCoeff.DEC(1:wltCoeff.L(1));

k=wltCoeff.L(1)+1;
j=2;
for i=level:-1:1
    wltCoeff.dc{i}=wltCoeff.DEC(k:k+wltCoeff.L(j)-1);
    k=k+wltCoeff.L(j);
    j=j+1;
end


