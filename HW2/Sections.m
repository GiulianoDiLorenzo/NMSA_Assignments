% vocal tract profile, non-dimensional [pos S] pairs
% /E/
S1 = [0 1;0.09 0.4;0.11 2.4;0.24 2.4;0.26 3.2;0.29 3.2;0.32 4.2;...
0.41 4.2;0.47 3.2;0.59 1.8;0.65 1.6;0.71 1.6;0.74 1;0.76 0.8;...
0.82 0.8;0.88 2;0.91 2;0.94 3.2;1 3.2];
% /A/
S2 = [0 1;0.03 0.60;0.09 0.4;0.12 1.6;0.18 0.6;0.29 0.2;0.35 0.4;...
    0.41 0.8;0.47 1;0.50 0.6;0.59 2;0.65 3.2;0.85 3.2;0.94 2;1 2];

figure(100)
subplot(2,1,1); plot(S2(:,1), sqrt(S2(:,2)),'k', S2(:,1), -sqrt(S2(:,2)),'k')
title('Vocal Tract Profile'); xlabel('x'); ylabel('sqrt(S)');
subplot(2,1,2); plot(S1(:,1), sqrt(S1(:,2)),'k', S1(:,1), -sqrt(S1(:,2)),'k')
title('Vocal Tract Profile'); xlabel('x'); ylabel('sqrt(S)');

