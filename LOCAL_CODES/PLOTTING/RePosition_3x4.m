
% Re-position 3x4 plot to eliminate dead-space
%   [left bottom width height]
hgt = 0.25;
wdth = 0.22;

r1b = 0.7;
r2b = r1b - hgt;
r3b = r1b - 2*hgt;

c1l = 0.1;
c2l = c1l + wdth;
c3l = c1l + 2*wdth;
c4l = c1l + 3*wdth;

% FIRST ROW
    subplot(3,4,1)
        set(gca, 'Position', [c1l r1b wdth hgt])
    subplot(3,4,2)
        set(gca, 'Position', [c2l r1b wdth hgt])
    
% SECOND ROW
    subplot(3,4,5)
        set(gca, 'Position', [c1l r2b wdth hgt])
    subplot(3,4,6)
        set(gca, 'Position', [c2l r2b wdth hgt])

% THIRD ROW
    subplot(3,4,9)
        set(gca, 'Position', [c1l r3b wdth hgt])
    subplot(3,4,10)
        set(gca, 'Position', [c2l r3b wdth hgt])
              
    subplot(3,4,4)
        set(gca, 'Position', [c4l r1b wdth hgt])    
    subplot(3,4,3)
        set(gca, 'Position', [c3l r1b wdth hgt])
    
    subplot(3,4,8)
        set(gca, 'Position', [c4l r2b wdth hgt])
    subplot(3,4,7)
        set(gca, 'Position', [c3l r2b wdth hgt])
    
    subplot(3,4,12)
        set(gca, 'Position', [c4l r3b wdth hgt])
    subplot(3,4,11)
        set(gca, 'Position', [c3l r3b wdth hgt])
    
