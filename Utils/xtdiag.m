% draw the x-t diagram of the exact solution 
% Input:
%   - QL, QR : left and right state vectors
%   - xi     : similarity parameter (xi=x/t)
% Output:
%   - Q      : solution Q=Q(xi)
function ierr = xtdiag(QL,QR,xL,xR,tend,time)
    global g fig2       % gravity constant
    hL = QL(1);         % left depth
    hR = QR(1);         % right depth
    uL = QL(2)/QL(1);   % left velocity
    uR = QR(2)/QR(1);   % right velocity
    psiL = QL(3)/QL(1); % left concentration
    psiR = QR(3)/QR(1); % right concentration
    hs = Newton(hL,hR,uL,uR); % compute h* using Newton
    us = uL - phi(hs,hL);     % compute u* by definition

    fig2=figure(2); 
    fig2.Position = [680   108   560   420]; 
    axis([xL xR 0 tend])
    clf  

    % "sampling" of the solution
    line([0,us*tend],[0,tend],'Color','black','LineStyle','--')
    if(hs>hL)
        % left shock
        s = uL - sqrt(0.5*g*hs/hL*(hL+hs));
        line([0,s*tend],[0,tend],'Color','blue','LineStyle','-')
    else
        % left rarefaction
        head = uL - sqrt(g*hL);
        tail = us - sqrt(g*hs);
        nr = 10; 
        ir = linspace(head,tail,nr); 
        for i=1:nr
            line([0,ir(i)*tend],[0,tend],'Color','blue','LineStyle','-');
        end
    end
    if(hs>hR)
        % right shock
        s = uR + sqrt(0.5*g*hs/hR*(hs+hR)); 
        line([0,s*tend],[0,tend],'Color',[0.9290 0.6940 0.1250],'LineStyle','-')
    else
        % right rarefaction
        tail = us + sqrt(g*hs); 
        head = uR + sqrt(g*hR); 
        nr = 10; 
        ir = linspace(head,tail,nr); 
        for i=1:nr
            line([0,ir(i)*tend],[0,tend],'Color',[0.9290 0.6940 0.1250],'LineStyle','-');
        end
    end
    line([xL,xR],[time,time],'Color','red')
    xlabel('x')
    ylabel('t')
    ierr = 0; 

end