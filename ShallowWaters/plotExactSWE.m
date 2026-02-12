%% Plot exact solution, left part in blue, right part in yellow, as in the experiment 
fig1=figure(1);
% Compute solution with the exact Riemann solver on a finer mesh 
EMAX = 10*IMAX; 
xe = linspace(xL,xR,EMAX); 
Qe = zeros(3,EMAX);
for i= 1:EMAX
    xi = xe(i)/time;
    Qe(:,i) = ExactRiemannSW(QL,QR,xi);
end
psie   = Qe(3,:)./Qe(1,:); 
ue     = Qe(2,:)./Qe(1,:); 
lefte  = psie>0.5;
righte = psie<=0.5;
nl     = 0; 
for i=1:EMAX
    if(lefte(i)==1)
        nl = nl + 1; 
    end
end
subplot(3,1,1)
hold on
plot(xe(lefte), Qe(1,lefte),'b-')              
plot(xe(righte),Qe(1,righte),Color=[0.9290 0.6940 0.1250],LineStyle='-')
subplot(3,1,2)
hold on
plot(xe(lefte), ue(1,lefte),'b-')              
plot(xe(righte),ue(1,righte),Color=[0.9290 0.6940 0.1250],LineStyle='-')
subplot(3,1,3)
hold on
plot([xe(lefte),xe(nl)], [psie(1,lefte),0.5],'b-')              
plot([xe(nl),xe(righte)],[0.5,psie(1,righte)],Color=[0.9290 0.6940 0.1250],LineStyle='-')
drawnow