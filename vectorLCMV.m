function Z = vectorLCMV(G2dU_red, Ca)
    [Nch, Nsites_red]= size(G2dU_red);
    Nsites_red = Nsites_red/2;
    iCa = tihinv(Ca, 0.01);
    range2d = 1:2;
    Z = zeros(1, Nsites_red);
    for i = 1:Nsites_red
        g = G2dU_red(:,range2d);
        m = inv(g'*iCa*g);
        [u ss v] = svd(m);
        Z(i) = ss(1,1);
        range2d = range2d+2;
    end
end