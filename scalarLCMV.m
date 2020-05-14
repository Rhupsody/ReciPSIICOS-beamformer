function Z = scalarLCMV(G2dU, Ca)
    [~, Nsites]= size(G2dU);
    Nsites = Nsites/2;
    iCa = tihinv(Ca, 0.01);
    range2d = 1:2;
    Z = zeros(1, Nsites);
    for i = 1:Nsites
        g = G2dU(:,range2d);
        m = g' * iCa * g;
        [~, ss, ~] = eig(m);
        if ss(1, 1) <= ss(2, 2)
            Z(i) = 1 / m(1, 1);
        else
            Z(i) = 1 / m(2, 2);
        range2d = range2d + 2;
    end
end