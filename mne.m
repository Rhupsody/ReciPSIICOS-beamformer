function Zmne = mne(Wmne, Ca)
        Zmne = zeros(1, size(Wmne, 1)/2);
        range2d = 1:2;
        for i = 1:(size(Wmne, 1)/2)
            w = Wmne(range2d,:);
            m = w*Ca*w';
            [~, ss, ~] = svd(m);
            Zmne(i) = ss(1,1);
            range2d = range2d+2;
        end
end