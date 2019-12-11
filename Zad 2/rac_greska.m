function greska = rac_greska(S, D)
    ID = [];
    greska = 0;
    for i = 1:length(S)
        tren = 1e18;
        id = 0;
        for j = 1:length(D)
            nrm = norm(S(i) - D(j));
            if tren > nrm && ~any(ID == j)
                id = j;
                tren = nrm;
            end
        end
        ID(i) = id;
        greska = greska + tren;
    end
    greska = greska / norm(D);
end