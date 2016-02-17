function [row, column, value, var, var_chip] = initial(chip, TSV, type)
    %in our modeling, each mesh node is an unknown;
    %for chip domain, there are chip.Nx*chip.Ny*die_num unknowns.
    var_chip = chip.N * chip.Nx * chip.Ny;    
    
    %we need to add #PAD unknowns for the current of flowing into each pad
    %we handle the power and ground pads/solution differently, thus we add
    %different number of unknowns
    if type == 1
        var = var_chip + TSV.P*chip.N;
    else
        var = var_chip + TSV.G*chip.N;
    end

    %for each unknown, it may have correlation with 7 nodes,
    %they are: left, right, north, south, bottom, top, itself
    %we use all the three arrays to store the non-zero element in the
    %future built conductance and capacitive/inductive MATRIX
    %keep in mind, conductance and capacitive/inductive MATRIX are highly
    %sparse, because for each unknown, we set it as one row, so the size of
    %such MATRIX is var*var. e.g. if the solution space is 100*100*2;
    %such MATRIX is (100*100*2)^2 large
    row = zeros(var*7, 1);
    column = zeros(var*7, 1);
    value = zeros(var*7, 2);
    
    %there are three matrixes to be calculated
    %capacitance/inductance matrix C;
    %conductive matrix    
end

