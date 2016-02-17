function display_T(chip, x)
    gridNx_chip = chip.Nx;
    gridNy_chip = chip.Ny;    
    const = gridNx_chip * gridNy_chip;
    drawP_die = zeros(gridNy_chip, gridNx_chip, chip.N);
    for k=1 : 1 : chip.N
        for i = 1:1:gridNx_chip
            for j=1:1:gridNy_chip
                id = i + (j-1)*gridNx_chip + const*(k-1);
                drawP_die(j,i,k) = x(id);
            end
        end
        disp(['Die ', num2str(k), '; Max: ', num2str(max(max(drawP_die(:,:,k))))]);
        disp(['Die ', num2str(k), '; Min: ', num2str(min(min(drawP_die(:,:,k))))]);
    end
end