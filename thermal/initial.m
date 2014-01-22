function [Cdt_L, Cdt_V, row, column, value, var] = initial(chip, pack, Layer, K_Material)
%initialization
    var = chip.Nx*chip.Ny*(Layer.N-1)+pack.Nx*pack.Ny*2;
    row = zeros(var*7,1);
    column = zeros(var*7,1);
    value = zeros(var*7,1);

%conductivity arrays
    Cdt_L = zeros(var-pack.Nx*pack.Ny,1);
    Cdt_V = zeros(var-pack.Nx*pack.Ny,1);
    
    Base_material = Layer.material;

%initial the conductivity arrays in chip domain
    for layer_id = 1:1:Layer.N-1
        left_id = 1+chip.Nx*chip.Ny*(layer_id-1);
        right_id = chip.Nx*chip.Ny*layer_id;
        if Base_material(layer_id) == 9
            Cdt_L(left_id:right_id, 1) = K_Material(Base_material(layer_id)+1);
        else
            Cdt_L(left_id:right_id, 1) = K_Material(Base_material(layer_id));
        end
        Cdt_V(left_id:right_id, 1) = K_Material(Base_material(layer_id) );
    end
    
%initial the conductivity arrays in package arrays
    for layer_id = Layer.N:1:Layer.N
        left_id = 1+chip.Nx*chip.Ny*(layer_id-1);
        right_id = pack.Nx*pack.Ny+left_id;
        Cdt_L(left_id:right_id, 1) = K_Material(Base_material(layer_id) );
        Cdt_V(left_id:right_id, 1) = K_Material(Base_material(layer_id) );
    end
end

