function Layer = stack_build(die, thick)
    %for each die, it has three layers: die & ILD & Under
    % TIM + interposer;
    Layer.N = (die.N*die.model) + 2;
    disp(['there are ' num2str(Layer.N) ' layers']);
    Layer.material = zeros(Layer.N,1);
    Layer.thick = zeros(Layer.N,1);
    
    Layer.material(1) = 1; %TIM;
    Layer.material(Layer.N) = 7; %interposer;
    
    for i=1:1:die.N     
        for j=1:1:die.model
            id = 1+(i-1)*die.model+j;
            if j==1
                Layer.material(id) = 2;
            elseif j == 2
                Layer.material(id) = 9;
            else
                Layer.material(id) = 3;
            end
        end
    end
    
    Layer.thick(1) = thick.tim; %TIM;
    Layer.thick(Layer.N) = thick.inter; %interposer;
    
    for i=1:1:die.N     
        for j=1:1:die.model
            id = 1+(i-1)*die.model+j;
            if j==1
                Layer.thick(id) = thick.die;
            elseif j == 2
                Layer.thick(id) = thick.ild;
            else
                Layer.thick(id) = thick.under;
            end
        end
    end
    Layer.thick(Layer.N-1) = thick.bump;
end

