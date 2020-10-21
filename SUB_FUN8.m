function [ Unow, CNTR, now_obj_fcn ] = SUB_FUN8( INP_IMG, N_CLTR )


if nargin < 2
    N_CLTR = 2;   
end

[row, col] = size(INP_IMG);
E_N = 2;      
EPLN_VL = 0.001;  
MX_ITR = 100;  

UP_VAL = rand(row, col, N_CLTR);
DEP_SUM = sum(UP_VAL, 3);
DEP_SUM = repmat(DEP_SUM, [1,1, N_CLTR]);
UP_VAL = UP_VAL./DEP_SUM;

CNTR = zeros(N_CLTR,1); 

for i=1:N_CLTR
    CNTR(i,1) = sum(sum(UP_VAL(:,:,i).*INP_IMG))/sum(sum(UP_VAL(:,:,i)));
end

pre_obj_fcn = 0;
for i=1:N_CLTR
    pre_obj_fcn = pre_obj_fcn + sum(sum((UP_VAL(:,:,i) .*INP_IMG - CNTR(i)).^2));
end
fprintf('Initial objective fcn = %f\n', pre_obj_fcn);

for ITR = 1:MX_ITR    

    Unow = zeros(size(UP_VAL));
    for i=1:row
        for j=1:col
            for uII = 1:N_CLTR
                tmp = 0;
                for uJJ = 1:N_CLTR
                    disUp = abs(INP_IMG(i,j) - CNTR(uII));
                    disDn = abs(INP_IMG(i,j) - CNTR(uJJ));
                    tmp = tmp + (disUp/disDn).^(2/(E_N-1));
                end
                Unow(i,j, uII) = 1/(tmp);
            end            
        end
    end   

    now_obj_fcn = 0;
    for i=1:N_CLTR
        now_obj_fcn = now_obj_fcn + sum(sum((Unow(:,:,i) .*INP_IMG - CNTR(i)).^2));
    end
    fprintf('Iter = %d, Objective = %f\n', ITR, now_obj_fcn);

    if max(max(max(abs(Unow-UP_VAL))))<EPLN_VL || abs(now_obj_fcn - pre_obj_fcn)<EPLN_VL %引入2个判定
        break;
    else
        UP_VAL = Unow.^E_N;
        for i=1:N_CLTR
            CNTR(i,1) = sum(sum(UP_VAL(:,:,i).*INP_IMG))/sum(sum(UP_VAL(:,:,i)));
        end
        pre_obj_fcn = now_obj_fcn;
    end
end
