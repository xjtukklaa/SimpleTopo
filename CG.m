%% CG求解器，可带预处理器使用
function [SolverFlag, Iter, Last_Value, Vector_Output] = CG(Matrix, RHS, Tol, MaxIter)
    % 设定初始猜测解
    Vector_Output = ones(length(RHS),1);
    % 计算残差向量r=Ax - b
    Epsi = Matrix * Vector_Output - RHS;
    % 初始化共轭梯度方向
    Direct = -Epsi;
    for Index = 1:MaxIter
        if max(abs(Epsi)) < Tol
            SolverFlag = true;
            Iter = Index;
            Last_Value = max(abs(Epsi));
            return;
        else
            % 计算步长
            Gamma = -Epsi' * Direct/(Direct' * Matrix * Direct);
            % 更新解向量
            Vector_Output = Vector_Output + Gamma * Direct;
            % 更新残差
            Epsi = Matrix * Vector_Output - RHS;
            % 更新共轭步长
            Beta = (Epsi' * Matrix * Direct)/(Direct' * Matrix * Direct);
            % 更新共轭方向
            Direct = -Epsi + Beta * Direct;
        end
    end
    SolverFlag = false;
    Iter = Index;
    Last_Value = 0;
    Vector_Output = ones(length(RHS),1);
end