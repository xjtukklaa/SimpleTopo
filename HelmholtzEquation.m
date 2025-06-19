%% 主要实现Helmholtz方程的矩阵和源项的计算
%% 此处的方程主要实现PDE过滤功能，只需要改变右端项
%% EvaluteFlag == 1表示只计算向量，EvaluteFlag == 2表示同时计算矩阵和向量
function HelmholtzStruct=HelmholtzEquation(CombineValue,Mesh,RhsValue,EvaluteFlag,HelmholtzStruct)
    % 计算标准单元上的形函数的梯度
    GradEvaluteFlag = 2;
    [GradValue,ShapeValue] = HGrad_Q4_Dim2(GradEvaluteFlag);
    % 创建稀疏矩阵
    MatrixMapI = zeros(Mesh.NumCells * 4 * 4, 1);
    MatrixMapJ = zeros(Mesh.NumCells * 4 * 4, 1);
    MatrixValueIJ = zeros(Mesh.NumCells * 4 * 4, 1);
    % 创建右手项
    HelmholtzRHS = zeros(Mesh.NumDofs,1);
    % 循环所有单元,计算对应的Jacobian
    CellDofsMap = zeros(1,4);
    for CellIter = 1:Mesh.NumCells
        % 拿到单元刚度矩阵
        CellDofsMap(:) = Mesh.DofsMap(CellIter,:);
        % 拿到网格的四个点
        CellPoints = Mesh.PointsMap(CellDofsMap,:);
        % 计算Jacobian和Det
        [CellJacobInv,CellJacoDet]=JacobianWeight(CellPoints);
        if EvaluteFlag == 2
            % 设定Map
            MatrixMapTmp = kron(ones(4,1),CellDofsMap);
            MatrixMapI(CellIter*16-15:CellIter*16,1) = MatrixMapTmp(:);
            MatrixMapTmp = kron(ones(1,4),CellDofsMap(:));
            MatrixMapJ(CellIter*16-15:CellIter*16,1) = MatrixMapTmp(:);
            % 设定当前单元的值
            CellMatrix = zeros(4,4);
            GradValueI = zeros(2,1);
            GradValueJ = zeros(2,1);
        end
        CellRhs = zeros(4,1);
        JacobianValue = zeros(2,2);
        % 循环所有正交点
        for QIndex = 1:4
            JacobianValue(:,:) = CellJacobInv(QIndex,:,:);
            % 循环所有自由度
            for DofI = 1:4
                GradValueI(:,1) = GradValue(DofI,QIndex,:);
                if EvaluteFlag == 2
                    for DofJ = 1:4
                    GradValueJ(:,1) = GradValue(DofJ,QIndex,:);
                    CellMatrix(DofI,DofJ) = CellMatrix(DofI,DofJ) + ...
                                            GradValueI' * ...
                                            JacobianValue' * ...
                                            CombineValue(CellIter) * ...
                                            JacobianValue * ...
                                            GradValueJ * ...
                                            CellJacoDet(QIndex) + ...
                                            ShapeValue(DofI,QIndex) * ...
                                            ShapeValue(DofJ,QIndex) * ...
                                            CellJacoDet(QIndex);
                    end
                end
                CellRhs(DofI,1) = CellRhs(DofI,1) + ...
                                  ShapeValue(DofI,QIndex) * ...
                                  RhsValue(CellIter) * ...
                                  CellJacoDet(QIndex);
            end
        end
        % 输入到矩阵和向量中
        if EvaluteFlag == 2
            MatrixValueIJ(CellIter*16-15:CellIter*16,1) = CellMatrix(:);
        end
        HelmholtzRHS(CellDofsMap,1) = HelmholtzRHS(CellDofsMap,1) + CellRhs;
    end
    % 创建稀疏矩阵
    if EvaluteFlag == 2
        HelmholtzStruct.Matrix = ...
         sparse(MatrixMapI,MatrixMapJ,MatrixValueIJ,Mesh.NumDofs,Mesh.NumDofs);
    end
    % 边界处理
    HelmholtzStruct.RHS = HelmholtzRHS;
    % 矩阵重排
    if EvaluteFlag == 2
        HelmholtzStruct.RCMMap = symrcm(HelmholtzStruct.Matrix);
        [~,HelmholtzStruct.IRCMMap] = sort(HelmholtzStruct.RCMMap);
        HelmholtzStruct.Matrix = HelmholtzStruct.Matrix(HelmholtzStruct.RCMMap,HelmholtzStruct.RCMMap);
    end
    HelmholtzStruct.RHS = HelmholtzStruct.RHS(HelmholtzStruct.RCMMap);
    % 计算完成
end