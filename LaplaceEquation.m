%% 主要实现Laplace方程的矩阵和源项的计算
function LaplaceStruct = LaplaceEquation(CombineValue,Mesh,Boundary)
    % 计算标准单元上的形函数的梯度
    GradEvaluteFlag = 2;
    [GradValue,ShapeValue] = HGrad_Q4_Dim2(GradEvaluteFlag);
    % 创建稀疏矩阵
    MatrixMapI = zeros(Mesh.NumCells * 4 * 4, 1);
    MatrixMapJ = zeros(Mesh.NumCells * 4 * 4, 1);
    MatrixValueIJ = zeros(Mesh.NumCells * 4 * 4, 1);
    % 创建右手项
    LaplaceRHS = zeros(Mesh.NumDofs,1);
    % 创建固定边界约束的计算向量
    DirichletRHS = zeros(Mesh.NumDofs,1);
    DirichletDofs  = Boundary("Dirichlet");
    DirichletDofs = cell2mat(DirichletDofs);
    DirichletValue = Boundary("Dirichlet_Value");
    DirichletValue = cell2mat(DirichletValue);
    DirichletRHS(DirichletDofs,1) = DirichletValue;
    % 循环所有单元,计算对应的Jacobian
    CellDofsMap = zeros(1,4);
    for CellIter = 1:Mesh.NumCells
        % 拿到单元自由度MAP
        CellDofsMap(:) = Mesh.DofsMap(CellIter,:);
        % 拿到网格的四个点
        CellPoints = Mesh.PointsMap(CellDofsMap,:);
        % 计算Jacobian和Det
        [CellJacobInv,CellJacoDet]=JacobianWeight(CellPoints);
        % 设定Map
        MatrixMapTmp = kron(ones(4,1),CellDofsMap);
        MatrixMapI(CellIter*16-15:CellIter*16,1) = MatrixMapTmp(:);
        MatrixMapTmp = kron(ones(1,4),CellDofsMap(:));
        MatrixMapJ(CellIter*16-15:CellIter*16,1) = MatrixMapTmp(:);
        % 设定当前单元的值
        CellMatrix = zeros(4,4);
        CellRhs = zeros(4,1);
        CellQPoitns = zeros(4,2);
        CellQPoitns(:,:) = Mesh.CellQPoints(CellIter,:,:);
        GradValueI = zeros(2,1);
        GradValueJ = zeros(2,1);
        JacobianValue = zeros(2,2);
        % 循环所有正交点
        for QIndex = 1:4
            JacobianValue(:,:) = CellJacobInv(QIndex,:,:);
            % 循环所有自由度
            for DofI = 1:4
                GradValueI(:,1) = GradValue(DofI,QIndex,:);
                for DofJ = 1:4
                GradValueJ(:,1) = GradValue(DofJ,QIndex,:);
                CellMatrix(DofI,DofJ) = CellMatrix(DofI,DofJ) + ...
                                        GradValueI' * ...
                                        JacobianValue' * ...
                                        CombineValue(CellIter) * ...
                                        JacobianValue * ...
                                        GradValueJ * ...
                                        CellJacoDet(QIndex);
                end
                CellRhs(DofI,1) = CellRhs(DofI,1) + ...
                                  ShapeValue(DofI,QIndex) * ...
                                  SourceLaplace(CellQPoitns(QIndex,:)) * ...
                                  CellJacoDet(QIndex);
            end
        end
        % 输入到矩阵和向量中
        MatrixValueIJ(CellIter*16-15:CellIter*16,1) = CellMatrix(:);
        LaplaceRHS(CellDofsMap,1) = LaplaceRHS(CellDofsMap,1) + CellRhs;
    end
    % 创建稀疏矩阵
    LaplaceStruct.Matrix = ...
     sparse(MatrixMapI,MatrixMapJ,MatrixValueIJ,Mesh.NumDofs,Mesh.NumDofs);
    % 边界处理
    LaplaceRHS = LaplaceRHS - LaplaceStruct.Matrix * DirichletRHS;
    LaplaceRHS(DirichletDofs,1) = DirichletValue;
    LaplaceStruct.RHS = LaplaceRHS;
    % 重置自由度
    for DDofs = 1:length(DirichletDofs)
        LaplaceStruct.Matrix(DirichletDofs(DDofs),:) = 0;
        LaplaceStruct.Matrix(:,DirichletDofs(DDofs)) = 0;
        LaplaceStruct.Matrix(DirichletDofs(DDofs),DirichletDofs(DDofs)) = 1;
    end
    % 矩阵重排
    LaplaceStruct.RCMMap = symrcm(LaplaceStruct.Matrix);
    [~,LaplaceStruct.IRCMMap] = sort(LaplaceStruct.RCMMap);
    LaplaceStruct.Matrix = LaplaceStruct.Matrix(LaplaceStruct.RCMMap,LaplaceStruct.RCMMap);
    LaplaceStruct.RHS = LaplaceStruct.RHS(LaplaceStruct.RCMMap);
    % 计算完成
end