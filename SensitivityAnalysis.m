%% 计算经典的U\frac{\partial K}{\partial \rho}U
function SensitivityVector=SensitivityAnalysis(SolutionVector,Simp_Rho_Filtered,Mesh,penality,KappaDelta)
    SensitivityVector = zeros(Mesh.NumCells,1);
    Simp_Rho_Sensitivity = penality * KappaDelta * Simp_Rho_Filtered.^ (penality - 1);
    % 计算标准单元上的形函数的梯度
    GradEvaluteFlag = 2;
    [GradValue,~] = HGrad_Q4_Dim2(GradEvaluteFlag);
    % 组合梯度获得当前单元上的解向量梯度
    % 循环所有单元,计算对应的Jacobian
    CellDofsMap = zeros(1,4);
    CellSolution = zeros(1,4);
    for CellIter = 1:Mesh.NumCells
        % 拿到单元自由度MAP
        CellDofsMap(1,:) = Mesh.DofsMap(CellIter,:);
        % 获取解向量
        CellSolution(1,:) = SolutionVector(CellDofsMap);
        % 拿到网格的四个点
        CellPoints = Mesh.PointsMap(CellDofsMap,:);
        % 计算Jacobian和Det
        [CellJacobInv,CellJacoDet]=JacobianWeight(CellPoints);
        % 获取解向量的梯度
        CellSolutionGrad = zeros(4,2);
        CellGradTmp = zeros(1,2);
        JacobianValue = zeros(2,2);
        for CellQPointIndex = 1:4
            JacobianValue(:,:) = CellJacobInv(CellQPointIndex,:,:);
            for CellDofIndex = 1:4
                CellGradTmp(1,:) = GradValue(CellDofIndex,CellQPointIndex,:);
                CellSolutionGrad(CellQPointIndex,:) = CellSolutionGrad(CellQPointIndex,:) + ...
                    CellGradTmp * JacobianValue * CellSolution(CellDofIndex) * CellJacoDet(CellQPointIndex);
            end
        end
        % 计算U\frac{\partial K}{\partial \rho}U
        SensitivityVector(CellIter) = -sum(CellSolutionGrad * ...
                                      Simp_Rho_Sensitivity(CellIter) * ...
                                      CellSolutionGrad',"all");
    end
end