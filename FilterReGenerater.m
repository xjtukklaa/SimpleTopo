function OutputVector=FilterReGenerater(Mesh,InputVector)
    CellDofsMap = zeros(1,4);
    CellSolution = zeros(1,4);
    [NumX,NumY] = size(InputVector);
    if NumX > 1
        OutputVector = zeros(Mesh.NumCells,1);
    else
        OutputVector = zeros(1,Mesh.NumCells);
    end
    for CellIter = 1:Mesh.NumCells
        % 拿到单元刚度矩阵
        CellDofsMap(1,:) = Mesh.DofsMap(CellIter,:);
        % 拿到当前单元的解向量
        CellSolution(1,:) = InputVector(CellDofsMap);
        % 计算均值返回向量
        OutputVector(CellIter) = mean(CellSolution);
    end
end