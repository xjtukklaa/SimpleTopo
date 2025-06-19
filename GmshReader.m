%% 网格的边界设定，1 -> Dirichlet，2 -> Neumann
function MeshStruct = GmshReader(Dim,MeshType)
    addpath("Mesh\");
    % 首先读取Gmsh文件
    % 第一个是只有实体网格
    MeshWithBoundary = MeshPart;
    % 这个是只有边界
    MeshBoundary = MeshSuf;
    % 设定网格维数
    MeshStruct.Dim = Dim;
    % 获得网格数量
    MeshStruct.NumCells = length(MeshWithBoundary.QUADS);
    % 获得自由度数量
    % 此处为节点个数,为向量场
    MeshStruct.NumDofs = MeshWithBoundary.nbNod;
    % 节点和自由度之间的联系
    if MeshType == "QUAD"
        MeshStruct.DofsMap = MeshWithBoundary.QUADS(:,1:4);
    % 三角形导入还未完全完成
    elseif MeshType == "Tri"
        MeshStruct.DofsMap = MeshWithBoundary.QUADS(:,1:3);
    else
        error("Error Mesh Type");
    end
    % 节点Map
    MeshStruct.PointsMap = MeshWithBoundary.POS(:,1:MeshStruct.Dim);
    % 固定约束的自由度
    DirichletPOS = MeshBoundary.POS;
    DirichletIndex = zeros(length(DirichletPOS),1);
    for Index = 1:length(DirichletPOS)
        DirichletCellPos = DirichletPOS(Index,1:MeshStruct.Dim);
        for PosIndex = 1:MeshWithBoundary.nbNod
            CellPos = MeshWithBoundary.POS(PosIndex,1:MeshStruct.Dim);
            if DirichletCellPos == CellPos
                DirichletIndex(Index) = PosIndex;
            end
        end
    end
    DirichletIndex = unique(DirichletIndex);
    MeshStruct.BoundaryDofsMap = dictionary("Dirichlet",{DirichletIndex});
    clear MeshWithBoundary DirichletPOS DirichletIndex;
    %% 附加信息
    % 当前网格面积
    MeshStruct.CellMeasure = zeros(MeshStruct.NumCells,1);
    CellPointsLoop = zeros(5, MeshStruct.Dim);
    for CellIter = 1:MeshStruct.NumCells
        CellPointsLoop(1:4,:) = MeshStruct.PointsMap(MeshStruct.DofsMap(CellIter,:),:);
        CellPointsLoop(5,:) = CellPointsLoop(1,:);
        MeshStruct.CellMeasure(CellIter) = polyarea(CellPointsLoop(:,1),CellPointsLoop(:,2));
    end
    % 所有单元的正交点坐标
    EvaluteShapeFlag = 1;
    [~,ShapeFunction] = HGrad_Q4_Dim2(EvaluteShapeFlag);
    CellPoints  = zeros(4, MeshStruct.Dim);
    MeshStruct.CellQPoints = zeros(MeshStruct.NumCells, 4, MeshStruct.Dim);
    for CellIter = 1:MeshStruct.NumCells
        CellPoints(:,:) = MeshStruct.PointsMap(MeshStruct.DofsMap(CellIter,:),:);
        MeshStruct.CellQPoints(CellIter,:,:) = ShapeFunction' * CellPoints;
    end
    % 计算完成
end