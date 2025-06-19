%% 创建长方形和正方形的均匀网格
%% NumX,NumY为个数,LengthX,LengthY为实际长度
function MeshStruct = MeshCreater(NumX,NumY,LengthX,LengthY)
    %% 基本信息
    % X方向上的网格长度
    MeshStruct.SizeX = LengthX/NumX;
    % X方向上的网格长度
    MeshStruct.SizeY = LengthY/NumY;
    % 
    MeshStruct.NumX = NumX;
    MeshStruct.NumY = NumY;
    %
    MeshStruct.MaxX = LengthX;
    MeshStruct.MaxY = LengthY;
    % 最大的单元个数
    MeshStruct.NumCells = NumX * NumY;
    % 单元网格的全局对照表
    CellMatrixAll = 1:MeshStruct.NumCells;
    CellMatrixAll = reshape(CellMatrixAll,MeshStruct.NumY,MeshStruct.NumX);
    % 最大的自由度个数
    MeshStruct.NumDofs = (NumX + 1) * (NumY + 1);
    % 单元的局部自由度和全局自由度对照表
    DofsMatrixAll = 1:MeshStruct.NumDofs;
    DofsMatrixAll = reshape(DofsMatrixAll,MeshStruct.NumY + 1,MeshStruct.NumX + 1);
    % 第四个自由度的全局自由度
    DofsMatrixEnd = DofsMatrixAll(1:end-1,1:end-1);
    % 自由度对照表
    MeshStruct.DofsMap = kron(ones(1,4),DofsMatrixEnd(:)) + ...
                        [1,MeshStruct.NumY+2,MeshStruct.NumY+1,0];
    clear DofsMatrixEnd;
    % 网格维度
    MeshStruct.Dim = 2;
    % 每个自由度对应的点坐标,默认另外一个点为(0,0)
    PointsX = 0:MeshStruct.SizeX:MeshStruct.MaxX;
    PointsXAxis = kron(ones(MeshStruct.NumY + 1,1),PointsX(:)');

    PointsY = MeshStruct.MaxY:-MeshStruct.SizeY:0;
    PointsYAxis = kron(ones(1,MeshStruct.NumX + 1),PointsY(:));

    MeshStruct.PointsMap = [PointsXAxis(:)';PointsYAxis(:)']';
    clear PointsX PointsXAxis PointsY PointsYAxis;
    %% 附加信息
    % 当前网格面积
    MeshStruct.CellMeasure = zeros(MeshStruct.NumCells,1);
    CellPointsLoop = zeros(5, MeshStruct.Dim);
    for CellIter = 1:MeshStruct.NumCells
        CellPointsLoop(1:4,:) = MeshStruct.PointsMap(MeshStruct.DofsMap(CellIter,:),:);
        CellPointsLoop(5,:) = CellPointsLoop(1,:);
        MeshStruct.CellMeasure(CellIter) = polyarea(CellPointsLoop(:,1),CellPointsLoop(:,2));
    end
    % 边界网格映射表
    BoudaryKeys = ["Upp","Down","Left","Right"];
    UppCells   = CellMatrixAll(1,:);
    DownCells  = CellMatrixAll(end,:);
    LeftCells  = CellMatrixAll(:,1);
    RightCells = CellMatrixAll(:,end);
    MeshStruct.BoundaryCellMap = dictionary(BoudaryKeys,{UppCells(:)',DownCells(:)',LeftCells(:)',RightCells(:)'});
    clear UppCells DownCells RightCells LeftCells;
    % 边界自由度映射表
    UppDofs   = DofsMatrixAll(1,:);
    DownDofs  = DofsMatrixAll(end,:);
    LeftDofs  = DofsMatrixAll(:,1);
    RightDofs = DofsMatrixAll(:,end);
    MeshStruct.BoundaryDofsMap = dictionary(BoudaryKeys,{UppDofs(:)',DownDofs(:)',LeftDofs(:)',RightDofs(:)'});
    clear UppDofs DownDofs RightDofs LeftDofs;
    % 所有单元的正交点坐标
    EvaluteShapeFlag = 1;
    [~,ShapeFunction] = HGrad_Q4_Dim2(EvaluteShapeFlag);
    CellPoints  = zeros(4, MeshStruct.Dim);
    MeshStruct.CellQPoints = zeros(MeshStruct.NumCells, 4, MeshStruct.Dim);
    for CellIter = 1:MeshStruct.NumCells
        CellPoints(:,:) = MeshStruct.PointsMap(MeshStruct.DofsMap(CellIter,:),:);
        MeshStruct.CellQPoints(CellIter,:,:) = ShapeFunction' * CellPoints;
    end
    % 创建完成
end