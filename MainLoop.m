clear;
clc;
%% 处理网格和边界
% 从外部读取网格
MyMesh = GmshReader(2,"QUAD");
% 自行创建网格
% NumX = 100;
% NumY = 100;
% LengthX = 1;
% LengthY = 1;
% MyMesh = MeshCreater(NumX,NumY,LengthX,LengthY);
% 创建边界条件
MyBoundary = BoundaryCreater(MyMesh);
%% 创建设计变量，优化问题
% 设定优化问题的基本参数
NumDesign = MyMesh.NumCells;
% 为了简单使用一组体积约束
NumConst = 1;
% 过滤半径此处以当前网格大小的倍数作为
% 过滤半径，需要除以(2*sqrt(3))^2
Rmin = 3;
Rmin = (Rmin^2)/12;
% SIMP法的材料属性
KappaMax = 100;
KappaMin = 1;
% 惩罚系数
penality = 3;
% flag
MatrixEvaluteFlag = 2;
VectorEvaluteFlag = 1;
% 全部网格的面积
MeshMeasure = sum(MyMesh.CellMeasure,"all");
volfrac = 0.3;
% 创建设计变量
Simp_Rho = ones(MyMesh.NumCells,1);
Simp_Rho_Filtered = ones(MyMesh.NumCells,1);
Simp_Rho_Max = ones(MyMesh.NumCells,1);
Simp_Rho_Min = 1e-3 * ones(MyMesh.NumCells,1);
% MMA相关的变量
MMAA = 1;
MMAC = 1e3;
Simp_Old1 = Simp_Rho;
Simp_Old2 = Simp_Old1;
MMALow =  Simp_Rho_Min;
MMAUpp =  Simp_Rho_Max;
% 创建目标函数的导数和约束函数的导数
Object_Diff_Value = zeros(MyMesh.NumCells,1);
Constraint_Diff_Value = zeros(MyMesh.NumCells,NumConst);
% 因为只有一个体积约束
Constraint_Diff_Value = MyMesh.CellMeasure;
%% 初始化PDE过滤以及过滤体积约束
% 初始化PDE过滤矩阵
HelmholtzStruct = struct();
HelmholtzStruct = HelmholtzEquation(max(MyMesh.CellMeasure) * Rmin * ones(MyMesh.NumCells,1), ...
                  MyMesh,Constraint_Diff_Value, ...
                  MatrixEvaluteFlag,HelmholtzStruct);
% 求解过滤方程
HelmholtzStruct.Solution = HelmholtzStruct.Matrix \ HelmholtzStruct.RHS;
HelmholtzStruct.Solution = HelmholtzStruct.Solution(HelmholtzStruct.IRCMMap);
% 重新生成过滤后的向量
Constraint_Diff_Value = FilterReGenerater(MyMesh,HelmholtzStruct.Solution);
% 归一化
Constraint_Diff_Max = max(abs(Constraint_Diff_Value));
Constraint_Diff_Value = Constraint_Diff_Value/Constraint_Diff_Max;
%% 拓扑优化循环
change = 1e10;
LoopIter = 0;
LoopMax = 50;
Constraint_Value = zeros(LoopMax,NumConst);
Object_Value = zeros(LoopMax,NumConst);
%% 输出结果
MyCellData = struct();
MyCellData.Size = 4;
MyCellData.Name = ["SimpOut","SimpOutFilter","ObjectDiff","ConstraintDiff"];
MyCellData.Data = zeros(MyMesh.NumCells,MyCellData.Size);
filename = "VTKOutput/Simp-Out-";
%% 
while (change > 1e-3 * MyMesh.NumCells && LoopIter <= LoopMax)
    LoopIter = LoopIter + 1;
    %% 过滤密度
    HelmholtzStruct = HelmholtzEquation(max(MyMesh.CellMeasure) * Rmin * ones(MyMesh.NumCells,1), ...
                      MyMesh,Simp_Rho, ...
                      VectorEvaluteFlag,HelmholtzStruct);
    % 求解过滤方程
    HelmholtzStruct.Solution = HelmholtzStruct.Matrix \ HelmholtzStruct.RHS;
    HelmholtzStruct.Solution = HelmholtzStruct.Solution(HelmholtzStruct.IRCMMap);
    % 重新生成过滤后的向量
    Simp_Rho_Filtered = FilterReGenerater(MyMesh,HelmholtzStruct.Solution);
    %% 物理场求解
    % 计算Laplace矩阵和右端项
    Simp_Rho_Input = KappaMin + (KappaMax - KappaMin) * Simp_Rho_Filtered .^ penality; 
    MyLaplace=LaplaceEquation(Simp_Rho_Input,MyMesh,MyBoundary);
    % 求解线性系统
    MyLaplace.Solution = MyLaplace.Matrix \ MyLaplace.RHS;
    SolutionIRCM = MyLaplace.Solution(MyLaplace.IRCMMap);
    % 计算目标函数的敏度
    Object_Diff_Value = SensitivityAnalysis(SolutionIRCM,Simp_Rho,MyMesh,penality,(KappaMax - KappaMin));
    %% 过滤目标函数的敏度
    HelmholtzStruct = HelmholtzEquation(MyMesh.CellMeasure * Rmin, ...
                      MyMesh,Object_Diff_Value, ...
                      VectorEvaluteFlag,HelmholtzStruct);
    % 求解过滤方程
    HelmholtzStruct.Solution = HelmholtzStruct.Matrix \ HelmholtzStruct.RHS;
    HelmholtzStruct.Solution = HelmholtzStruct.Solution(HelmholtzStruct.IRCMMap);
    % 重新生成过滤后的向量
    Object_Diff_Value = FilterReGenerater(MyMesh,HelmholtzStruct.Solution);
    Object_Diff_Value = Object_Diff_Value/max(abs(Object_Diff_Value)); 
    %% 计算目标函数和约束函数的值
    Constraint_Value(LoopIter) = dot(Simp_Rho,MyMesh.CellMeasure) - volfrac * MeshMeasure;
    Constraint_Value(LoopIter) = Constraint_Value(LoopIter)/Constraint_Diff_Max;
    Object_Value(LoopIter) = abs(MyLaplace.Solution' * MyLaplace.RHS);
    %% MMA更新设计变量
    [Simp_New,~]=MMA_Dual_ZXL(Simp_Rho,Constraint_Value(LoopIter),Object_Diff_Value,Constraint_Diff_Value', ...
                 NumConst,NumDesign,Simp_Rho_Min,Simp_Rho_Max,1,MMAA,MMAC,1,LoopIter, ...
                 Simp_Old1,Simp_Old2,MMALow,MMAUpp);
    % 计算变动量
    change = sum(abs(Simp_Rho - Simp_New));
    Simp_Old2 = Simp_Old1;
    Simp_Old1 = Simp_Rho;
    Simp_Rho = Simp_New;
    % 输出结果
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',LoopIter,Object_Value(LoopIter), ...
    Constraint_Value(LoopIter)*Constraint_Diff_Max/MeshMeasure + volfrac,change);
    % 写入文件
    MyCellData.Data(:,1) = Simp_Rho;
    MyCellData.Data(:,2) = Simp_Rho_Filtered;
    MyCellData.Data(:,3) = Object_Diff_Value;
    MyCellData.Data(:,4) = Constraint_Diff_Value;
    VTKWriter(filename + num2str(LoopIter) + ".vtk",MyMesh,MyCellData);
end
%% 画图
Index = 1:LoopIter;
figure("Name","Object Constraint Values")
Object_Value = Object_Value/max(abs(Object_Value));
plot(Index,Object_Value(1:LoopIter),"Color","black","Marker","+");
grid on
hold on 
Constraint_Value = Constraint_Value*Constraint_Diff_Max/MeshMeasure + volfrac;
plot(Index,Constraint_Value(1:LoopIter),"Color","blue","Marker","o");
hold off
%
fclose all;