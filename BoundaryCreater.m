function BoundaryStruct = BoundaryCreater(Mesh)
    BoundaryKey = ["Dirichlet","Dirichlet_Value"];
    % 内部网格创建
    % BoundaryDirichletDofs = cell2mat(Mesh.BoundaryDofsMap("Down"));
    % BoundaryPoints = zeros(1,2);
    % BoundaryDirichletDofsChange = zeros(length(BoundaryDirichletDofs),1);
    % for Index = 1:length(BoundaryDirichletDofs)
    %     BoundaryPoints(:) = Mesh.PointsMap(BoundaryDirichletDofs(Index),:);
    %     if BoundaryPoints(1) >= 0.45 && BoundaryPoints(1) <= 0.55
    %         BoundaryDirichletDofsChange(Index) =  1;
    %     end
    % end
    % BoundaryDirichletDofs = BoundaryDirichletDofs(find(BoundaryDirichletDofsChange == 1));
    % 外部网格输入
    BoundaryDirichletDofs = cell2mat(Mesh.BoundaryDofsMap("Dirichlet"));
    % 
    BoundaryDirichletValue = 1;
    BoundaryStruct = dictionary(BoundaryKey,[{BoundaryDirichletDofs},BoundaryDirichletValue]);
end