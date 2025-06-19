function VTKWriter(filename,Mesh,CellData)
    VTK = fopen(filename, 'w'); 
    % VTK版本
    fprintf(VTK, '# vtk DataFile Version 2.0\n');
    fprintf(VTK, 'VTK Output Mesh Data\n');
    % VTK类型
    fprintf(VTK, 'ASCII\n');
    % 网格信息
    fprintf(VTK, 'DATASET UNSTRUCTURED_GRID\n');
    % 网格点坐标
    fprintf(VTK, ['POINTS ' , num2str(Mesh.NumDofs) , ' double\n']);
    if Mesh.Dim == 2
        for CellIndex = 1:Mesh.NumDofs
            for DimIndex = 1:2
                fprintf(VTK, num2str(Mesh.PointsMap(CellIndex,DimIndex)));
                fprintf(VTK, ' ');
            end
            fprintf(VTK, num2str(0));
            fprintf(VTK, '\n');
        end
        fprintf(VTK, '\n');
    else
        for CellIndex = 1:Mesh.NumDofs
            for DimIndex = 1:Mesh.Dim
                fprintf(VTK, num2str(Mesh.PointsMap(CellIndex,DimIndex)));
                fprintf(VTK, ' ');
            end
            fprintf(VTK, '\n');
        end
        fprintf(VTK, '\n');
    end
    % 网格单元
    fprintf(VTK, ['CELLS ' , num2str(Mesh.NumCells) , ' ' , num2str(Mesh.NumCells * 5) , '\n']);
    for CellIndex = 1:Mesh.NumCells
        fprintf(VTK, '4');
        fprintf(VTK, ' ');
        for DofIndex = 1:3
            fprintf(VTK, num2str(Mesh.DofsMap(CellIndex,DofIndex)-1));
            fprintf(VTK, ' ');
        end
        fprintf(VTK, num2str(Mesh.DofsMap(CellIndex,4)-1));
        fprintf(VTK, '\n');
    end
    fprintf(VTK, '\n');
    % 网格类型
    fprintf(VTK, ['CELL_TYPES ' , num2str(Mesh.NumCells) , '\n']);
    for CellIndex = 1:Mesh.NumCells
        fprintf(VTK, '9');
        fprintf(VTK, '\n');
    end
    fprintf(VTK, '\n');
    % 网格上的数据
    fprintf(VTK, ['CELL_DATA ' , num2str(Mesh.NumCells) , '\n']);
    for DataIndex = 1:CellData.Size
        fprintf(VTK, ['SCALARS ', num2str(CellData.Name(DataIndex)), ' float 1 \n']);
        fprintf(VTK,'LOOKUP_TABLE default\n');
        for CellIndex = 1:Mesh.NumCells
            fprintf(VTK, num2str(CellData.Data(CellIndex,DataIndex),1e-8));
            fprintf(VTK, '\n');
        end
        fprintf(VTK, '\n');
    end
    % 
end
    