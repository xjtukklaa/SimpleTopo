%% 用于计算不规则的网格的JacobianInverse和JacobianDet
function [JacobianInverse,JacobianDet]=JacobianWeight(OrderedPoints)
    % Jacobian矩阵的逆
    JacobianInverse = zeros(4,2,2);
    % Jacobian矩阵的行列式的逆
    JacobianDet = zeros(4,1);
    % 权重都是1
    % 正交点,旋转顺序为(-1,-1)->(1,-1)->(1,1)->(-1,1)
    Q_Points = [-1/sqrt(3),-1/sqrt(3);
                 1/sqrt(3),-1/sqrt(3);
                 1/sqrt(3), 1/sqrt(3);
                -1/sqrt(3), 1/sqrt(3)];
    for QIndex = 1:4
        zeta = Q_Points(QIndex, 1);
        eta  = Q_Points(QIndex, 2);
        % 形函数对zeta的导数
        dNdzeta = [-(1 - eta);
                    (1 - eta);
                    (1 + eta);
                   -(1 + eta)]/4;
        % 形函数对eta的导数
        dNdeta  = [-(1 - zeta);
                   -(1 + zeta);
                    (1 + zeta);
                    (1 - zeta)]/4;
        % 这里是正常的Jacobian
        JacobianInverse(QIndex,:,:) = [dNdzeta';dNdeta'] * OrderedPoints;
        % Jacobian矩阵的行列式
        JacobianDet(QIndex) = JacobianInverse(QIndex,1,1) * JacobianInverse(QIndex,2,2) -...
                              JacobianInverse(QIndex,1,2) * JacobianInverse(QIndex,2,1);
        if abs(JacobianDet(QIndex)) < 1e-10
            error("Error : Jacobian NullSpace !!");
        end
        % 计算Jacobian矩阵的逆
        JacobianTmp = JacobianInverse(QIndex,1,1);
        JacobianInverse(QIndex,1,1) = JacobianInverse(QIndex,2,2);
        JacobianInverse(QIndex,2,2) = JacobianTmp;
        JacobianInverse(QIndex,2,1) = -JacobianInverse(QIndex,2,1);
        JacobianInverse(QIndex,1,2) = -JacobianInverse(QIndex,2,1);
        JacobianInverse(QIndex,:,:) = JacobianInverse(QIndex,:,:) / ...
                                      JacobianDet(QIndex);
    end
    % 计算结束
end