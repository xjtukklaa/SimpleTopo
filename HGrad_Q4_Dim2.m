%% OrderedPoints排序后的节点,Flag=1为计算形函数,Flag=2为计算形函数梯度和值
%% 返回当前单元上的形函数的值,形函数的梯度值,(Dofs,QPoints)
function [GradValue, ShapeValue] = HGrad_Q4_Dim2(EvaluteFlags)
    % 4个形函数在正交点的值
    ShapeValue = zeros(4, 4);
    % 4个形函数在正交点的梯度值,[dN_dzeta,dN_deta]'
    GradValue  = zeros(4, 4, 2);
    % 正交点,旋转顺序为(-1,-1)->(1,-1)->(1,1)->(-1,1)
    Q_Points = [-1/sqrt(3),-1/sqrt(3);
                 1/sqrt(3),-1/sqrt(3);
                 1/sqrt(3), 1/sqrt(3);
                -1/sqrt(3), 1/sqrt(3)];
    % 循环所有正交点,计算形函数的值
    for QIndex = 1:4
        zeta = Q_Points(QIndex, 1);
        eta  = Q_Points(QIndex, 2);
        ShapeValue(1, QIndex) = (1 - zeta) * (1 - eta)/4;
        ShapeValue(2, QIndex) = (1 + zeta) * (1 - eta)/4;
        ShapeValue(3, QIndex) = (1 + zeta) * (1 + eta)/4;
        ShapeValue(4, QIndex) = (1 - zeta) * (1 + eta)/4;
    end
    % 循环所有正交点,计算形函数梯度的值
    if EvaluteFlags == 2
        for QIndex = 1:4
            zeta = Q_Points(QIndex, 1);
            eta  = Q_Points(QIndex, 2);
            GradValue(1, QIndex, :) = [-(1 - eta), -(1 - zeta)]/4;
            GradValue(2, QIndex, :) = [ (1 - eta), -(1 + zeta)]/4;
            GradValue(3, QIndex, :) = [ (1 + eta),  (1 + zeta)]/4;
            GradValue(4, QIndex, :) = [-(1 + eta),  (1 - zeta)]/4;
        end
    end
    % 计算完成
end