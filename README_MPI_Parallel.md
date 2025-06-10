# 基于MPI并行计算的热传递拓扑优化代码

## 概述

这是一个针对问题三（带热传递的管道弯曲问题）进行MPI并行优化的MATLAB代码。该代码使用MATLAB的Parallel Computing Toolbox来加速有限元计算，特别是针对热传递-流体耦合问题的计算密集型部分。

## 主要特性

### 1. 并行计算增强
- **元素级并行化**：使用`parfor`循环并行计算所有有限元操作
- **残差和雅可比矩阵并行计算**：牛顿求解器中的核心计算部分
- **目标函数并行评估**：并行计算每个元素的目标函数贡献
- **敏感度分析并行化**：并行计算设计敏感度

### 2. 热传递问题专用优化
- **温度场求解**：支持带热源的稳态热传导方程
- **流体-热耦合**：考虑对流换热和热传导
- **热敏感度计算**：包含温度方程对设计变量的敏感度

### 3. 自适应优化算法
- **自适应移动限制**：根据收敛特性动态调整优化步长
- **β投影连续化**：自动调整Heaviside投影参数
- **多阶段延拓策略**：渐进式优化求解

## 系统要求

### 必需软件
1. **MATLAB** (R2020a或更高版本)
2. **MATLAB Parallel Computing Toolbox**
3. **足够的系统内存**（推荐16GB+）
4. **多核处理器**（推荐4核+）

### 依赖文件
确保以下文件在MATLAB路径中：
- `RES.m` - 残差函数
- `JAC.m` - 雅可比矩阵函数
- `PHI.m` - 目标函数
- `dPHIds.m` - 目标函数对状态变量的导数
- `dPHIdg.m` - 目标函数对设计变量的导数
- `dRESdg.m` - 残差对设计变量的导数
- `dRESTemperature.m` - 温度方程对设计变量的导数
- `problems.m` - 问题定义文件
- `postproc.m` - 后处理文件

## 使用方法

### 1. 启动并行计算
```matlab
% 启动MATLAB并运行
topFlow_mpi
```

### 2. 并行池配置
代码会自动检测并启动并行池：
```matlab
% 检查并启动并行池
p = gcp('nocreate');
if isempty(p)
    parpool('local');  % 使用本地工作器
    p = gcp;
end
```

### 3. 自定义并行设置
如需自定义并行设置，可在代码开始前手动启动：
```matlab
% 指定工作器数量（可选）
parpool('local', 4);  % 使用4个工作器

% 然后运行主程序
topFlow_mpi
```

## 性能优化建议

### 1. 硬件配置
- **CPU核心数**：更多核心 = 更好的并行性能
- **内存容量**：确保足够内存避免交换分区使用
- **存储速度**：SSD比HDD有更好的I/O性能

### 2. 网格尺寸调整
```matlab
% 在topFlow_mpi.m中调整网格密度
nely = 60;  % 增大此值可提高精度，但增加计算成本
nelx = nely*Lx/Ly;
```

### 3. 优化参数调整
```matlab
% 调整收敛容差
nltol = 1e-6;  % 牛顿求解器容差
chlim = 1e-3;  % 设计变化容差

% 调整最大迭代次数
maxiter = qanum*conit + beta_num*beta_conit;
```

## 关键并行化部分

### 1. 残差计算并行化
```matlab
% 并行计算所有元素的残差
parfor i = 1:neltot
    sR(:,i) = RES(dxv(i),dyv(i),muv(i),rhov(i),kv(i),Cpv(i),Qv(i),alpha(i),...
        uVars(i,1:8), pVars(i,1:4), TVars(i,1:4));
end
```

### 2. 雅可比矩阵并行化
```matlab
% 并行计算雅可比矩阵
parfor i = 1:neltot
    sJ(:,i) = JAC(dxv(i),dyv(i),muv(i),rhov(i),kv(i),Cpv(i),Qv(i),alpha(i),...
        uVars(i,1:8), pVars(i,1:4), TVars(i,1:4));
end
```

### 3. 敏感度分析并行化
```matlab
% 并行计算设计敏感度
parfor i = 1:neltot
    drdgVals(:,i) = dRESdg(...);
    dphidgVals(i) = dPHIdg(...);
end
```

## 监控和调试

### 1. 性能监控
```matlab
% 查看并行池信息
p = gcp;
fprintf('使用%d个工作器\n', p.NumWorkers);

% 监控内存使用
memory
```

### 2. 调试模式
```matlab
% 如需调试，可暂时关闭并行计算
% 将parfor改为for进行串行调试
for i = 1:neltot  % 替代 parfor i = 1:neltot
    % 计算代码
end
```

## 输出和结果

### 1. 计算进度显示
代码会显示详细的计算进度：
```
=========================================================
      Problem number:  3 - Reynolds number: 1.00e+03
      Heat transfer problem with uniform heat source
      Parallel Computing Enabled
      Number of workers: 4
=========================================================
```

### 2. 性能统计
最终会显示并行计算性能统计：
```
=========================================================
      Number of design iterations:   25
      Final objective: 1.23e+02
      Total time taken:  12.34 min
      Parallel Performance:
      Workers used: 4
      Total processes: 5
=========================================================
```

## 故障排除

### 1. 并行池启动失败
```matlab
% 检查Parallel Computing Toolbox许可证
license('test', 'Distrib_Computing_Toolbox')

% 手动启动并行池
delete(gcp('nocreate'));  % 删除现有池
parpool('local');         % 重新启动
```

### 2. 内存不足
- 减小网格尺寸（降低`nely`值）
- 增加系统虚拟内存
- 减少并行工作器数量

### 3. 函数文件缺失
确保所有必需的`.m`文件在MATLAB路径中：
```matlab
% 检查文件是否存在
which RES
which JAC
which PHI
```

## 参考和引用

基于以下原始代码开发：
- J. Alexandersen, "A detailed introduction to density-based topology optimisation of fluid flow problems including implementation in MATLAB", SMO 2022

并行化增强：
- 使用MATLAB Parallel Computing Toolbox
- 元素级并行化策略
- 针对热传递问题的优化

## 许可证

遵循原始代码的BSD 3-Clause License。 