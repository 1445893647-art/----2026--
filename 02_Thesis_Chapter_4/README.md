# 大论文第四章：面向空中动态目标的无人机反隐蔽通信方案支撑代码

## 1. 简介

本代码是学位论文第四章的支撑代码。主要实现了在动态目标场景下，基于**扩展卡尔曼滤波 (EKF)** 进行非法接收者 (Bob) 轨迹跟踪，并结合交替优化 (AO/BCD) 或连续凸近似 (SCA) 算法，对无人机的**三维飞行轨迹、干扰功率及发射波束成形**进行联合优化的隐蔽通信方案。

## 2. 环境配置

* **编程语言**：MATLAB 2021a 或更高版本
* **核心依赖工具箱**：
  * **CVX Toolbox**：必需。用于求解轨迹优化（`solve_trajectory_opt`）和波束成形优化（基于 LMI/Schur 补的 `solve_beamforming_opt`）子问题。建议使用 Mosek 求解器以提高稳定性。
  * **Sensor Fusion and Tracking Toolbox / Phased Array System Toolbox**：建议预装，辅助雷达测量模型与阵列处理计算。

## 3. 文件结构说明

### 主控仿真脚本 (Main Scripts)
* `main_ekf_sim.m`：**核心推荐脚本**。基于 EKF + 交替优化 (AO) 的联合优化完整流程。
* `main_ekf_sim_damped_3iter.m`：**阻尼版 AO**。引入阻尼系数 `alpha` 并放宽阈值，抑制优化振荡，提升早期时隙的可行性。
* `main_ekf_sim_SCA_with_iter_plot.m`：**SCA 联合优化版**。基于信赖域和回溯机制进行 SCA 迭代，含 Oracle 与 Lag-1 对比。
* `main_dep_sim.m`：**参数扫描脚本**。用于遍历检测错误概率阈值 $\delta$（代码中为 `xi_req`），统计平均中断概率。

### 核心优化子模块 (Optimization Solvers)
* `solve_power_opt.m`：功率优化子问题。基于检测约束计算隐蔽约束下的干扰功率上界（闭式解）。
* `solve_beamforming_opt.m`：波束成形优化子问题。调用 CVX 求解，包含 soft PCRB 约束（松弛变量）。
* `solve_trajectory_opt.m`：轨迹优化子问题。调用 CVX 处理无人机速度、高度、到 Alice 距离及目标跟踪 PCRB 等约束。

### 目标跟踪与测量模型 (EKF & Measurements)
* `measurement_function.m` / `calculate_jacobian.m`：计算雷达回波、时延、多普勒的非线性测量函数及其对 EKF 的雅可比矩阵计算。
* `get_angle_dist.m` / `get_steering_vec.m` / `get_steering_deriv.m`：几何距离计算、UPA 阵列导向矢量及关于角度偏导的生成函数。

### 全局配置与基准工具
* `config_param.m` & `config_param_sca.m`：全局参数入口（物理参数、通道参数、SCA 信赖域参数等）。
* `get_gamma_req.m` / `calc_metrics.m`：基于 LambertW 函数反解检测阈值，并计算中断概率 $P_{out}$ 等性能指标。
* `run_noekf_opt.m` / `build_fixed_trajectory.m`：不包含 EKF 的基准对比（Oracle/Lag-1）与固定直线轨迹生成工具。

---

## 4. 论文图表复现指南

运行以下脚本以复现论文中的对应图表：

| 论文图题                                    | 对应主脚本与运行说明                                         |
| :------------------------------------------ | :----------------------------------------------------------- |
| **算法收敛性能图**                          | 运行 `main_ekf_sim_damped_3iter.m` (观察交替优化的收敛) 或 `main_ekf_sim_SCA_with_iter_plot.m` (生成 `SCA Convergence` 图表) |
| **无人机三维飞行轨迹**                      | 运行 `main_ekf_sim.m` 直接生成三维跟踪避让轨迹图             |
| **EKF 跟踪效果图**                          | 运行 `main_ekf_sim.m`，自动输出真实轨迹与 EKF 预测/更新轨迹的对比 |
| **$P_{out}$ 随时隙变化图**                  | 运行 `main_ekf_sim.m` (生成 Proposed 曲线)。对比方案操作见下文第 5 节。 |
| **平均 $P_{out}$ 随检测阈值 $\delta$ 变化** | 运行 `main_dep_sim.m`。若需加入 Oracle/Lag-1 基准，取消脚本末尾被注释的代码即可。 |
| **平均 $P_{out}$ 随 Alice 功率变化**        | 修改 `main_dep_sim.m` 外层循环（将扫描 $\delta$ 改为扫描 `P.P_a_dB`），详见下方第 5 节。 |

---

## 5. 如何运行对比基准与扩展仿真

为了生成全面的对比数据，论文设计了三类退化基线方案。您可以在主脚本的基础上通过“冻结变量”的方式实现：

### 5.1 三类基线方案配置 (OTOP-FB / OBOP-FT / OBOT-FP)

1. **OTOP-FB (优化轨迹+功率，波束固定)**
   * **操作**：在主循环中，跳过调用 `solve_beamforming_opt`。直接使用 `get_steering_vec` 固定计算波束向量 `w_wb`（使其始终机械对准 Bob 的预测/真实位置），并仅更新轨迹和功率。
2. **OBOP-FT (优化波束+功率，轨迹固定)**
   * **操作**：在时隙循环前，调用 `q_w_fixed = build_fixed_trajectory(...)` 预先生成完整固定轨迹。
   * 在时隙循环内部，强制赋值 `q_w_curr = q_w_fixed`。仅交替调用 `solve_power_opt` 和 `solve_beamforming_opt`。
3. **OBOT-FP (优化波束+轨迹，功率固定)**
   * **操作**：强制设定干扰功率为常数 `P_J_curr = P.P_J`（查阅 `config_param` 中的初始设定值）。
   * 在主循环中仅调用 `solve_trajectory_opt` 和 `solve_beamforming_opt`，不再更新功率。

### 5.2 扩展仿真：扫描 Alice 发射功率 ($P_a$)

若需生成**“平均 $P_{out}$ 随 Alice 发射功率变化”**图，请基于 `main_dep_sim.m` 进行修改。将原外层遍历检测阈值的部分替换为如下代码：

```matlab
Pa_dB_list = -10:2:20;               % 定义发射功率扫描范围
avg_pout = zeros(size(Pa_dB_list)); 

for i = 1:length(Pa_dB_list)
    P.P_a_dB = Pa_dB_list(i);
    P.P_a    = 10^(P.P_a_dB/10);     % 转换为线性值
    
    % --- 在此调用单次完整时隙的仿真主体(如 main_ekf_sim 的逻辑) ---
    % avg_pout(i) = ... (获取当前配置下的平均中断概率)
end

figure;
plot(Pa_dB_list, avg_pout, '-o', 'LineWidth', 1.5);
xlabel('Alice Transmit Power P_a (dBm)');
ylabel('Average Outage Probability');
grid on;
```