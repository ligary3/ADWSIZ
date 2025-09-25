import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mgwr.gwr import GWR # 导入GWR类，用于地理加权回归，支持岭回归
from mgwr.sel_bw import Sel_BW # 导入Sel_BW类，用于选择最佳带宽
from sklearn.preprocessing import StandardScaler # 导入StandardScaler，用于数据标准化
from sklearn.metrics import mean_squared_error # 导入mean_squared_error，用于计算均方误差
import warnings # 导入warnings模块，用于控制警告信息

warnings.filterwarnings("ignore", category=UserWarning) # 忽略无关警告，使输出更简洁

# 设置SCI论文绘图参数
plt.rcParams.update({
    'font.family': 'serif', # 设置字体族为衬线字体
    'font.serif': ['Times New Roman'], # 具体指定衬线字体为Times New Roman
    'font.size': 10, # 全局字体大小
    'axes.labelsize': 10, # 轴标签字体大小
    'axes.titlesize': 12, # 轴标题字体大小
    'legend.fontsize': 8, # 图例字体大小
})

# 加载数据
df = pd.read_csv(r"E:\30m\景观格局指数\圆1-62\站点.csv")#包含景观格局指数结果
# 提取坐标（经度和纬度），mgwr库会根据这些坐标计算空间权重
coords = df[['经度', '纬度']].values
# 定义目标变量（因变量）
target = 'NO2_24h'

# 定义不同缓冲区尺度
scales = ['500m', '1000m', '1500m', '2000m', '2500m', '3000m', '3500m', '4000m', '4500m', '5000m', '5500m', '6000m']

# 存储模型结果的DataFrame
results_df = pd.DataFrame(columns=['Scale', 'R2', 'Adj_R2', 'AICc', 'RMSE', 'Best_Lambda', 'Used_Vars'])

# 遍历每个缓冲区尺度进行分析
for scale in scales:
    print(f"\n==== 处理 {scale} 缓冲区 ====")
    
    # 筛选当前尺度对应的自变量（景观格局指数）
    scale_vars = [col for col in df.columns if col.endswith(f"_{scale}")]
    
    # 数据预处理：删除包含NaN值的行
    df_clean = df.dropna(subset=scale_vars + [target])
    # 提取自变量X和因变量y
    X = df_clean[scale_vars].values # 将自变量转换为NumPy数组
    y = df_clean[target].values.reshape(-1, 1) # 将因变量转换为NumPy数组，并重塑为列向量
    # 提取清理后的坐标
    coords_clean = df_clean[['经度', '纬度']].values
    # 获取样本数量和特征数量
    n_samples, n_features = X.shape
    
    # 最终使用的变量列表
    used_vars = scale_vars
    X = df_clean[used_vars].values # 再次确认X使用清理后的变量
    n_features = X.shape[1] # 再次确认特征数量
    print(f"使用变量: {used_vars}, 特征数: {n_features}")
    
    # 标准化自变量，以消除量纲影响
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X) # 对每个特征列单独进行标准化
    
    # 1. 选择GWR的最佳带宽 (Best Bandwidth)
    # Sel_BW用于选择最优带宽，这里使用AICc准则
    # fixed=False 表示带宽是适应性的（每个点有自己的带宽），而不是固定距离
    selector = Sel_BW(coords_clean, y, X_scaled, fixed=False)
    best_bw = selector.search(criterion='AICc', verbose=False) # 执行带宽搜索
    print(f"最佳带宽: {best_bw:.2f}")
    
    # 2. 定义岭回归参数λ的搜索范围
    lambda_grid = np.logspace(-3, 3, 20) # 在1e-3到1e3之间生成20个对数均匀分布的λ值
    
    # 3. 网格搜索最优λ（基于全局RMSE）
    best_lambda = None # 初始化最优λ
    best_rmse = np.inf # 初始化最小RMSE为无穷大
    
    # 遍历每个λ值，找到使模型RMSE最小的λ
    for lambda_val in lambda_grid:
        # 初始化GWR模型，并传入当前λ值作为岭回归参数
        # mgwr.gwr.GWR 会在内部处理空间权重和岭回归的结合
        model = GWR(coords_clean, y, X_scaled, bw=best_bw, fixed=False, ridge=lambda_val)
        # 拟合模型
        results = model.fit()
        
        # 计算当前λ值下的RMSE
        current_rmse = np.sqrt(mean_squared_error(y, results.predictions))
        
        # 更新最优λ和最小RMSE
        if current_rmse < best_rmse:
            best_rmse = current_rmse
            best_lambda = lambda_val
    
    print(f"最优λ: {best_lambda:.4f}, 对应的RMSE: {best_rmse:.3f}")
    
    # 4. 使用最优λ重新训练最终的GWRR模型
    final_model = GWR(coords_clean, y, X_scaled, bw=best_bw, fixed=False, ridge=best_lambda)
    final_results = final_model.fit()
    
    # 5. 模型评估
    # 直接从mgwr模型结果中获取R2、调整R2和AICc
    r2 = final_results.R2
    adj_r2 = final_results.adj_R2
    aicc = final_results.aicc
    # 从预测结果计算RMSE
    rmse = np.sqrt(mean_squared_error(y, final_results.predictions))
    
    # 存储当前尺度的模型结果
    temp_df = pd.DataFrame({
        'Scale': [scale],
        'R2': [r2],
        'Adj_R2': [adj_r2],
        'AICc': [aicc],
        'RMSE': [rmse],
        'Best_Lambda': [best_lambda],
        'Used_Vars': [','.join(used_vars)] # 将使用的变量列表转换为字符串
    })
    # 将当前尺度的结果添加到总结果DataFrame中
    results_df = pd.concat([results_df, temp_df], ignore_index=True)

# 可视化模型性能对比
# 创建一个2x2的子图布局
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
# 定义要可视化的指标及其颜色
metrics = [('R2', 'skyblue'), ('Adj_R2', 'lightgreen'), ('AICc', 'salmon'), ('RMSE', 'gold')]

# 遍历每个子图和指标
for ax, (metric, color) in zip(axes.flatten(), metrics):
    ax.bar(results_df['Scale'], results_df[metric], color=color) # 绘制柱状图
    ax.set_title(f'{metric} Comparison') # 设置子图标题
    # 根据指标类型设置Y轴范围
    if metric in ['R2', 'Adj_R2']:
        ax.set_ylim(0, 1) # R2和Adj_R2范围在0到1之间
    elif metric == 'RMSE':
        ax.set_ylim(bottom=0) # RMSE从0开始

plt.tight_layout() # 自动调整子图布局，防止重叠
plt.show() # 显示图形

# 输出最终结果
print("\nGWRR模型性能汇总:")
# 打印关键模型性能指标，并保留3位小数
print(results_df[['Scale', 'R2', 'Adj_R2', 'AICc', 'RMSE', 'Best_Lambda']].round(3))
# 找到调整R²最大的模型（即最优模型）
best_idx = results_df['Adj_R2'].idxmax()
# 打印最优模型的详细信息
print(f"\n最优模型: 尺度={results_df.loc[best_idx,'Scale']}, 调整R²={results_df.loc[best_idx,'Adj_R2']:.3f}")
print(f"使用变量: {results_df.loc[best_idx,'Used_Vars']}")
