#---------------------------------半变异函数分析------------------------
import rasterio # 导入rasterio库，用于处理栅格数据
import numpy as np # 导入numpy库，用于数值计算
import matplotlib.pyplot as plt # 导入matplotlib.pyplot库，用于绘图
from skgstat import Variogram # 从skgstat库导入Variogram类，用于半变异函数计算
import skgstat as skg # 导入skgstat库，用于地统计分析
from tabulate import tabulate # 导入tabulate库，用于美化表格输出
import gc # 导入gc库，用于垃圾回收，管理内存

# 配置SCI论文绘图参数 - 进一步增大字体和调整线条
plt.rcParams.update({
    'figure.dpi': 300, # 设置图形分辨率为300 DPI
    'font.family': 'Times New Roman', # 设置字体为Times New Roman
    'font.size': 12,  # 全局基础字体增大
    'axes.labelsize': 14, # 轴标签字体增大
    'axes.labelweight': 'bold', # 轴标签字体加粗
    'xtick.labelsize': 12, # X轴刻度标签字体增大
    'ytick.labelsize': 12, # Y轴刻度标签字体增大
    'axes.spines.top': True, # 显示顶部轴脊
    'axes.spines.right': True # 显示右侧轴脊
})

def load_geotiff(path):
    """
    加载GeoTIFF文件并返回数据、仿射变换矩阵和NoData值。
    参数:
        path (str): GeoTIFF文件的路径。
    返回:
        tuple: 包含栅格数据 (numpy.ndarray), 仿射变换矩阵 (rasterio.Affine), 和NoData值。
    """
    with rasterio.open(path) as src: # 以只读模式打开GeoTIFF文件
        data = src.read(1) # 读取第一个波段的数据
        transform = src.transform # 获取仿射变换矩阵
        nodata = src.nodata # 获取NoData值
    return data, transform, nodata

def preprocess_data(data, nodata):
    """
    数据预处理：将NoData值替换为NaN，并提取有效数据的行、列索引和值。
    参数:
        data (numpy.ndarray): 栅格数据。
        nodata (int/float): 栅格数据的NoData值。
    返回:
        tuple: 包含有效数据的行索引 (numpy.ndarray), 列索引 (numpy.ndarray), 和值 (numpy.ndarray)。
    """
    data = np.where(data == nodata, np.nan, data) # 将NoData值替换为NaN
    mask = ~np.isnan(data) # 创建一个布尔掩膜，标记非NaN值
    rows, cols = np.where(mask) # 获取非NaN值的行和列索引
    values = data[rows, cols] # 提取非NaN值
    return rows, cols, values

def calculate_coordinates(rows, cols, transform):
    """
    根据行、列索引和仿射变换矩阵计算UTM坐标。
    参数:
        rows (numpy.ndarray): 数据的行索引。
        cols (numpy.ndarray): 数据的列索引。
        transform (rasterio.Affine): 栅格的仿射变换矩阵。
    返回:
        numpy.ndarray: 包含X和Y坐标的二维数组。
    """
    x_coords, y_coords = transform * (cols + 0.5, rows + 0.5) # 计算每个像素中心点的UTM坐标
    return np.column_stack((x_coords, y_coords)) # 将X和Y坐标堆叠成一个二维数组

def downsample_data(coordinates, values, sample_size=10000):
    """
    数据下采样：如果数据量超过指定样本大小，则随机抽取子集。
    参数:
        coordinates (numpy.ndarray): 坐标数据。
        values (numpy.ndarray): 值数据。
        sample_size (int): 目标样本大小。
    返回:
        tuple: 下采样后的坐标和值。
    """
    if len(coordinates) > sample_size: # 如果数据点数量大于样本大小
        np.random.seed(5) # 设置随机种子，以便结果可复现
        indices = np.random.choice(len(coordinates), sample_size, replace=False) # 随机选择指定数量的索引
        return coordinates[indices], values[indices] # 返回下采样后的坐标和值
    return coordinates, values # 如果数据量不大，则直接返回原始数据

def fit_variogram_models(coordinates, values):
    """
    拟合不同类型的半变异函数模型（球状、指数、高斯、Matern）。
    参数:
        coordinates (numpy.ndarray): 坐标数据。
        values (numpy.ndarray): 值数据。
    返回:
        list: 包含每个模型拟合结果的字典列表。
    """
    models = ['spherical', 'exponential', 'gaussian', 'matern'] # 定义要拟合的模型列表
    results = [] # 初始化存储结果的列表
    
    for model in models: # 遍历每个模型
        try:
            # 拟合半变异函数模型
            V = Variogram(
                coordinates=coordinates, # 坐标
                values=values, # 值
                n_lags=50, # 滞后距离的数量
                maxlag=0.8, # 最大滞后距离（相对于数据范围的百分比）
                model=model, # 半变异函数模型类型
                estimator='matheron', # 估计器类型
            )
            
            # 获取模型参数 (范围, 基台, 块金)
            params = V.parameters
            
            # 计算模型的预测值
            pred = V.model(V.bins, *params)
            
            # 计算R²和SSE (残差平方和)
            exp = V.experimental # 实验半变异函数值
            mask = ~np.isnan(exp) # 创建掩膜，排除NaN值
            exp = exp[mask] # 过滤掉NaN值
            pred = pred[mask] # 过滤掉NaN值
            ss_res = np.sum((exp - pred) ** 2)  # 计算残差平方和 (SSE)
            ss_tot = np.sum((exp - np.mean(exp)) ** 2) # 计算总平方和
            r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else np.nan # 计算R²
            
            # 收集结果
            result = {
                'model': model, # 模型名称
                'nugget': params[2], # 块金值
                'range': params[0], # 范围值
                'sill': params[1], # 基台值
                'r_squared': r_squared, # R²值
                'variogram': V, # 半变异函数对象
                'sse': ss_res # 残差平方和
            }
            
            # 处理Matern模型的额外参数 (nu)
            if model == 'matern' and len(params) > 3:
                result['nu'] = params[3] # Matern模型的nu参数
                
            results.append(result) # 将结果添加到列表中
            
        except Exception as e: # 捕获异常
            print(f"Error fitting {model} model: {str(e)}") # 打印错误信息
            continue # 继续下一个模型
    
    return results

def plot_comparison(results, N):
    """
    绘制不同半变异函数模型的拟合比较图。
    参数:
        results (list): 包含模型拟合结果的字典列表。
        N (int): 样本数量。
    """
    # 调整figsize以适应双栏中的双栏宽度（整页宽度），并给下半部分留出空间
    # 通常整页宽度在6.5到7英寸，这里设置为6.8，高度设置为5.5，以便与下半部分组合
    plt.figure(figsize=(6.8, 5.5)) # 设置图形大小
    
    # 绘制实验值
    V = results[0]['variogram'] # 获取第一个模型的半变异函数对象（用于获取实验值）
    plt.scatter(V.bins, V.experimental, 
                color='k', alpha=0.5, label='Experimental', zorder=5, s=20) # 绘制散点图，表示实验半变异函数值
    
    # 生成平滑曲线的X轴范围
    x_smooth = np.linspace(0, V.bins.max(), 300)
    
    # 定义绘图样式
    colors = ['red','blue','green','purple'] # 曲线颜色
    linestyles = ['-', '--', '-.', ':'] # 曲线线型
    
    # 绘制拟合曲线
    for i, res in enumerate(results): # 遍历每个模型结果
        V = res['variogram'] # 获取当前模型的半变异函数对象
        y_smooth = V.model(x_smooth, *V.parameters) # 计算模型在平滑X轴范围上的预测值
        plt.plot(x_smooth, y_smooth,
                 color=colors[i], # 颜色
                 linestyle=linestyles[i], # 线型
                 linewidth=2.5, # 线宽
                 label=f"{res['model'].capitalize()}\nR²={res['r_squared']:.2f}, SSE={res['sse']:.2f}") # 图例标签
    
    plt.xlabel('Lag Distance (m)', fontweight='bold', fontsize=12) # 设置X轴标签
    plt.ylabel('Semivariance', fontweight='bold', fontsize=12) # 设置Y轴标签
    plt.legend(frameon=False, loc='lower right', fontsize=10) # 设置图例
    plt.grid(True, linestyle='--', color='gray', alpha=0.5) # 添加网格线
    plt.axvline(x=1098.22, color='gray', linestyle='--', alpha=0.5) # 添加垂直参考线
    plt.text(1015, plt.ylim()[1]*0.65, '1098.22m', rotation=90, fontsize=10) # 在参考线旁添加文本
    plt.tick_params(axis='both', which='major', labelsize=11) # 设置刻度标签字体大小
    plt.tight_layout() # 自动调整布局，防止标签重叠
    plt.savefig('semiva_plot_a.png', dpi=300, bbox_inches='tight') # 保存图形为PNG文件
    plt.show() # 显示图形

def generate_results_table(results):
    """
    生成并打印半变异函数模型拟合结果的表格。
    参数:
        results (list): 包含模型拟合结果的字典列表。
    """
    table = [] # 初始化表格数据列表
    headers = ['Model', 'Nugget', 'Range (m)', 'Sill', 'R²', 'SSE', 'Nu (Matern)'] # 定义表格标题
    
    for res in results: # 遍历每个模型结果
        row = [
            res['model'].capitalize(), # 模型名称（首字母大写）
            f"{res['nugget']:.2f}", # 块金值（保留两位小数）
            f"{res['range']:.2f}", # 范围值（保留两位小数）
            f"{res['sill']:.2f}", # 基台值（保留两位小数）
            f"{res['r_squared']:.2f}", # R²值（保留两位小数）
            f"{res['sse']:.2f}", # SSE值（保留两位小数）
            f"{res.get('nu', '-')}" # Matern模型的nu参数，如果不存在则显示'-'
        ]
        table.append(row) # 将行数据添加到表格列表中
    
    print("\nVariogram Model Comparison") # 打印表格标题
    print(tabulate(table, headers=headers, tablefmt="grid", stralign="center", numalign="center")) # 打印格式化后的表格

# 主程序入口
if __name__ == "__main__":
    # 加载数据 
    tif_path = r"E:\数据\v2.3\random_crops\stratified_sample.tif" # 分层采样后的重组成tif文件
    data, transform, nodata = load_geotiff(tif_path) # 加载GeoTIFF数据
    
    # 预处理数据
    rows, cols, values = preprocess_data(data, nodata) # 预处理数据，获取有效数据的行、列和值
    coordinates = calculate_coordinates(rows, cols, transform) # 计算有效数据的UTM坐标
    
    # 下采样数据
    coordinates, values = downsample_data(coordinates, values, sample_size=10000) # 对数据进行下采样
    
    # 清理内存
    del data, rows, cols # 删除不再需要的变量以释放内存
    gc.collect() # 执行垃圾回收
    
    # 拟合半变异函数模型
    results = fit_variogram_models(coordinates, values) # 拟合不同模型并获取结果
    
    # 生成结果表格
    generate_results_table(results) # 打印模型拟合结果表格
    
    # 绘制图形（传递样本数量N）
    plot_comparison(results, N=len(coordinates)) # 绘制模型比较图

#------------------------------------方向半变异函数-----------------------------------------------------
import rasterio # 导入rasterio库，用于处理栅格数据
import numpy as np # 导入numpy库，用于数值计算
import matplotlib.pyplot as plt # 导入matplotlib.pyplot库，用于绘图
from skgstat import Variogram, DirectionalVariogram # 从skgstat库导入Variogram和DirectionalVariogram类，用于半变异函数计算
import skgstat as skg # 导入skgstat库，用于地统计分析
from tabulate import tabulate # 导入tabulate库，用于美化表格输出
import gc # 导入gc库，用于垃圾回收，管理内存
import pandas as pd # 导入pandas库，用于数据处理，特别是Dunn's多重比较的结果
from scipy.stats import kruskal # 从scipy.stats导入kruskal函数，用于Kruskal-Wallis非参数检验
from scikit_posthocs import posthoc_dunn # 从scikit_posthocs导入posthoc_dunn函数，用于Dunn's多重比较
import seaborn as sns # 导入seaborn库，用于统计图形绘制

# 配置SCI论文绘图参数 - 进一步增大字体和调整线条
plt.rcParams.update({
    'figure.dpi': 300, # 设置图形分辨率为300 DPI
    'font.family': 'Times New Roman', # 设置字体为Times New Roman
    'font.size': 12,  # 全局基础字体增大
    'axes.labelsize': 14, # 轴标签字体增大
    'axes.labelweight': 'bold', # 轴标签字体加粗
    'xtick.labelsize': 12, # X轴刻度标签字体增大
    'ytick.labelsize': 12, # Y轴刻度标签字体增大
    'axes.spines.top': True, # 显示顶部轴脊
    'axes.spines.right': True # 显示右侧轴脊
})

def load_geotiff(path):
    """
    加载GeoTIFF文件并返回数据、仿射变换矩阵和NoData值。
    参数:
        path (str): GeoTIFF文件的路径。
    返回:
        tuple: 包含栅格数据 (numpy.ndarray), 仿射变换矩阵 (rasterio.Affine), 和NoData值。
    """
    with rasterio.open(path) as src: # 以只读模式打开GeoTIFF文件
        data = src.read(1) # 读取第一个波段的数据
        transform = src.transform # 获取仿射变换矩阵
        nodata = src.nodata # 获取NoData值
    return data, transform, nodata

def preprocess_data(data, nodata):
    """
    数据预处理：将NoData值替换为NaN，并提取有效数据的行、列索引和值。
    参数:
        data (numpy.ndarray): 栅格数据。
        nodata (int/float): 栅格数据的NoData值。
    返回:
        tuple: 包含有效数据的行索引 (numpy.ndarray), 列索引 (numpy.ndarray), 和值 (numpy.ndarray)。
    """
    data = np.where(data == nodata, np.nan, data) # 将NoData值替换为NaN
    mask = ~np.isnan(data) # 创建一个布尔掩膜，标记非NaN值
    rows, cols = np.where(mask) # 获取非NaN值的行和列索引
    values = data[rows, cols] # 提取非NaN值
    return rows, cols, values

def calculate_coordinates(rows, cols, transform):
    """
    根据行、列索引和仿射变换矩阵计算UTM坐标。
    参数:
        rows (numpy.ndarray): 数据的行索引。
        cols (numpy.ndarray): 数据的列索引。
        transform (rasterio.Affine): 栅格的仿射变换矩阵。
    返回:
        numpy.ndarray: 包含X和Y坐标的二维数组。
    """
    x_coords, y_coords = transform * (cols + 0.5, rows + 0.5) # 计算每个像素中心点的UTM坐标
    return np.column_stack((x_coords, y_coords)) # 将X和Y坐标堆叠成一个二维数组

def downsample_data(coordinates, values, sample_size=10000):
    """
    数据下采样：如果数据量超过指定样本大小，则随机抽取子集。
    参数:
        coordinates (numpy.ndarray): 坐标数据。
        values (numpy.ndarray): 值数据。
        sample_size (int): 目标样本大小。
    返回:
        tuple: 下采样后的坐标和值。
    """
    if len(coordinates) > sample_size: # 如果数据点数量大于样本大小
        np.random.seed(5) # 设置随机种子，以便结果可复现
        indices = np.random.choice(len(coordinates), sample_size, replace=False) # 随机选择指定数量的索引
        return coordinates[indices], values[indices] # 返回下采样后的坐标和值
    return coordinates, values # 如果数据量不大，则直接返回原始数据

def fit_variogram_models(coordinates, values):
    """
    拟合不同类型的半变异函数模型（球状、指数、高斯、Matern）。
    参数:
        coordinates (numpy.ndarray): 坐标数据。
        values (numpy.ndarray): 值数据。
    返回:
        list: 包含每个模型拟合结果的字典列表。
    """
    models = ['spherical', 'exponential', 'gaussian', 'matern'] # 定义要拟合的模型列表
    results = [] # 初始化存储结果的列表
    
    for model in models: # 遍历每个模型
        try:
            # 拟合半变异函数模型
            V = Variogram(
                coordinates=coordinates, # 坐标
                values=values, # 值
                n_lags=50, # 滞后距离的数量
                maxlag=0.8, # 最大滞后距离（相对于数据范围的百分比）
                model=model, # 半变异函数模型类型
                estimator='matheron', # 估计器类型
            )
            
            # 获取模型参数 (范围, 基台, 块金)
            params = V.parameters
            
            # 计算模型的预测值
            pred = V.model(V.bins, *params)
            
            # 计算R²和SSE (残差平方和)
            exp = V.experimental # 实验半变异函数值
            mask = ~np.isnan(exp) # 创建掩膜，排除NaN值
            exp = exp[mask] # 过滤掉NaN值
            pred = pred[mask] # 过滤掉NaN值
            ss_res = np.sum((exp - pred) ** 2)  # 计算残差平方和 (SSE)
            ss_tot = np.sum((exp - np.mean(exp)) ** 2) # 计算总平方和
            r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else np.nan # 计算R²
            
            # 收集结果
            result = {
                'model': model, # 模型名称
                'nugget': params[2], # 块金值
                'range': params[0], # 范围值
                'sill': params[1], # 基台值
                'r_squared': r_squared, # R²值
                'variogram': V, # 半变异函数对象
                'sse': ss_res # 残差平方和
            }
            
            # 处理Matern模型的额外参数 (nu)
            if model == 'matern' and len(params) > 3:
                result['nu'] = params[3] # Matern模型的nu参数
                
            results.append(result) # 将结果添加到列表中
            
        except Exception as e: # 捕获异常
            print(f"Error fitting {model} model: {str(e)}") # 打印错误信息
            continue # 继续下一个模型
    
    return results

def plot_comparison(results, N):
    """
    绘制不同半变异函数模型的拟合比较图。
    参数:
        results (list): 包含模型拟合结果的字典列表。
        N (int): 样本数量。
    """
    # 调整figsize以适应双栏中的双栏宽度（整页宽度），并给下半部分留出空间
    # 通常整页宽度在6.5到7英寸，这里设置为6.8，高度设置为5.5，以便与下半部分组合
    plt.figure(figsize=(6.8, 5.5)) # 设置图形大小
    
    # 绘制实验值
    V = results[0]['variogram'] # 获取第一个模型的半变异函数对象（用于获取实验值）
    plt.scatter(V.bins, V.experimental, 
                color='k', alpha=0.5, label='Experimental', zorder=5, s=20) # 绘制散点图，表示实验半变异函数值
    
    # 生成平滑曲线的X轴范围
    x_smooth = np.linspace(0, V.bins.max(), 300)
    
    # 定义绘图样式
    colors = ['red','blue','green','purple'] # 曲线颜色
    linestyles = ['-', '--', '-.', ':'] # 曲线线型
    
    # 绘制拟合曲线
    for i, res in enumerate(results): # 遍历每个模型结果
        V = res['variogram'] # 获取当前模型的半变异函数对象
        y_smooth = V.model(x_smooth, *V.parameters) # 计算模型在平滑X轴范围上的预测值
        plt.plot(x_smooth, y_smooth,
                 color=colors[i], # 颜色
                 linestyle=linestyles[i], # 线型
                 linewidth=2.5, # 线宽
                 label=f"{res['model'].capitalize()}\nR²={res['r_squared']:.2f}, SSE={res['sse']:.2f}") # 图例标签
    
    plt.xlabel('Lag Distance (m)', fontweight='bold', fontsize=12) # 设置X轴标签
    plt.ylabel('Semivariance', fontweight='bold', fontsize=12) # 设置Y轴标签
    plt.legend(frameon=False, loc='lower right', fontsize=10) # 设置图例
    plt.grid(True, linestyle='--', color='gray', alpha=0.5) # 添加网格线
    plt.axvline(x=1098.22, color='gray', linestyle='--', alpha=0.5) # 添加垂直参考线
    plt.text(1015, plt.ylim()[1]*0.65, '1098.22m', rotation=90, fontsize=10) # 在参考线旁添加文本
    plt.tick_params(axis='both', which='major', labelsize=11) # 设置刻度标签字体大小
    plt.tight_layout() # 自动调整布局，防止标签重叠
    plt.savefig('semiva_plot_a.png', dpi=300, bbox_inches='tight') # 保存图形为PNG文件
    plt.show() # 显示图形

def generate_results_table(results):
    """
    生成并打印半变异函数模型拟合结果的表格。
    参数:
        results (list): 包含模型拟合结果的字典列表。
    """
    table = [] # 初始化表格数据列表
    headers = ['Model', 'Nugget', 'Range (m)', 'Sill', 'R²', 'SSE', 'Nu (Matern)'] # 定义表格标题
    
    for res in results: # 遍历每个模型结果
        row = [
            res['model'].capitalize(), # 模型名称（首字母大写）
            f"{res['nugget']:.2f}", # 块金值（保留两位小数）
            f"{res['range']:.2f}", # 范围值（保留两位小数）
            f"{res['sill']:.2f}", # 基台值（保留两位小数）
            f"{res['r_squared']:.2f}", # R²值（保留两位小数）
            f"{res['sse']:.2f}", # SSE值（保留两位小数）
            f"{res.get('nu', '-')}" # Matern模型的nu参数，如果不存在则显示'-'
        ]
        table.append(row) # 将行数据添加到表格列表中
    
    print("\nVariogram Model Comparison") # 打印表格标题
    print(tabulate(table, headers=headers, tablefmt="grid", stralign="center", numalign="center")) # 打印格式化后的表格

# 主程序入口
if __name__ == "__main__":
    # 加载数据 - 
    tif_path = r"E:\数据\v2.3\random_crops\stratified_sample.tif" # 分层采样后的tif文件
    
    # 加载GeoTIFF数据
    data, transform, nodata = load_geotiff(tif_path)
    
    # 预处理数据
    rows, cols, values = preprocess_data(data, nodata)
    coordinates = calculate_coordinates(rows, cols, transform)
    
    # 下采样数据
    coordinates, values = downsample_data(coordinates, values, sample_size=10000)
    
    # 清理内存
    del data, rows, cols
    gc.collect()
    
    # 拟合半变异函数模型
    results = fit_variogram_models(coordinates, values)
    
    # 生成结果表格
    generate_results_table(results)
    
    # 绘制图形（传递样本数量N）
    plot_comparison(results, N=len(coordinates))

    # 计算四个方向的半变异函数
    # DirectionalVariogram用于计算指定方向上的半变异函数
    Vnorth = DirectionalVariogram(coordinates, values, azimuth=90, tolerance=10, maxlag=0.8, n_lags=20) # 正北方向（90度）
    Veast = DirectionalVariogram(coordinates, values, azimuth=0, tolerance=10, maxlag=0.8, n_lags=20) # 正东方向（0度）
    Vnortheast = DirectionalVariogram(coordinates, values, azimuth=45, tolerance=10, maxlag=0.8, n_lags=20) # 东北方向（45度）
    Vsoutheast = DirectionalVariogram(coordinates, values, azimuth=135, tolerance=10, maxlag=0.8, n_lags=20) # 东南方向（135度）

    # 提取1000m处的半变异函数值（或最接近的bin）
    def get_value_at_distance(variogram, target_distance=1000):
        """
        获取最接近目标距离处的半变异函数值。
        参数:
            variogram (skgstat.Variogram): 半变异函数对象。
            target_distance (float): 目标距离。
        返回:
            float: 最接近目标距离处的实验半变异函数值。
        """
        bins = variogram.bins # 滞后距离分箱
        values = variogram.experimental # 实验半变异函数值
        if not bins.size:  # 检查bins是否为空
            return np.nan # 如果为空则返回NaN
        idx = np.abs(bins - target_distance).argmin() # 找到最接近目标距离的索引
        return values[idx] # 返回对应的值

    north_value = get_value_at_distance(Vnorth) # 获取北方向1000m处的半变异函数值
    east_value = get_value_at_distance(Veast) # 获取东方向1000m处的半变异函数值
    northeast_value = get_value_at_distance(Vnortheast) # 获取东北方向1000m处的半变异函数值
    southeast_value = get_value_at_distance(Vsoutheast) # 获取东南方向1000m处的半变异函数值

    # 为每个方向创建多个样本点（通过块重采样）
    def create_directional_samples(coordinates, values, azimuth, tolerance, n_samples=30, sample_size=500):
        """
        创建方向特异性的多个样本点集，用于统计检验。
        通过随机抽取子集并计算其方向半变异函数在1000m处的值。
        参数:
            coordinates (numpy.ndarray): 原始坐标数据。
            values (numpy.ndarray): 原始值数据。
            azimuth (float): 方向（方位角）。
            tolerance (float): 方向容差。
            n_samples (int): 要生成的样本集数量。
            sample_size (int): 每个样本集的大小。
        返回:
            numpy.ndarray: 包含多个样本半变异函数值的数组。
        """
        samples = [] # 初始化样本列表
        for _ in range(n_samples): # 循环生成指定数量的样本
            # 随机选择一部分点进行重采样
            indices = np.random.choice(len(coordinates), sample_size, replace=False) # 随机选择索引
            subset_coords = coordinates[indices] # 提取子集坐标
            subset_values = values[indices] # 提取子集值
            
            # 计算方向半变异函数
            vario = DirectionalVariogram(subset_coords, subset_values, 
                                         azimuth=azimuth, tolerance=tolerance, 
                                         maxlag=0.8, n_lags=20)
            
            # 提取1000m处的值
            value = get_value_at_distance(vario) # 获取1000m处的半变异函数值
            if not np.isnan(value): # 如果值不是NaN
                samples.append(value) # 添加到样本列表
        
        return np.array(samples) # 返回样本数组

    # 为每个方向创建30个样本
    np.random.seed(42)  # 设置随机种子保证结果可重现
    north_samples = create_directional_samples(coordinates, values, 90, 10) # 北方向样本
    east_samples = create_directional_samples(coordinates, values, 0, 10) # 东方向样本
    northeast_samples = create_directional_samples(coordinates, values, 45, 10) # 东北方向样本
    southeast_samples = create_directional_samples(coordinates, values, 135, 10) # 东南方向样本

    # 执行Kruskal-Wallis检验，用于比较多个独立样本的中位数是否存在显著差异
    statistic, p_value = kruskal(north_samples, east_samples, northeast_samples, southeast_samples)

    # 如果Kruskal-Wallis检验结果显著 (p < 0.05)，则执行Dunn's多重比较
    if p_value < 0.05:
        # Dunn's多重比较，p_adjust='bonferroni'进行Bonferroni校正
        posthoc = posthoc_dunn([north_samples, east_samples, northeast_samples, southeast_samples], 
                               p_adjust='bonferroni')
        # 重命名Dunn's比较结果的列和索引，使其更具可读性
        posthoc.columns = ['North', 'East', 'Northeast', 'Southeast']
        posthoc.index = ['North', 'East', 'Northeast', 'Southeast']
    else:
        posthoc = None # 如果不显著，则不执行多重比较

    # 可视化四个方向的半变异函数
    # 调整figsize以适应双栏中的双栏宽度（整页宽度），并给上半部分留出空间
    plt.figure(figsize=(6.8, 5.5)) # 与上图保持相同宽度，调整高度

    # 绘制各个方向的实验半变异函数曲线
    plt.plot(Vnorth.bins, Vnorth.experimental, 'o-', color='red', label='North-South', markersize=6, linewidth=1.5) # 增大标记和线宽
    plt.plot(Veast.bins, Veast.experimental, 's-', color='blue', label='East-West', markersize=6, linewidth=1.5)
    plt.plot(Vnortheast.bins, Vnortheast.experimental, 'd-', color='green', label='Northeast', markersize=6, linewidth=1.5)
    plt.plot(Vsoutheast.bins, Vsoutheast.experimental, '^-', color='purple', label='Southeast', markersize=6, linewidth=1.5)

    # 添加1000m参考线
    plt.axvline(x=1000, color='gray', linestyle='--', alpha=0.5)
    plt.text(1010, plt.ylim()[1]*0.9, '1000m', rotation=90, fontsize=10) # 文本字体增大

    plt.xlabel('Lag Distance (m)', fontsize=12) # 轴标签字体增大
    plt.ylabel('Semivariance', fontsize=12) # 轴标签字体增大
    #plt.title('Directional Semivariograms') # 通常在组合图时不加子图标题，在总标题中说明
    plt.legend(frameon=False, facecolor='none', loc='lower right', fontsize=10) # 图例字体增大
    plt.grid(True, linestyle='--', alpha=0.7) # 添加网格线
    plt.tick_params(axis='both', which='major', labelsize=11) # 刻度标签字体增大
    plt.tight_layout() # 自动调整布局
    plt.savefig('semiva_plot_b1.png', dpi=300, bbox_inches='tight') # 保存为独立文件
    plt.show() # 显示图形

    # 可视化1000m处的差异 - 如果这个图不需要出现在最终组合图中，可以考虑移除
    plt.figure(figsize=(6.8, 5.0)) # 调整figsize，可能需要比上面两个图矮一点
    # 使用seaborn绘制箱线图，比较四个方向在1000m处半变异函数值的分布
    sns.boxplot(data=[north_samples, east_samples, northeast_samples, southeast_samples])
    plt.xticks([0, 1, 2, 3], ['North', 'East', 'Northeast', 'Southeast'], fontsize=12) # 设置X轴刻度标签
    plt.ylabel('Semivariance at 1000m', fontsize=12) # 设置Y轴标签
    plt.title(f'Comparison of Semivariance at 1000m\nKruskal-Wallis p-value: {p_value:.4f}', fontsize=14) # 设置标题，包含Kruskal-Wallis p值
    plt.grid(True, linestyle='--', alpha=0.7) # 添加网格线

    # 添加显著性标注
    if p_value < 0.05: # 如果p值小于0.05，表示显著
        plt.text(1.5, max([max(s) for s in [north_samples, east_samples, northeast_samples, southeast_samples]])*1.05,
                 f'Kruskal-Wallis检验: p={p_value:.4f}*', ha='center', fontsize=12) # 添加带星号的显著性文本
    else: # 如果p值不小于0.05，表示不显著
        plt.text(1.5, max([max(s) for s in [north_samples, east_samples, northeast_samples, southeast_samples]])*1.05,
                 f'Kruskal-Wallis检验: p={p_value:.4f}', ha='center', fontsize=12) # 添加不带星号的文本
    plt.tick_params(axis='both', which='major', labelsize=11) # 刻度标签字体增大
    plt.tight_layout() # 自动调整布局
    plt.savefig('semiva_plot_b2.png', dpi=300, bbox_inches='tight') # 保存为独立文件
    plt.show() # 显示图形

    # 输出统计检验结果
    print(f"Kruskal-Wallis检验统计量: {statistic:.4f}") # 打印Kruskal-Wallis统计量
    print(f"p值: {p_value:.4f}") # 打印p值

    if posthoc is not None: # 如果Dunn's多重比较结果存在
        print("\nDunn's多重比较结果（调整后p值）:")
        print(posthoc) # 打印Dunn's多重比较结果
        
        # 找出显著差异的方向对
        significant_pairs = [] # 初始化显著差异对列表
        directions = ['North', 'East', 'Northeast', 'Southeast'] # 定义方向名称
        for i in range(len(directions)): # 遍历所有方向对
            for j in range(i + 1, len(directions)):
                if posthoc.iloc[i, j] < 0.05: # 如果p值小于0.05，表示显著差异
                    significant_pairs.append(f"{directions[i]} vs {directions[j]}: p={posthoc.iloc[i, j]:.4f}") # 添加到显著差异对列表
        
        if significant_pairs: # 如果存在显著差异对
            print("\n显著差异的方向对:")
            for pair in significant_pairs: # 打印所有显著差异对
                print(pair)
        else:
            print("\n各方向对之间无显著差异") # 否则打印无显著差异
