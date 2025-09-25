import numpy as np # 导入numpy库，用于数值计算
import pandas as pd # 导入pandas库，用于数据处理和分析
import rasterio # 导入rasterio库，用于处理栅格数据
from rasterio.mask import mask # 从rasterio.mask模块导入mask函数，用于按矢量边界裁剪栅格
import geopandas as gpd # 导入geopandas库，用于处理地理空间矢量数据
from shapely.geometry import Polygon, Point # 从shapely.geometry导入Polygon和Point类，用于创建几何对象
from shapely.ops import transform, unary_union # 从shapely.ops导入transform和unary_union，用于几何操作（transform用于坐标转换，unary_union用于合并几何对象）
from pyproj import Transformer # 导入pyproj的Transformer，用于坐标系转换（虽然代码中未直接使用，但通常用于地理坐标转换）
import math # 导入math库，用于数学运算，如三角函数和平方根
import os # 导入os库，用于操作系统相关功能，如路径操作和文件夹创建
from tqdm import tqdm # 导入tqdm库，用于显示循环进度条，提供可视化进度反馈

# 1. 读取LCZ栅格数据
lcz_path = r"E:\数据\v2.3\lcz_yrd2_30m补0.tif" # 定义LCZ（局部气候区）栅格数据的路径
with rasterio.open(lcz_path) as src: # 以只读模式打开LCZ栅格文件
    lcz_data = src.read(1) # 读取栅格的第一个波段数据（通常LCZ分类结果只有一个波段）
    lcz_meta = src.meta.copy() # 复制栅格的元数据，包括CRS、变换矩阵等
    lcz_crs = src.crs # 获取栅格的坐标参考系统（CRS）
    lcz_transform = src.transform # 获取栅格的仿射变换矩阵，用于将行列索引转换为地理坐标
    pixel_width = lcz_meta['transform'].a # 获取栅格像素的宽度（x方向分辨率）
    pixel_height = abs(lcz_meta['transform'].e) # 获取栅格像素的高度（y方向分辨率，取绝对值因为e通常为负）

# 2. 读取GWRR系数数据（使用UTF-8编码）
# 假设此CSV文件包含了每个监测点在4000m尺度下，各LCZ类型PLAND（景观百分比）的GWRR回归系数。
# encoding='gbk' 用于处理可能存在的中文编码问题，确保文件正确读取。
gwrr_data = pd.read_csv(r"E:\30m\景观格局指数\自适应\GWRR_4000m_PLAND_coefficients.csv", encoding='gbk')

# 3. 转换监测点坐标至EPSG:32651 (UTM WGS84 Zone 51N)
# 将GWRR系数数据转换为GeoDataFrame，以便进行地理空间操作。
gdf = gpd.GeoDataFrame(
    gwrr_data, # 基于GWRR系数数据创建GeoDataFrame
    geometry=[Point(lon, lat) for lon, lat in zip(gwrr_data['经度'], gwrr_data['纬度'])], # 从经度、纬度列创建Shapely Point几何对象
    crs="EPSG:4326" # 原始坐标系为WGS84地理坐标系
)
# 将GeoDataFrame的坐标系转换为UTM投影坐标系 (EPSG:32651)，以便进行距离和面积的精确计算。
gdf_utm = gdf.to_crs("EPSG:32651")

# 4. 定义8个方向的角度范围（以东为0°，逆时针增加）
# 每个方向由名称、起始角度和结束角度定义。
directions = [
    ("东", 337.5, 22.5),   # 东方向（0°±22.5°，即337.5°到22.5°，跨越0度）
    ("东北", 22.5, 67.5), # 东北方向
    ("北", 67.5, 112.5),  # 北方向
    ("西北", 112.5, 157.5), # 西北方向
    ("西", 157.5, 202.5),  # 西方向
    ("西南", 202.5, 247.5), # 西南方向
    ("南", 247.5, 292.5),  # 南方向
    ("东南", 292.5, 337.5)  # 东南方向
]

# 5. 创建扇形函数（修正角度转换逻辑）
def create_sector(center, radius, start_angle, end_angle):
    """
    生成扇形多边形。
    修正了角度转换逻辑，确保角度在0-360范围内，并正确处理跨0°的情况。
    参数:
        center (shapely.geometry.Point): 扇形中心点。
        radius (float): 扇形半径。
        start_angle (float): 扇形起始角度（度）。
        end_angle (float): 扇形结束角度（度）。
    返回:
        shapely.geometry.Polygon: 生成的扇形多边形。
    """
    # 确保角度在0-360范围内
    start_angle = start_angle % 360
    end_angle = end_angle % 360
    
    # 处理跨0°的情况（例如东方向：337.5-22.5）
    if start_angle > end_angle:
        # 分割为两部分：start_angle到360度 和 0度到end_angle
        points1 = []
        # 使用np.linspace生成角度序列，确保平滑的弧线，这里生成50个点
        for angle in np.linspace(start_angle, 360, 50):
            rad = math.radians(angle) # 将角度转换为弧度
            points1.append((center.x + radius * math.cos(rad), center.y + radius * math.sin(rad)))
        
        points2 = []
        for angle in np.linspace(0, end_angle, 50):
            rad = math.radians(angle)
            points2.append((center.x + radius * math.cos(rad), center.y + radius * math.sin(rad)))
        
        points = points1 + points2 + [(center.x, center.y)] # 合并两部分点并加上中心点以闭合多边形
        return Polygon(points) # 返回合并后的多边形
    
    # 普通情况（不跨0°）
    points = []
    # 使用np.linspace生成角度序列，确保平滑的弧线，这里生成100个点
    for angle in np.linspace(start_angle, end_angle, 100):
        rad = math.radians(angle) # 将角度转换为弧度
        x = center.x + radius * math.cos(rad)
        y = center.y + radius * math.sin(rad)
        points.append((x, y))
    points.append((center.x, center.y)) # 添加中心点以闭合多边形
    return Polygon(points) # 返回扇形多边形

# 6. 计算扇形内LCZ面积
def calculate_lcz_area(sector, lcz_data, transform):
    """
    计算扇形内各LCZ类型的精确面积。
    通过遍历扇形覆盖的像素，计算每个像素与扇形相交的面积。
    参数:
        sector (shapely.geometry.Polygon): 扇形多边形。
        lcz_data (numpy.ndarray): LCZ栅格数据。
        transform (rasterio.Affine): LCZ栅格的仿射变换矩阵。
    返回:
        dict: 包含各LCZ类型面积的字典。
    """
    lcz_area = {i: 0 for i in range(1, 18)} # 初始化LCZ类型面积字典（LCZ类型1-17，共17种）
    
    # 获取扇形边界框，并将其转换为栅格行列索引的范围
    minx, miny, maxx, maxy = sector.bounds
    row_min, col_min = rasterio.transform.rowcol(transform, minx, maxy) # 左上角行列
    row_max, col_max = rasterio.transform.rowcol(transform, maxx, miny) # 右下角行列
    
    # 遍历扇形边界框内的所有像素
    for row in range(row_min, row_max + 1):
        for col in range(col_min, col_max + 1):
            # 检查像素索引是否在栅格数据范围内，防止越界
            if 0 <= row < lcz_data.shape[0] and 0 <= col < lcz_data.shape[1]:
                # 获取像素中心点的UTM坐标
                x, y = rasterio.transform.xy(transform, row, col)
                # 创建代表当前像素的矩形多边形
                pixel = Polygon([
                    (x - pixel_width/2, y - pixel_height/2),
                    (x + pixel_width/2, y - pixel_height/2),
                    (x + pixel_width/2, y + pixel_height/2),
                    (x - pixel_width/2, y + pixel_height/2)
                ])
                
                # 计算当前像素与扇形的交集
                intersection = pixel.intersection(sector)
                if not intersection.is_empty: # 如果存在交集（像素部分或全部在扇形内）
                    lcz_type = lcz_data[row, col] # 获取当前像素的LCZ类型
                    if 1 <= lcz_type <= 17: # 确保LCZ类型在有效范围内
                        lcz_area[lcz_type] += intersection.area # 将交集面积累加到对应LCZ类型的总面积中
    
    return lcz_area

# 7. 计算方向影响量并调整扇形半径（修正版）
def calculate_directional_impact(station_row, center, base_radius, directions, lcz_data, transform):
    """
    计算各方向影响量并调整扇形半径，确保总面积守恒，允许半径扩展和缩小。
    此函数是自适应方向加权扇形影响区（ADWSIZ）的核心逻辑。
    参数:
        station_row (pd.Series): 当前监测点的GWRR系数数据行，包含各LCZ类型的PLAND系数。
        center (shapely.geometry.Point): 监测点中心坐标（UTM投影坐标）。
        base_radius (float): 基础圆形影响区的半径（例如4000m），用于计算原始扇形面积。
        directions (list): 包含方向名称和角度范围的列表。
        lcz_data (numpy.ndarray): LCZ栅格数据。
        transform (rasterio.Affine): LCZ栅格的仿射变换矩阵。
    返回:
        list: 包含调整后扇形信息的字典列表，每个字典包含方向名称、调整后的扇形几何、影响量、调整后半径等。
    """
    circle_area = math.pi * base_radius ** 2 # 计算基础圆形影响区的总面积
    direction_impacts = {} # 初始化存储各方向影响量及相关数据的字典
    
    # 计算各方向的原始影响量 (Id)
    for dir_name, start_angle, end_angle in directions:
        # 创建一个基础半径的扇形（用于计算该方向的原始LCZ组成）
        sector = create_sector(center, base_radius, start_angle, end_angle)
        # 计算该扇形内各LCZ类型的面积
        lcz_area = calculate_lcz_area(sector, lcz_data, transform)
        
        # 计算方向影响量 (Id)
        impact = 0
        for lcz_type in range(1, 18): # 遍历所有LCZ类型
            if lcz_area[lcz_type] > 0: # 如果该LCZ类型在扇形内有面积
                # 从GWRR系数数据中获取对应LCZ类型的PLAND系数
                # 假设GWRR_4000m_PLAND_coefficients.csv中列名为 'PLAND_LCZ类型_4000m'
                pland_coef = station_row[f'PLAND_{lcz_type}_4000m']
                # 累加影响量：系数 * (LCZ类型面积 / 基础圆形总面积)
                # 注意：这里除以circle_area是为了将LCZ面积转换为在整个圆形区域内的比例，再乘以系数
                impact += pland_coef * (lcz_area[lcz_type] / circle_area)
        
        direction_impacts[dir_name] = {
            'impact': impact, # 该方向的总影响量
            'lcz_area': lcz_area, # 该方向各LCZ类型的面积（用于后续统计）
            'sector': sector, # 原始扇形几何（用于后续调试或验证）
            'start_angle': start_angle, # 扇形起始角度
            'end_angle': end_angle # 扇形结束角度
        }
    
    # 计算总影响量（所有方向影响量绝对值之和）
    # 使用绝对值是为了避免正负抵消，确保影响大的方向能获得更大的调整。
    total_impact = sum(abs(dir_data['impact']) for dir_data in direction_impacts.values())
    
    # 防止除零错误：如果总影响量为0或负（理论上不应为负），则设置为一个极小正值
    if total_impact <= 0:
        total_impact = 1e-10 # 极小值，避免除以零
    
    # 计算各方向的归一化影响量（绝对值）
    # 归一化影响量表示该方向影响在所有方向总影响中的比例。
    norm_impacts = {
        dir_name: abs(dir_data['impact']) / total_impact
        for dir_name, dir_data in direction_impacts.items()
    }
    
    # 计算调整后扇形的目标面积
    # 目标面积是根据归一化影响量，在保持总面积不变的前提下，重新分配给每个方向的面积。
    target_areas = {
        dir_name: norm_impacts[dir_name] * circle_area
        for dir_name in direction_impacts.keys()
    }
    
    # 计算调整后的半径
    # 原始每个扇形的面积（圆形影响区被8个扇形均分）
    original_sector_area = circle_area / 8
    adjusted_sectors = [] # 初始化存储调整后扇形信息的列表

    # 遍历每个方向，计算调整后的半径并创建新的扇形
    for dir_name, dir_data in direction_impacts.items():
        target_area = target_areas[dir_name]
        
        # 关键修正：正确计算调整后半径，允许大于base_radius或小于base_radius
        # 根据面积比例调整半径：新半径 = 基础半径 * sqrt(目标面积 / 原始扇形面积)
        # 这样可以确保新扇形的面积与目标面积匹配，同时保持扇形形状（角度不变）。
        adjusted_radius = base_radius * math.sqrt(target_area / original_sector_area)
        
        # 创建调整后的扇形几何对象
        start_angle, end_angle = dir_data['start_angle'], dir_data['end_angle']
        adjusted_sector = create_sector(center, adjusted_radius, start_angle, end_angle)
        
        adjusted_sectors.append({
            'name': dir_name, # 方向名称
            'sector': adjusted_sector, # 调整后的扇形几何对象
            'impact': dir_data['impact'], # 原始影响量
            'radius': adjusted_radius, # 调整后的半径
            'target_area': target_area, # 目标面积
            'actual_area': adjusted_sector.area, # 实际调整后扇形的面积
            'lcz_area': dir_data['lcz_area'] # 原始扇形内各LCZ类型面积（用于后续统计）
        })
    
    # 验证调整后所有扇形的总面积是否与原始圆形面积近似
    total_adjusted_area = sum(s['actual_area'] for s in adjusted_sectors)
    print(f"原始圆形面积: {circle_area/1e6:.4f} km², 调整后总面积: {total_adjusted_area/1e6:.4f} km²")
    
    return adjusted_sectors

# 8. 裁剪并保存TIFF函数
def clip_and_save_tiff(merged_sector, station_code, lcz_path, output_dir):
    """
    裁剪合并后的扇形区域（自适应影响区）并保存为TIFF文件。
    参数:
        merged_sector (shapely.geometry.Polygon): 合并后的扇形多边形（代表自适应影响区）。
        station_code (str): 监测点编码，用于命名输出文件。
        lcz_path (str): LCZ栅格数据的路径。
        output_dir (str): 输出TIFF文件的目录。
    返回:
        str: 保存的TIFF文件路径，如果失败则返回None。
    """
    from shapely.geometry import mapping # 导入mapping函数，用于将shapely几何对象转换为GeoJSON格式，rasterio.mask需要此格式
    
    # 验证多边形有效性，如果无效则尝试修复（例如自相交问题）
    if not merged_sector.is_valid:
        merged_sector = merged_sector.buffer(0) # 使用buffer(0)尝试修复无效多边形，生成一个有效的副本
    
    sector_geojson = [mapping(merged_sector)] # 将合并后的扇形转换为GeoJSON格式列表，作为裁剪的掩膜
    
    with rasterio.open(lcz_path) as src: # 以只读模式打开LCZ栅格文件
        try:
            # 使用rasterio.mask裁剪栅格数据
            # crop=True 表示裁剪到几何对象的精确边界
            # all_touched=True 表示只要像素与几何对象有任何接触，就包含在裁剪结果中（而不是只包含像素中心在内的）
            out_image, out_transform = mask(src, sector_geojson, crop=True, all_touched=True)
            
            if out_image.size == 0: # 如果裁剪结果为空（没有像素数据）
                bounds = merged_sector.bounds # 获取扇形边界
                print(f"警告: {station_code}裁剪结果为空，边界: {bounds}") # 打印警告信息
                raise ValueError("裁剪结果为空") # 抛出ValueError，中断当前监测点的处理
                
            out_meta = src.meta.copy() # 复制原始栅格的元数据
            # 更新输出栅格的元数据，以反映裁剪后的影像尺寸和变换
            out_meta.update({
                'driver': 'GTiff', # 驱动格式为GeoTIFF
                'height': out_image.shape[1], # 输出影像的高度（行数）
                'width': out_image.shape[2], # 输出影像的宽度（列数）
                'transform': out_transform, # 输出影像的仿射变换矩阵
                'crs': src.crs # 输出影像的坐标参考系统（与原始栅格相同）
            })
            
            output_path = os.path.join(output_dir, f"{station_code}_扇形.tif") # 构造输出文件路径
            with rasterio.open(output_path, "w", **out_meta) as dst: # 以写入模式打开输出文件
                dst.write(out_image) # 写入裁剪后的影像数据
            
            return output_path # 返回保存的TIFF文件路径
            
        except Exception as e: # 捕获裁剪过程中的其他异常
            print(f"裁剪失败 - 监测点: {station_code}, 错误: {str(e)}") # 打印错误信息
            return None # 返回None表示裁剪失败

# 9. 主处理流程（使用UTF-8编码保存文件）
output_dir = r"E:\30m\景观格局指数\扇形TIFF结果_调整" # 定义输出目录，用于保存结果文件
os.makedirs(output_dir, exist_ok=True) # 创建输出目录，如果目录已存在则不会报错

processed_count = 0 # 初始化已成功处理的监测点计数
# 遍历每个监测点（使用tqdm显示进度条，提供可视化反馈）
for i, row in tqdm(gdf_utm.iterrows(), total=len(gdf_utm), desc="处理监测点"):
    station_code = row['监测点'] # 获取当前监测点的编码
    center = row.geometry # 获取当前监测点的几何对象（UTM坐标）
    base_radius = 4000 # 定义基础圆形半径为4000米，这是论文中确定的最佳过程尺度
    
    try:
        # 调用函数计算各方向影响量并调整扇形半径，生成自适应影响区
        adjusted_sectors = calculate_directional_impact(
            row, center, base_radius, directions, lcz_data, lcz_transform
        )
        
        # 输出东方向的调试信息，方便检查每个监测点的处理情况
        east_sector = next(s for s in adjusted_sectors if s['name'] == '东') # 找到东方向的扇形数据
        print(f"监测点 {station_code}: 东方向影响量={east_sector['impact']:.6f}, 半径={east_sector['radius']:.2f}m")
        
        # 保存详细统计信息到CSV文件
        stats = {
            '监测点': station_code, # 监测点编码
            '原始圆形面积(km²)': math.pi * base_radius ** 2 / 1e6, # 原始圆形影响区的面积（平方公里）
            '总影响量': sum(abs(s['impact']) for s in adjusted_sectors) # 所有方向影响量的绝对值之和
        }
        
        # 记录每个方向的详细信息
        for s in adjusted_sectors:
            stats[f"{s['name']}_影响量"] = s['impact'] # 各方向影响量
            stats[f"{s['name']}_归一化影响量"] = abs(s['impact']) / max(stats['总影响量'], 1e-10) # 归一化影响量（防止除以零）
            stats[f"{s['name']}_半径(m)"] = s['radius'] # 调整后的半径
            stats[f"{s['name']}_面积(km²)"] = s['actual_area'] / 1e6 # 实际调整后扇形的面积（平方公里）
            
            # 记录各方向的LCZ组成
            for lcz_type in range(1, 18):
                lcz_area = s['lcz_area'].get(lcz_type, 0) # 获取LCZ类型面积，如果不存在则为0
                stats[f"{s['name']}_LCZ{lcz_type}_面积(km²)"] = lcz_area / 1e6 # LCZ类型面积（平方公里）
                stats[f"{s['name']}_LCZ{lcz_type}_占比"] = lcz_area / s['actual_area'] if s['actual_area'] > 0 else 0 # LCZ类型占比（防止除以零）
        
        # 计算并记录整个合并扇形区域的LCZ组成
        for lcz_type in range(1, 18):
            lcz_area = 0
            for s in adjusted_sectors:
                lcz_area += s['lcz_area'].get(lcz_type, 0) # 累加所有调整后扇形中该LCZ类型的面积
            stats[f"LCZ{lcz_type}_面积(km²)"] = lcz_area / 1e6 # 总LCZ类型面积（平方公里）
            stats[f"LCZ{lcz_type}_占比"] = lcz_area / sum(s['actual_area'] for s in adjusted_sectors) if sum(s['actual_area'] for s in adjusted_sectors) > 0 else 0 # 总LCZ类型占比（防止除以零）
        
        # 合并所有调整后的扇形为一个多边形（unary_union用于合并重叠或相邻的几何对象）
        merged_sector = unary_union([s['sector'] for s in adjusted_sectors])
        stats['合并后总面积(km²)'] = merged_sector.area / 1e6 # 合并后总面积（平方公里）
        
        # 保存统计信息到CSV文件（使用UTF-8编码）
        pd.DataFrame([stats]).to_csv(
            os.path.join(output_dir, f"{station_code}_扇形统计.csv"),
            index=False, # 不写入行索引
            encoding='utf-8' # 指定UTF-8编码，确保中文兼容性
        )
        
        # 保存影响量统计到CSV文件（使用UTF-8编码）
        impact_stats = {f"{s['name']}_impact": s['impact'] for s in adjusted_sectors}
        impact_stats['监测点'] = station_code
        pd.DataFrame([impact_stats]).to_csv(
            os.path.join(output_dir, f"{station_code}_影响量.csv"),
            index=False, # 不写入行索引
            encoding='utf-8' # 指定UTF-8编码
        )
        
        # 裁剪并保存TIFF文件
        output_path = clip_and_save_tiff(merged_sector, station_code, lcz_path, output_dir)
        
        if output_path: # 如果TIFF文件成功保存
            processed_count += 1 # 增加成功处理的监测点计数
    
    except Exception as e: # 捕获处理单个监测点时的任何异常
        print(f"处理监测点{station_code}时出错: {str(e)}") # 打印错误信息

print(f"全部处理完成！共生成{processed_count}个TIFF文件，所有CSV文件以UTF-8编码保存")


