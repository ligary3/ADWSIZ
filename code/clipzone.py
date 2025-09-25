# 裁剪小栅格
import rasterio # 导入rasterio库，用于处理栅格数据
from rasterio.mask import mask # 从rasterio.mask模块导入mask函数，用于按矢量边界裁剪栅格
import geopandas as gpd # 导入geopandas库，用于处理地理空间矢量数据
import shapely.geometry as sg # 导入shapely.geometry模块，用于几何对象的创建和操作
import os # 导入os库，用于操作系统相关功能，如路径操作和文件夹创建

# 定义输入文件路径
shp_path = r"E:\数据\v2\lcz\站点shp2024年均2.shp" # 监测站点shp文件的路径
landuse_data = r"E:\数据\v2.3\lcz_yrd2_30m.tif" # 土地利用（LCZ）栅格数据的路径

# 确保输出文件夹存在
output_folder = r"E:\30m\yuan1-62" # 裁剪结果的输出文件夹路径
if not os.path.exists(output_folder): # 如果输出文件夹不存在
    os.makedirs(output_folder) # 则创建该文件夹

# 读取矢量数据
gdf = gpd.read_file(shp_path) # 使用geopandas读取shp文件到GeoDataFrame
count = 0 # 初始化裁剪完成文件的计数器

# 读取栅格数据
with rasterio.open(landuse_data) as src: # 以只读模式打开土地利用栅格数据
    for index, row in gdf.iterrows(): # 遍历GeoDataFrame中的每一行（即每个监测站点）
        site = row['监测点'] # 获取当前监测点的名称
        geometry = row.geometry # 获取当前监测点的几何对象（通常是点）
        # 遍历不同的缓冲区半径
        for radius in [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000]:
            # 生成缓冲区
            buffer_geometry = geometry.buffer(radius) # 以当前监测点为中心，生成指定半径的圆形缓冲区

            # 定义裁剪后的土地利用数据输出路径和文件名
            landuse_clip = os.path.join(output_folder, f"{site}_{radius}m.tif")
            try:
                # 裁剪土地利用数据
                # mask函数用于根据矢量几何对象裁剪栅格数据
                # [sg.mapping(buffer_geometry)] 将shapely几何对象转换为rasterio可接受的格式
                # crop=True 表示裁剪到几何对象的精确边界
                out_image, out_transform = mask(src, [sg.mapping(buffer_geometry)], crop=True)
                out_meta = src.meta.copy() # 复制原始栅格的元数据
                # 更新元数据，包括驱动、高度、宽度和变换信息（新的仿射变换矩阵）
                out_meta.update({"driver": "GTiff",
                                 "height": out_image.shape[1], # 裁剪后影像的高度
                                 "width": out_image.shape[2], # 裁剪后影像的宽度
                                 "transform": out_transform}) # 裁剪后影像的仿射变换矩阵
                # 将裁剪后的影像写入新的GeoTIFF文件
                with rasterio.open(landuse_clip, "w", **out_meta) as dest:
                    dest.write(out_image) # 写入裁剪后的影像数据
                    count += 1 # 裁剪成功计数加1
                    print(f"裁剪 {site}_{radius}m.tif 完成。", count) # 打印裁剪完成信息

            except ValueError: # 捕获ValueError异常，通常发生在裁剪区域没有有效数据时
                print(f"裁剪 {site}_{radius}m.tif 时未找到有效数据。") # 打印未找到有效数据的提示

