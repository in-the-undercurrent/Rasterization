import pandas as pd
import numpy as np
from numpy.f2py.rules import options
from shapely.geometry import Polygon, box
import os
import h5py
import time
from osgeo import gdal, osr

import multiprocessing
import warnings


# 忽略 FutureWarning

def is_crossing_prime_meridian(row):
    longitudes = [row['LongitudeA'], row['LongitudeB'], row['LongitudeC'], row['LongitudeD']]
    min_lon = min(longitudes)
    max_lon = max(longitudes)
    return (min_lon < 30 and max_lon >= 330)


# 对经度进行对称变换
def vector_mirror_longitude(row):
    # 从 row 中提取经纬度
    lons = [row['LongitudeA'], row['LongitudeB'], row['LongitudeC'], row['LongitudeD']]

    # 创建一个新的列表来存储修改后的经纬度
    mirrored_lons = []
    for lon in lons:
        # 如果经度在 180° 到 360° 范围内，则对称到负值
        if 180 <= lon < 360:
            mirrored_lon = lon - 360  # 变换到 -180 到 -1 范围内
        else:
            mirrored_lon = lon  # 其他经度不变

        mirrored_lons.append(mirrored_lon)
    # 更新 row 中的经度
    row['LongitudeA'] = mirrored_lons[0]
    row['LongitudeB'] = mirrored_lons[1]
    row['LongitudeC'] = mirrored_lons[2]
    row['LongitudeD'] = mirrored_lons[3]
    return row


# global_crossing_count = 0
def process_single_row(row):
    # global global_crossing_count
    grid_inter_split = pd.DataFrame(columns=['x', 'y', 'no2', 'percent'])

    if is_crossing_prime_meridian(row):
        # 打印原始坐标
        # print(row[['LongitudeA', 'LongitudeB', 'LongitudeC', 'LongitudeD']])

        # 处理跨越子午线的情况
        row = vector_mirror_longitude(row)

        # global_crossing_count +=1

        x_left_down = min(row['LatitudeA'], row['LatitudeB'], row['LatitudeC'], row['LatitudeD'])
        y_left_down = min(row['LongitudeA'], row['LongitudeB'], row['LongitudeC'], row['LongitudeD'])
        x_right_up = max(row['LatitudeA'], row['LatitudeB'], row['LatitudeC'], row['LatitudeD'])
        y_right_up = max(row['LongitudeA'], row['LongitudeB'], row['LongitudeC'], row['LongitudeD'])

        grid_range_x_left_down = int((x_left_down + 90) // 0.25)
        grid_range_y_left_down = int((y_left_down + 180) // 0.25)
        grid_range_x_right_up = int(((x_right_up + 90) // 0.25) + 1)
        grid_range_y_right_up = int(((y_right_up + 180) // 0.25) + 1)
        if (y_right_up-y_left_down)>30:
            print("面状矢量经度跨度过大")
            print(y_left_down,y_right_up)
            return grid_inter_split
        poly = Polygon([
            (row['LatitudeA'], row['LongitudeA']),
            (row['LatitudeB'], row['LongitudeB']),
            (row['LatitudeC'], row['LongitudeC']),
            (row['LatitudeD'], row['LongitudeD'])
        ])
        for grid_x in range(grid_range_x_left_down, grid_range_x_right_up):
            for grid_y in range(grid_range_y_left_down, grid_range_y_right_up):
                grid_cell = box(-90 + 0.25 * grid_x, -180 + 0.25 * grid_y, -90 + 0.25 * grid_x + 0.25,
                                -180 + 0.25 * grid_y + 0.25)
                if grid_cell.intersects(poly):
                    intersection = grid_cell.intersection(poly.buffer(0))
                    inter_area = intersection.area
                    new_row = pd.Series([grid_x, grid_y, row['no2_concentration'], inter_area / grid_cell.area],
                                        index=grid_inter_split.columns)
                    if new_row['percent'] < 0.001:
                        continue
                    if 0 <= new_row['y'] <= 719:
                        new_row['y'] = 720 + new_row['y']
                    else:
                        new_row['y'] -= 720
                    grid_inter_split = pd.concat([grid_inter_split, new_row.to_frame().T], ignore_index=True)


    else:
        # 不跨越子午线的情况
        x_left_down = min(row['LatitudeA'], row['LatitudeB'], row['LatitudeC'], row['LatitudeD'])
        y_left_down = min(row['LongitudeA'], row['LongitudeB'], row['LongitudeC'], row['LongitudeD'])
        x_right_up = max(row['LatitudeA'], row['LatitudeB'], row['LatitudeC'], row['LatitudeD'])
        y_right_up = max(row['LongitudeA'], row['LongitudeB'], row['LongitudeC'], row['LongitudeD'])

        grid_range_x_left_down = int((x_left_down + 90) // 0.25)
        grid_range_y_left_down = int(y_left_down // 0.25)
        grid_range_x_right_up = int(((x_right_up + 90) // 0.25) + 1)
        grid_range_y_right_up = int((y_right_up // 0.25) + 1)

        poly = Polygon([
            (row['LatitudeA'], row['LongitudeA']),
            (row['LatitudeB'], row['LongitudeB']),
            (row['LatitudeC'], row['LongitudeC']),
            (row['LatitudeD'], row['LongitudeD'])
        ])

        for grid_x in range(grid_range_x_left_down, grid_range_x_right_up):
            for grid_y in range(grid_range_y_left_down, grid_range_y_right_up):
                grid_cell = box(-90 + 0.25 * grid_x, 0.25 * grid_y, -90 + 0.25 * grid_x + 0.25, 0.25 * grid_y + 0.25)
                if grid_cell.intersects(poly):
                    intersection = grid_cell.intersection(poly.buffer(0))
                    inter_area = intersection.area
                    new_row = pd.Series([grid_x, grid_y, row['no2_concentration'], inter_area / grid_cell.area],
                                        index=grid_inter_split.columns)
                    if new_row['percent'] < 0.001:
                        continue
                    grid_inter_split = pd.concat([grid_inter_split, new_row.to_frame().T], ignore_index=True)

    return grid_inter_split


def process_hdf5_file(hdf5_path):
    # 打开HDF5文件
    with h5py.File(hdf5_path, "r") as hdf5_file:
        # 读取NO2浓度数据
        no2_data = hdf5_file['TOTAL_COLUMNS/NO2Tropo'][:]
        # 将负值替换为空值（NaN）
        no2_data = np.where(no2_data < 0, np.nan, no2_data)
        # 单位换算因子
        molecules_to_DU = 2.6867e16
        # 计算NO2浓度转换为DU
        no2_DU = no2_data / molecules_to_DU

        # 检测空值
        mask = np.isnan(no2_DU)
        # 获取非空值的索引
        valid_indices = np.where(~mask)[0]

        # 读取地理坐标数据
        latitudes = hdf5_file['GEOLOCATION/LatitudeCentre'][:]
        latitudeA = hdf5_file['GEOLOCATION/LatitudeA'][:]
        latitudeB = hdf5_file['GEOLOCATION/LatitudeB'][:]
        latitudeC = hdf5_file['GEOLOCATION/LatitudeC'][:]
        latitudeD = hdf5_file['GEOLOCATION/LatitudeD'][:]
        longitudes = hdf5_file['GEOLOCATION/LongitudeCentre'][:]
        longitudeA = hdf5_file['GEOLOCATION/LongitudeA'][:]
        longitudeB = hdf5_file['GEOLOCATION/LongitudeB'][:]
        longitudeC = hdf5_file['GEOLOCATION/LongitudeC'][:]
        longitudeD = hdf5_file['GEOLOCATION/LongitudeD'][:]

        # 根据有效索引筛选数据
        no2_DU = no2_DU[valid_indices]
        latitudes = latitudes[valid_indices]
        latitudeA = latitudeA[valid_indices]
        latitudeB = latitudeB[valid_indices]
        latitudeC = latitudeC[valid_indices]
        latitudeD = latitudeD[valid_indices]
        longitudes = longitudes[valid_indices]
        longitudeA = longitudeA[valid_indices]
        longitudeB = longitudeB[valid_indices]
        longitudeC = longitudeC[valid_indices]
        longitudeD = longitudeD[valid_indices]

    # 创建DataFrame存储数据
    data = {
        'no2_concentration': no2_DU.flatten(),
        'Latitude': latitudes.flatten(),
        'LatitudeA': latitudeA.flatten(),
        'LatitudeB': latitudeB.flatten(),
        'LatitudeC': latitudeC.flatten(),
        'LatitudeD': latitudeD.flatten(),
        'Longitude': longitudes.flatten(),
        'LongitudeA': longitudeA.flatten(),
        'LongitudeB': longitudeB.flatten(),
        'LongitudeC': longitudeC.flatten(),
        'LongitudeD': longitudeD.flatten(),
    }
    df = pd.DataFrame(data)

    # 创建一个空的 DataFrame，只包含列名
    grid_inter = pd.DataFrame(columns=['x', 'y', 'no2', 'percent'])
    # 初始化总的 DataFrame 来存储所有行的处理结果

    # global_crossing_count
    # 遍历每一行并调用 process_single_row 函数

    with multiprocessing.Pool(processes=16) as pool:
        # 使用pool.map并行处理每一行
        results = pool.map(process_single_row, [row for _, row in df.iterrows()])

    # 合并所有子进程的处理结果
    grid_inter_total = pd.concat(results, ignore_index=True)


    return grid_inter_total


def process_value(value, dup):
    # 筛选出具有相同 grid_id 的数据
    # 从 `dup` DataFrame 中筛选出 `grid_id` 等于当前 `value` 的子集
    subset = dup[dup['grid_id'] == value]

    # 计算 percent 列的和
    # 用于归一化每个 `percent` 权重以便后续计算加权平均值
    total_percent = subset['percent'].sum()

    # 初始化加权值
    weight_value = 0

    # 计算每一行的加权值
    # 对于每一行，根据 `percent` 和 `no2` 计算加权值，并归一化到 `total_percent` 上
    for index, row in subset.iterrows():
        weight_value += row['percent'] * row['no2'] / total_percent

    # 拆分 grid_id 得到 x 和 y
    # `grid_id` 是用逗号分隔的值（如 "x,y"），将其拆解为 x 和 y
    xy = value.split(',')

    # 创建一个新的 Pandas Series 来表示新的一行数据
    # 包含 x 坐标、y 坐标、加权计算后的 NO2 值、固定值 '10' 和原始的 grid_id
    # 返回这个新的数据行，以便后续在主进程中汇总
    return pd.Series([xy[0], xy[1], weight_value, '10', value], index=dup.columns)


def remove_and_weight_duplicates(df):
    # 创建 grid_id 列
    df['grid_id'] = df['x'].astype(str) + ',' + df['y'].astype(str)

    # 遮罩出重复的行
    dup_mask = df.duplicated(subset=['grid_id'], keep=False)
    dup = df[dup_mask]
    df_unique = df[~dup_mask]

    # 获取 'grid_id' 的唯一值
    unique_values = dup['grid_id'].unique()

    with multiprocessing.Pool(processes=16) as pool:
        # 使用 starmap 来传递多个参数
        new_rows = pool.starmap(process_value, [(value, dup) for value in unique_values])
    df_unique = df_unique._append(new_rows, ignore_index=True)
    df_unique = df_unique.drop(columns=['percent', 'grid_id'])
    df_unique['x'] = -90 + 0.25 * df_unique['x'].astype(float)
    df_unique['y'] = 0.25 * df_unique['y'].astype(float)
    return df_unique


def create_geotiff_from_dataframe(data, output_tiff):
    lat_min, lat_max = -90, 90  # 纬度范围
    # 修改经度范围
    lon_min, lon_max = -180, 180  # 改为 -180 到 180 的范围
    pixel_size = 0.25

    # 将输入数据的经度转换为 -180 到 180
    data['y'] = data['y'].apply(lambda lon: lon - 360 if lon >= 180 else lon)

    # 计算栅格的行列数
    nlat = int((lat_max - lat_min) / pixel_size)
    nlon = int((lon_max - lon_min) / pixel_size)

    # 创建一个空的栅格数据
    grid = np.full((nlat, nlon), np.nan)

    # 填充栅格数据
    for index, row in data.iterrows():
        lat_idx = int((row['x'] - lat_min) / pixel_size)
        lon_idx = int((row['y'] - lon_min) / pixel_size)
        if 0 <= lat_idx < nlat and 0 <= lon_idx < nlon:
            grid[lat_idx, lon_idx] = row['no2']

    # 在写入之前颠倒 grid 的纬度方向
    grid = np.flipud(grid)

    # 创建GeoTIFF文件
    driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(output_tiff, nlon, nlat, 1, gdal.GDT_Float32,options=['COMPRESS=LZW'])

    # 设置地理变换（左上角x, 像素宽度, 0, 左上角y, 0, 像素高度）
    dataset.SetGeoTransform((lon_min, pixel_size, 0, lat_max, 0, -pixel_size))

    # 设置投影信息
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)  # EPSG 4326 对应于 WGS 84
    dataset.SetProjection(srs.ExportToWkt())

    # 写入数据
    band = dataset.GetRasterBand(1)
    band.WriteArray(grid)

    # 设置 NoData 值
    band.SetNoDataValue(np.nan)

    # 关闭文件
    dataset.FlushCache()


warnings.simplefilter(action='ignore', category=FutureWarning   )
import os
import time
import pandas as pd

# 假设 process_hdf5_file、remove_and_weight_duplicates、create_geotiff_from_dataframe 已经定义
# 不再重复定义这些函数

if __name__ == '__main__':
    start_time = time.time()
    print(f"程序开始时间: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}")
    print(f"代码版本2.3")
    base_folder_path = "/home/yourname/project/lv2/"  # 根目录路径
    output_folder_path = "/home/yourname/project/lv3/A"  # 输出文件夹路径

    year = '2012'
    year_path = os.path.join(base_folder_path, year)

    # 检查该目录是否存在
    if os.path.isdir(year_path):
        for month in os.listdir(year_path):
            month_path = os.path.join(year_path, month)
            if not os.path.isdir(month_path):
                continue  # 跳过非文件夹的内容

            for day in os.listdir(month_path):
                day_path = os.path.join(month_path, day)
                if not os.path.isdir(day_path):
                    continue  # 跳过非文件夹的内容

                # 创建每日总表
                grid_total = pd.DataFrame(columns=['x', 'y', 'no2', 'percent'])

                # 遍历指定文件夹中的所有 HDF5 文件
                for file_name in os.listdir(day_path):
                    file_path = os.path.join(day_path, file_name)

                    # 检查文件是否是 HDF5 文件
                    if file_name.lower().endswith('.hdf5'):
                        try:
                            processed_df = process_hdf5_file(file_path)
                        except Exception as e:
                            # 打印错误的序号和详细错误信息
                            print(f"处理文件出错，文件: {file_name}, 错误: {e}")
                            continue
                        grid_total = pd.concat([grid_total, processed_df], ignore_index=True)
                # 去重和加权处理
                grid_final = remove_and_weight_duplicates(grid_total)
                #print(f"权重结束: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))}")
                # 生成输出 TIFF 文件路径
                date_str = f"{year}{month}{day}"
                out_tiff_path = os.path.join(output_folder_path, f"{year}",f"{month}",f"{date_str}.tif")
                os.makedirs(os.path.dirname(out_tiff_path), exist_ok=True)
                # 创建 GeoTIFF 文件
                create_geotiff_from_dataframe(grid_final, out_tiff_path)
                #print(f"写栅格结束: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))}")
                del grid_total, grid_final


    print(f"程序结束时间: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))}")
