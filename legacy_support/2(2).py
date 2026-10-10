import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# 读取CSV文件
df = pd.read_csv('data1.csv', header=None, names=['X', 'Y','V','Probability'])
df['Probability'] = pd.to_numeric(df['Probability'], errors='coerce')


# 检查并处理无效值
if df['Probability'].isnull().any():
    print(f"警告：发现 {df['Probability'].isnull().sum()} 个无效的概率值")
    df = df.dropna(subset=['Probability'])  # 移除无效值

# 创建3D图形
fig = plt.figure(figsize=(14, 10), dpi=100)
ax = fig.add_subplot(111, projection='3d')

# 创建对数归一化对象
norm = colors.LogNorm(vmin=df['Probability'].min(), vmax=df['Probability'].max())

# 创建3D散点图
scatter = ax.scatter(
    df['X'],
    df['Y'],
    df['V'],
    c=df['Probability'],
    cmap='plasma',  # 使用等离子色图
    norm=norm,      # 对数标准化
    s=80,           # 点的大小
    alpha=0.8,      # 透明度
    edgecolor='w',  # 白色边框
    linewidth=0.5
)



# 设置轴标签
ax.set_xlabel('X Coordinate', fontsize=12, labelpad=15)
ax.set_ylabel('Y Coordinate', fontsize=12, labelpad=15)
ax.set_zlabel('V Coordinate', fontsize=12, labelpad=15)
# 在所有坐标轴设置之后添加以下代码
ax.xaxis.set_tick_params(labelsize=10)  # X轴刻度标签大小
ax.yaxis.set_tick_params(labelsize=10)  # Y轴刻度标签大小
ax.zaxis.set_tick_params(labelsize=10)  # Z轴刻度标签大小
# 设置标题
ax.set_title('3D Hit Probability Distribution', fontsize=16, pad=20)

# 优化视角
ax.view_init(elev=9, azim=-2)  # 设置初始视角

# 设置背景颜色
ax.xaxis.pane.set_edgecolor('w')
ax.yaxis.pane.set_edgecolor('w')
ax.zaxis.pane.set_edgecolor('w')
ax.set_facecolor('#f8f8f8')



# 调整布局
plt.tight_layout()
plt.savefig('3d_probability_heatmap.png', dpi=150, bbox_inches='tight')
plt.show()
#