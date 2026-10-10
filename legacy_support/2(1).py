import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

# 读取CSV文件（假设文件名为data.csv）
df = pd.read_csv('data.csv', header=None, names=['X', 'Y', 'Probability'])
df['Probability'] = pd.to_numeric(df['Probability'], errors='coerce')

# 创建图形和坐标轴
plt.figure(figsize=(6, 6), dpi=100)
scatter = plt.scatter(
    df['X'],
    df['Y'],
    c=df['Probability'],  # 颜色映射到概率值
    cmap='viridis',       # 使用viridis色图（黄-蓝渐变）
    norm=colors.LogNorm(),# 对数标准化（因概率值跨度大）
    s=100,               # 点的大小
    edgecolor='w',        # 白色边界
    linewidth=0.5,
    alpha=0.8
)



# 添加标题和标签
plt.title('Hit Probability Distribution', fontsize=14, pad=20)
plt.xlabel('X Coordinate', fontsize=12)
plt.ylabel('Y Coordinate', fontsize=12)

# 背景优化
plt.gca().set_facecolor('#f0f0f0')

# 自动调整布局
plt.tight_layout()
plt.savefig('heat_scatter.png', bbox_inches='tight')  # 保存图像
plt.show()