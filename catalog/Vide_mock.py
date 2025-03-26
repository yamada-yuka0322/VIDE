import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
    
    # データ部分を読み込む
df = pd.read_csv("/mnt/data/yuka/dataset/cola/fof00100a.txt", sep=r"\s+", header=None)

header_lines = [
    "300.0",
    "0.273",
    "0.705",
    "0.73",
    str(len(df))
    ]

# インデックスを 0,1,2,... に変更
df.iloc[:, 0] = range(len(df))

# 4列目（3番目のインデックス）をすべて 1.0 にする
df.iloc[:, 4] = 1.0

# 5,6列目を削除
df = df.drop(columns=[5,6])

# ファイルに書き出す
with open("/mnt/data/yuka/dataset/vide_data/cola_test.dat", "w") as f:
    for line in header_lines:
        f.write(line + "\n")  # ヘッダー部分を縦に書き出し
    
    df.to_csv(f, sep=" ", index=False, header=False)