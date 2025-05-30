# MD_tutorial_miyabi

# 1. Module installation

## CPU version
```sh
# conda がインストール済みであることを前提にしています

# 新規環境を作成し、Python や OpenMM をインストール
conda create -n openmm_cpu_env python=3.9
conda activate openmm_cpu_env

# CPU 版の OpenMM をインストール (conda-forge チャンネルを利用)
conda install -c conda-forge openmm
```

## GPU version
```sh
conda create -n openmm_gpu_env python=3.9
conda activate openmm_gpu_env
# OpenMM インストール（conda-forge）
conda install -c conda-forge openmm
# CUDA ランタイムも必要な場合
conda install -c conda-forge cudatoolkit
```

## checking 
```sh
$ conda info -e
...
openmm_cpu_env         /home/t01045/.conda/envs/openmm_cpu_env
openmm_gpu_env       * /home/t01045/.conda/envs/openmm_gpu_env
...
```


# 2. トイケース用の小さな RNA 構造を用意

2.1 PDB の用意

ここでは例として、PDB ID “1ZIH” (短い RNA ヘアピン) を使ってみます。
	1.	RCSB PDB などから 1ZIH をダウンロードし、1zih.pdb とリネーム
	2.	同じフォルダに配置

あるいは独自に RNA の 4-6 塩基程度の断片を用意した PDB を置いてもOKです。

## CPU version で実行
```sh
conda activate openmm_cpu_env
python run_min_md_cpu.py
```


# RNA 用の力場の installation
```sh
conda install -c conda-forge parmed
find ~/.conda/envs/openmm_cpu_env/ -name '*rna*'

```