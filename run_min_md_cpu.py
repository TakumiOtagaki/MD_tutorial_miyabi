#!/usr/bin/env python

import sys
from openmm import app
import openmm as mm
from openmm import unit

# -- 入力PDBファイル --
id = '1zih'
pdb_file = f'data/{id}.pdb'

# PDB読み込み
pdb = app.PDBFile(pdb_file)

# 力場読み込み
# もし 'amber14/rna.ff14SB.xml' が含まれていない場合は 'amber14/rna.xml' を使用。
# 実際に環境に存在するファイル名を必ず確認してください。
# forcefield = app.ForceField('amber14-all.xml', 'amber14/rna.xml')
forcefield = app.ForceField('amber14-all.xml', 'amber14/RNA.OL3.xml')


# システムを構築
# 本来は水中を想定しますが、ここでは真空(NoCutoff)設定に
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=None,
    removeCMMotion=True
)

# 温度やタイムステップなどを設定 (Langevin dynamics)
temperature = 300 * unit.kelvin
friction = 1.0 / unit.picosecond
timestep = 0.002 * unit.picoseconds

integrator = mm.LangevinIntegrator(temperature, friction, timestep)

# CPU を指定
platform = mm.Platform.getPlatformByName('CPU')

# Simulation オブジェクトを生成
simulation = app.Simulation(pdb.topology, system, integrator, platform)

# 初期座標セット
simulation.context.setPositions(pdb.positions)

# エネルギー最小化
print("Energy Minimization...")
simulation.minimizeEnergy(maxIterations=500)

# 最小化した構造を保存
min_positions = simulation.context.getState(getPositions=True).getPositions()
with open('minimized_cpu.pdb', 'w') as f:
    app.PDBFile.writeFile(simulation.topology, min_positions, f)

# レポーター(ログ出力)の設定
simulation.reporters.append(app.DCDReporter('trajectory.dcd', 100))
simulation.reporters.append(app.StateDataReporter(
    sys.stdout, 100, step=True, potentialEnergy=True, temperature=True
))

# 短いMD実行 (1000ステップ)
print("MD Running...")
simulation.step(10000)

# 最終構造を保存
final_positions = simulation.context.getState(getPositions=True).getPositions()
with open('final_cpu.pdb', 'w') as f:
    app.PDBFile.writeFile(simulation.topology, final_positions, f)

print("Done.")
