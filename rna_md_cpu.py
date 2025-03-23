#!/usr/bin/env python

"""
# RNA MD simulation script with optional solvent models.
# 例1: 真空 (NoCutoff) モデルで 300K, 拘束付き 5000 ステップ, 拘束なし 20000 ステップを 2 fs で
python argparse.py \
  --pdb_file data/1zih.pdb \
  --solvent vacuum \
  --temperature 300 \
  --friction 1.0 \
  --timestep 2.0 \
  --restraint_k 5.0 \
  --restraint_steps 5000 \
  --production_steps 10000 \
  --output_frames 200
"""

import sys
import argparse

import openmm as mm
from openmm import app, unit


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="RNA MD simulation script with optional solvent models."
    )

    # ======== 入力ファイル系 ========
    parser.add_argument(
        "--pdb_file", type=str, default="data/1zih.pdb",
        help="Input PDB file path. (default: data/1zih.pdb)"
    )

    # ======== 溶媒モデルの選択 (vacuum, implicit, explicit) ========
    parser.add_argument(
        "--solvent", type=str,
        choices=["vacuum", "implicit", "explicit"],
        default="explicit",
        help="Solvent model choice: vacuum, implicit (GBn2), or explicit (TIP3P+PME). (default: explicit)"
    )

    # ======== シミュレーション設定 ========
    parser.add_argument(
        "--temperature", type=float, default=300.0,
        help="Simulation temperature in Kelvin. (default: 300 K)"
    )
    parser.add_argument(
        "--friction", type=float, default=1.0,
        help="Langevin friction coefficient (1/ps). (default: 1.0 /ps)"
    )
    parser.add_argument(
        "--timestep", type=float, default=2.0,
        help="Integration timestep in femtoseconds. (default: 2.0 fs)"
    )

    # ======== 段階的(ステージング) MDステップ数 ========
    parser.add_argument(
        "--restraint_k", type=float, default=5.0,
        help="Positional restraint force constant on heavy atoms (kcal/mol/A^2) for initial stage. (default: 5.0)"
    )
    parser.add_argument(
        "--restraint_steps", type=int, default=5000,
        help="Number of MD steps under positional restraint. (default: 5000)"
    )
    parser.add_argument(
        "--production_steps", type=int, default=10000,
        help="Number of MD steps in production (unrestrained) stage. (default: 10000)"
    )

    # ======== 出力設定 ========
    parser.add_argument(
        "--output_frames", type=int, default=200,
        help="Number of total frames to be saved in the trajectory. (default: 200)"
    )

    args = parser.parse_args()
    return args


def print_configuration(args):
    """設定内容の表示。"""
    print("Input PDB file:", args.pdb_file)
    print("Solvent model:", args.solvent)
    print("Temperature:", args.temperature, "K")
    print("Friction:", args.friction, "/ps")
    print("Timestep:", args.timestep, "fs")
    print("Positional restraint force constant:", args.restraint_k, "kcal/mol/A^2")
    print("Positional restraint steps:", args.restraint_steps)
    print("Production steps:", args.production_steps)
    print("Output frames:", args.output_frames)


def create_system_and_integration(pdb, args):
    """
    ForceField やモデラー(溶媒追加)を適切に設定して、
    OpenMM の System, Integrator, Topology, Positions を作成する。
    """

    # --- ForceField の読み込み ---
    # RNA 用 ForceField (amber14 + RNA.OL3)
    # 明示的溶媒の場合は TIP3P (amber14/tip3p.xml) も必要
    if args.solvent == "implicit":
        # Implicit (GBn2) 用の追加ファイル
        forcefield = app.ForceField(
            'amber14-all.xml',
            'amber14/RNA.OL3.xml',
            'amber14/implicit/gbn2.xml'
        )
    else:
        # vacuum または explicit
        forcefield = app.ForceField(
            'amber14-all.xml',
            'amber14/RNA.OL3.xml'
        )
        if args.solvent == "explicit":
            forcefield.loadFile('amber14/tip3p.xml')

    # --- 溶媒設定 ---
    if args.solvent == "explicit":
        # PBC + PME 用に PDBFile -> Modeller に変換して水追加
        modeller = app.Modeller(pdb.topology, pdb.positions)
        # ボックスに 1.0 nm のマージンを取り、0.15M のイオンを加える例
        modeller.addSolvent(
            forcefield,
            model='tip3p',  # 大文字に変更
            padding=1.0 * unit.nanometer,
            ionicStrength=0.15 * unit.molar
        )
        topology = modeller.topology
        positions = modeller.positions

        # システム構築 (PME)
        system = forcefield.createSystem(
            topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0 * unit.nanometer,
            constraints=app.HBonds,
            removeCMMotion=True
        )

    elif args.solvent == "implicit":
        # 真空 + GB
        topology = pdb.topology
        positions = pdb.positions
        # constraints=HBonds で 2 fs が無難
        system = forcefield.createSystem(
            topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds,
            removeCMMotion=True
        )
    else:
        # vacuum
        topology = pdb.topology
        positions = pdb.positions
        # 真空状態でも、ここでは constraints=HBonds にして 2 fs を許容
        system = forcefield.createSystem(
            topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds,
            removeCMMotion=True
        )

    # --- Integrator 設定 (Langevin) ---
    temperature = args.temperature * unit.kelvin
    friction = args.friction / unit.picosecond  # 1/ps
    timestep = args.timestep * unit.femtoseconds

    integrator = mm.LangevinIntegrator(temperature, friction, timestep)

    return system, integrator, topology, positions


def add_positional_restraints(system, topology, positions, force_k_kcal):
    """
    指定した force_k (kcal/mol/A^2) で重原子に対する位置拘束 (positional restraints) を追加。
    OpenMM のエネルギー単位 (kJ/mol/nm^2) に変換して CustomExternalForce を適用します。
    """
    # kcal/mol/A^2 => kJ/mol/nm^2 に変換
    # 1 kcal/mol = 4.184 kJ/mol, 1 Å = 0.1 nm なので
    # force_k (kcal/mol/A^2) * 4.184 / (0.1^2) = force_k * 418.4
    force_k = force_k_kcal * 418.4  # [kJ/mol/nm^2]

    restraint = mm.CustomExternalForce(
        "0.5 * k * ((x - x0)^2 + (y - y0)^2 + (z - z0)^2)"
    )
    restraint.addPerParticleParameter("k")
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    for atom in topology.atoms():
        # 重原子のみ拘束 (RNA なので H 以外)
        if atom.element.symbol != 'H':
            pos = positions[atom.index]  # positions リストから座標を取得
            restraint.addParticle(
                atom.index,
                [force_k, pos.x, pos.y, pos.z]
            )

    system.addForce(restraint)


def run_minimization(simulation, max_iterations=500):
    """簡易エネルギー最小化。"""
    print("Energy Minimization...")
    simulation.minimizeEnergy(maxIterations=max_iterations)


def run_md(simulation, steps, report_interval, dcd_name, state_stdout):
    """MD を run する。レポーターを設定してステップを実行。"""
    print("MD Running...")
    simulation.reporters.clear()
    simulation.reporters.append(app.DCDReporter(dcd_name, report_interval))
    simulation.reporters.append(
        app.StateDataReporter(
            state_stdout, report_interval,
            step=True, potentialEnergy=True, temperature=True,
            progress=True, remainingTime=True, speed=True, totalSteps=steps
        )
    )
    simulation.step(steps)


def main():
    args = parse_arguments()
    print_configuration(args)

    # --- 入力PDBファイル読み込み ---
    pdb = app.PDBFile(args.pdb_file)
    pdb_id = args.pdb_file.split('/')[-1].split('.')[0]

    # --- ForceField & System 構築 ---
    system, integrator, topology, positions = create_system_and_integration(pdb, args)

    # --- シミュレーションオブジェクト準備 ---
    platform = mm.Platform.getPlatformByName('CPU')  # 環境に合わせて変更
    simulation = app.Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(positions)

    # --- 最小化 ---
    run_minimization(simulation, max_iterations=1000)

    # --- 最小化構造を PDB に保存 ---
    min_state = simulation.context.getState(getPositions=True)
    with open(f"result/{pdb_id}_minimized.pdb", "w") as f:
        app.PDBFile.writeFile(topology, min_state.getPositions(), f)

    # --- 拘束付きプレMD (stage 1) ---
    if args.restraint_steps > 0 and args.restraint_k > 0.0:
        print(f"Adding positional restraints (k={args.restraint_k} kcal/mol/A^2)...")
        add_positional_restraints(system, topology, positions, args.restraint_k)
        simulation.context.reinitialize(preserveState=True)

        # 出力フレーム数を固定: interval = max(1, steps // output_frames)
        interval_stage1 = max(1, args.restraint_steps // args.output_frames)

        print("Running restrained MD...")
        run_md(simulation, args.restraint_steps, interval_stage1,
               dcd_name=f"result/{pdb_id}_restrained.dcd",
               state_stdout=sys.stdout)

        # 拘束付き段階終了後、拘束力を除去
        system_forces = system.getForces()
        for f in system_forces:
            if isinstance(f, mm.CustomExternalForce):
                system.removeForce(system_forces.index(f))
                break

        simulation.context.reinitialize(preserveState=True)

    # --- 拘束なし生産MD (stage 2) ---
    if args.production_steps > 0:
        interval_prod = max(1, args.production_steps // args.output_frames)
        print("Running production MD (unrestrained)...")
        run_md(simulation, args.production_steps, interval_prod,
               dcd_name=f"result/{pdb_id}_production.dcd",
               state_stdout=sys.stdout)

    # --- 終了時の座標保存 ---
    final_state = simulation.context.getState(getPositions=True)
    with open(f"result/{pdb_id}_final.pdb", "w") as f:
        app.PDBFile.writeFile(topology, final_state.getPositions(), f)

    print("All Done.")


if __name__ == "__main__":
    main()