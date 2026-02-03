import os, glob, time
import sys
import pandas as pd
from multiprocessing import Pool
import pyrosetta
from pyrosetta import *
from pyrosetta.teaching import *

#Protocol Includes
#from rosetta.protocols import relax as rel
#from rosetta.protocols.antibody.residue_selector import CDRResidueSelector
#from rosetta.protocols.antibody import *
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

# PyRosetta 初始化
init('-use_input_sc -input_ab_scheme AHo_Scheme -ignore_unrecognized_res \
     -ignore_zero_occupancy false -load_PDB_components false  -relax:default_repeats 2 -no_fconfig')

def bool_type(bool_str: str):
    bool_str_lower = bool_str.lower()
    if bool_str_lower in ('false', 'f', 'no', 'n', '0'):
        return False
    elif bool_str_lower in ('true', 't', 'yes', 'y', '1'):
        return True
    else:
        raise ValueError(f'Cannot interpret {bool_str} as bool')

# file_name 函数在 main 中未使用，保持不变
def file_name(file_dir):
    L = []
    for root, dirs, files in os.walk(file_dir):
        for file in files:
            if os.path.splitext(file)[1] == '.pdb':
                L.append(os.path.join(root , file))
        return L

def create_score_functions():
    """Create all score functions needed for relaxation and scoring.

    Must be called per-worker since PyRosetta C++ objects are not picklable.
    """
    sfxns = {}
    sfxns['scorefxn'] = get_score_function()
    sfxns['fa_atr'] = ScoreFunction(); sfxns['fa_atr'].set_weight(fa_atr, 1.0)
    sfxns['fa_rep'] = ScoreFunction(); sfxns['fa_rep'].set_weight(fa_rep, 1.0)
    sfxns['fa_intra_rep'] = ScoreFunction(); sfxns['fa_intra_rep'].set_weight(fa_intra_rep, 1.0)
    sfxns['fa_sol'] = ScoreFunction(); sfxns['fa_sol'].set_weight(fa_sol, 1.0)
    sfxns['lk_ball_wtd'] = ScoreFunction(); sfxns['lk_ball_wtd'].set_weight(lk_ball_wtd, 1.0)
    sfxns['fa_intra_sol'] = ScoreFunction(); sfxns['fa_intra_sol'].set_weight(fa_intra_sol, 1.0)
    sfxns['fa_elec'] = ScoreFunction(); sfxns['fa_elec'].set_weight(fa_elec, 1.0)
    sfxns['hbond_lr_bb'] = ScoreFunction(); sfxns['hbond_lr_bb'].set_weight(hbond_lr_bb, 1.0)
    sfxns['hbond_sr_bb'] = ScoreFunction(); sfxns['hbond_sr_bb'].set_weight(hbond_sr_bb, 1.0)
    sfxns['hbond_bb_sc'] = ScoreFunction(); sfxns['hbond_bb_sc'].set_weight(hbond_bb_sc, 1.0)
    sfxns['hbond_sc'] = ScoreFunction(); sfxns['hbond_sc'].set_weight(hbond_sc, 1.0)
    sfxns['dslf_fa13'] = ScoreFunction(); sfxns['dslf_fa13'].set_weight(dslf_fa13, 1.0)
    sfxns['rama_prepro'] = ScoreFunction(); sfxns['rama_prepro'].set_weight(rama_prepro, 1.0)
    sfxns['p_aa_pp'] = ScoreFunction(); sfxns['p_aa_pp'].set_weight(p_aa_pp, 1.0)
    sfxns['fa_dun'] = ScoreFunction(); sfxns['fa_dun'].set_weight(fa_dun, 1.0)
    sfxns['omega'] = ScoreFunction(); sfxns['omega'].set_weight(omega, 1.0)
    sfxns['pro_close'] = ScoreFunction(); sfxns['pro_close'].set_weight(pro_close, 1.0)
    sfxns['yhh_planarity'] = ScoreFunction(); sfxns['yhh_planarity'].set_weight(yhh_planarity, 1.0)
    sfxns['ref'] = ScoreFunction(); sfxns['ref'].set_weight(ref, 1.0)
    sfxns['rg'] = ScoreFunction(); sfxns['rg'].set_weight(rg, 1)
    return sfxns


def process_single_pdb(task):
    """Worker function to process a single PDB file.

    Args:
        task: dict with keys 'pdb_path', 'args', 'interface' (str or None)

    Returns:
        (pdb_name, temp_dict) on success, None on failure.
    """
    pdb_path = task['pdb_path']
    args = task['args']
    interface = task['interface']

    try:
        sfxns = create_score_functions()
        scorefxn = sfxns['scorefxn']
        fixed_chain = args.fixed_chain.split("_")

        start = time.time()
        pdb_name = os.path.basename(pdb_path).replace(".pdb", "")
        pose = pose_from_pdb(pdb_path)
        original_pose = pose.clone()

        # FastRelax
        if args.relax:
            fr = FastRelax()
            fr.set_scorefxn(scorefxn)
            fr.max_iter(int(args.max_iter))
            movemap = MoveMap()
            movemap.set_bb(True)
            movemap.set_chi(True)

            if args.fixbb:
                for i in range(1, pose.total_residue() + 1):
                    chain = pose.pdb_info().chain(i)
                    if chain in fixed_chain:
                        movemap.set_bb(i, False)
                    else:
                        movemap.set_bb(i, True)
            fr.set_movemap(movemap)

        if (not os.getenv("DEBUG")) and (args.relax):
            if 'fr' in locals():
                fr.apply(pose)

        # Interface analysis
        ia = InterfaceAnalyzerMover()
        if interface is not None:
            ia.set_interface(interface)
        else:
            print(f"Warning: No interface information provided for {pdb_name}. Using default interface.")
        ia.set_skip_reporting(True)
        ia.apply(pose)

        time_consumed = time.time() - start

        # Collect scores - preserve exact key names including 'hbond_sr_bb(pose)' for CSV compatibility
        temp_dict = {
            'relaxed': scorefxn(pose),
            'interface_score': ia.get_interface_dG(),
            'original': scorefxn(original_pose),
            'delta': scorefxn(pose) - scorefxn(original_pose),
            'fa_atr': sfxns['fa_atr'](pose),
            'fa_rep': sfxns['fa_rep'](pose),
            'fa_intra_rep': sfxns['fa_intra_rep'](pose),
            'fa_sol': sfxns['fa_sol'](pose),
            'lk_ball_wtd': sfxns['lk_ball_wtd'](pose),
            'fa_intra_sol': sfxns['fa_intra_sol'](pose),
            'fa_elec': sfxns['fa_elec'](pose),
            'hbond_lr_bb': sfxns['hbond_lr_bb'](pose),
            'hbond_sr_bb(pose)': sfxns['hbond_sr_bb'](pose),
            'hbond_bb_sc': sfxns['hbond_bb_sc'](pose),
            'hbond_sc': sfxns['hbond_sc'](pose),
            'dslf_fa13': sfxns['dslf_fa13'](pose),
            'rama_prepro': sfxns['rama_prepro'](pose),
            'p_aa_pp': sfxns['p_aa_pp'](pose),
            'fa_dun': sfxns['fa_dun'](pose),
            'omega': sfxns['omega'](pose),
            'pro_close': sfxns['pro_close'](pose),
            'yhh_planarity': sfxns['yhh_planarity'](pose),
            'ref': sfxns['ref'](pose),
            'get_complexed_sasa': ia.get_complexed_sasa(),
            'get_interface_delta_sasa': ia.get_interface_delta_sasa(),
            'time_consumed': time_consumed,
        }

        if args.dump_pdb:
            pose.dump_pdb(os.path.join(args.output_dir, 'relax_' + os.path.basename(pdb_path)))

        return (pdb_name, temp_dict)

    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")
        return None


def main(args):
    # Build task list
    tasks = []

    if args.csv_path != "":
        # Mode 1: CSV file with PDB paths and optional interface info
        df = pd.read_csv(args.csv_path)

        if len(df) == 0:
            print(f"CSV file at {args.csv_path} is empty.")
            return

        print(f"Processing {len(df)} pdb from CSV")
        for idx, row in df.iterrows():
            pdb_path = row["pdb"]
            interface = None
            if ("ligand" in df.columns) and ("receptor" in df.columns):
                interface = row["ligand"].replace(",", "") + "_" + row["receptor"].replace(",", "")
            else:
                print(f"Warning: 'ligand' or 'receptor' column missing for {os.path.basename(pdb_path)}. Using default interface.")
            tasks.append({'pdb_path': pdb_path, 'args': args, 'interface': interface})
    else:
        # Mode 2: PDB directory
        all_pdb_name = sorted(glob.glob(os.path.join(args.pdb_dir, '*.pdb')))

        if len(all_pdb_name) > 0:
            print(f"Processing {len(all_pdb_name)} pdb in {args.pdb_dir}")
        else:
            print(f"No pdb available in {args.pdb_dir}")
            return

        for pdb_path in all_pdb_name:
            tasks.append({'pdb_path': pdb_path, 'args': args, 'interface': None})

    # Process tasks
    if args.num_workers > 1:
        print(f"Using {args.num_workers} parallel workers")
        with Pool(args.num_workers) as pool:
            results = pool.map(process_single_pdb, tasks)
    else:
        results = [process_single_pdb(t) for t in tasks]

    # Collect results, filter failures
    results = [r for r in results if r is not None]

    if len(results) > 0:
        output = [[pdb_name] + [v for v in temp_dict.values()] for pdb_name, temp_dict in results]
        columns = ['pdb_name'] + list(results[0][1].keys())
        score_df = pd.DataFrame(output, columns=columns)
        score_df.to_csv(os.path.join(args.output_dir, f'rosetta_complex_{args.batch_idx}.csv'), index=False)
    else:
        print("No results to write. Output list is empty.")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_dir',type=str, default="", help="input pdb path.")
    parser.add_argument('--csv_path',type=str, default="", help="save absolute path of pdb file in the 'pdb' column")
    parser.add_argument('--output_dir',type=str, help="output relaxed pdb path", required=True)
    parser.add_argument('--dump_pdb', default= False, type=bool_type, help="dump pdb or not")
    parser.add_argument('--batch_idx', type=int, default= False, help="batch index for output filename")
    parser.add_argument('--relax', default= False, type=bool_type, help="run relax or not")
    parser.add_argument('--fixbb', default= False, type=bool_type, help="fix backbone or not")
    parser.add_argument('--fixed_chain', default="", type=str, help="chains whose backbone should be fixed (e.g., 'A_B')")
    parser.add_argument('--max_iter', default=1, type=int, help="max iteration for FastRelax")
    parser.add_argument('--num_workers', default=len(os.sched_getaffinity(0)), type=int, help="number of parallel workers for processing PDBs (default: all available CPUs)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    main(args)
