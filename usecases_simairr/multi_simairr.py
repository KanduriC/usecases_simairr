import random
import os
from itertools import product
from time import sleep
import argparse
import logging
from usecases_simairr.util import makedir_if_not_exists, write_yaml_file, parse_user_yaml

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--config_file', help='path to YAML specification file describing the desired multi '
                                                'simulation parameters', required=True)
args = parser.parse_args()


def gen_output_path(super_path, usecase, run, phenotype_burden):
    return os.path.join(super_path, usecase, run, "phenotype_burden_" + str(phenotype_burden), "data")


def run_simarirr(yaml_file):
    command = 'sim_airr -i ' + yaml_file
    exit_code = os.system(command)
    if exit_code != 0:
        raise RuntimeError(f"Running simAIRR failed:{command}.")


def multi_simairr(super_path, params_tuple, base_params_dict):
    makedir_if_not_exists(os.path.join(super_path, "simulation_scripts"), fail_if_exists=False)
    signal_sequence_files, replication_space, burden_and_pool_space = params_tuple
    for run, burden_and_pool_size in product(replication_space, burden_and_pool_space):
        burden, pool_size = burden_and_pool_size
        outpath = gen_output_path(super_path, signal_sequence_files["usecase_name"], run, burden)
        new_params = {"output_path": outpath, "seed": random.randint(9, 999),
                      "signal_sequences_file": signal_sequence_files["file_path"],
                      "phenotype_burden": burden, "phenotype_pool_size": pool_size}
        base_params_dict.update(new_params)
        yaml_file_path = os.path.join(super_path, "simulation_scripts",
                                      str(signal_sequence_files["usecase_name"]) +
                                      "_" + "phenotype_burden_" + str(burden) + "_" + str(run) + ".yaml")
        write_yaml_file(base_params_dict, yaml_file_path)
        sleep(10)
        logging.info(f'simulating the following config file {yaml_file_path}')
        run_simarirr(yaml_file_path)


def execute():
    config = parse_user_yaml(args.config_file)
    logging.basicConfig(format='%(asctime)s: %(levelname)s: %(message)s', level=logging.DEBUG,
                        filename=os.path.join(config['super_path'], "multi_simairr_log.txt"), filemode='a')
    burden_and_pools = config['burden_and_pools']
    reps = ["run1", "run2", "run3"]
    signal_file_params = {"usecase_name": config['usecase_name'],
                          "file_path": config['signal_file']}
    super_path = config['super_path']
    rem_list = ['burden_and_pools', 'usecase_name', 'signal_file', 'super_path']
    for key in rem_list:
        del config[key]
    multi_simairr(super_path=super_path,
                  params_tuple=(signal_file_params, reps, burden_and_pools), base_params_dict=config)
