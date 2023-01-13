import argparse
import glob
import os.path
import copy
from time import sleep
from usecases_simairr.util import parse_user_yaml, put_values_in_target_key, makedir_if_not_exists, write_yaml_file

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--sim_data_path', help='super path containing all simulated_repertoire directories', required=True)
parser.add_argument('-m', '--ml_yaml_file', help='ML config file for immune-ml', required=True)
args = parser.parse_args()


def run_immuneml(yaml_file, output_folder):
    command = f'immune-ml {yaml_file} {output_folder}'
    exit_code = os.system(command)
    if exit_code != 0:
        raise RuntimeError(f"Running immune-ml failed:{command}.")


def execute():
    dir_list = glob.glob(f"{args.sim_data_path}/**/simulated_repertoires", recursive=True)
    metadata_list = [os.path.join(simdata_path, "metadata.csv") for simdata_path in dir_list]
    ml_config = parse_user_yaml(args.ml_yaml_file)
    ml_configs_list = [put_values_in_target_key(copy.deepcopy(ml_config), "path", simdata_path) for simdata_path in
                       dir_list]
    mod_ml_configs = [put_values_in_target_key(ml_config_dict, "metadata_file", meta_path)
                      for ml_config_dict, meta_path in zip(ml_configs_list, metadata_list)]
    ml_out_paths = [os.path.join(simdata_path.replace("/data/simulated_repertoires", ""), "ml_output") for simdata_path
                    in dir_list]
    ml_config_out_fns = [os.path.join(ml_out_path, "ml_config_file.yaml") for ml_out_path in ml_out_paths]
    for outpath in ml_out_paths:
        makedir_if_not_exists(outpath)
    for mod_ml_config, file_path in zip(mod_ml_configs, ml_config_out_fns):
        write_yaml_file(mod_ml_config, file_path)
        sleep(10)
        run_immuneml(file_path, os.path.join(os.path.dirname(file_path), "immuneml_output"))

