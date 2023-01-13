import os
import sys
import yaml


def parse_user_yaml(yaml_file):
    with open(yaml_file, "r") as yaml_file:
        try:
            yaml_obj = yaml.load(yaml_file, Loader=yaml.FullLoader)
            assert yaml_obj is not None, "The supplied yaml file," + yaml_file + ", is empty"
        except Exception as e:
            print(
                "Error: that looks like an invalid yaml file. Consider validating your yaml file using one of the "
                "online yaml validators; for instance: https://jsonformatter.org/yaml-validator")
            print("Exception: %s" % str(e))
            sys.exit(1)
    return yaml_obj


def put_values_in_target_key(yaml_dict, target_key, insert_value):
    for key, value in yaml_dict.items():
        if isinstance(value, dict):
            put_values_in_target_key(value, target_key, insert_value)
        elif key == target_key and isinstance(value, str):
            yaml_dict[key] = insert_value
    return yaml_dict


def write_yaml_file(yaml_dict, out_file_path):
    with open(out_file_path, "w+") as yaml_file:
        yaml.dump(yaml_dict, yaml_file)


def makedir_if_not_exists(some_path, fail_if_exists=False):
    if os.path.exists(some_path):
        files = [fn for fn in os.listdir(some_path) if
                 os.path.isfile(os.path.join(some_path, fn))]
        files_exists = f"Output folder may already contain relevant files: {some_path}"
        if fail_if_exists:
            assert len(files) == 0, files_exists
    else:
        os.makedirs(some_path)


