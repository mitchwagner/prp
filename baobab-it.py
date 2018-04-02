import yaml
import copy
import argparse
import itertools

import src.external.utils.baobab_utils as baobab

def main():
    opts = parse_arguments()

    config_map = None

    # 1) Read the yaml file
    with open(opts.config, 'r') as handle:
        config_map = yaml.load(handle)

    # 2) Create many yaml files
    algorithms = config_map["input_settings"]["algorithms"]
    #print(algorithms)

    del config_map["input_settings"]["algorithms"]

    #for algorithm in algorithms:
    #    print(algorithm)

    for algorithm in algorithms:
        name = algorithm["name"]
        combos = [dict(zip(algorithm["params"], val))
            for val in itertools.product(
                *(algorithm["params"][param]
                    for param in algorithm["params"]))]

        for combo in combos:
            for item in combo:
                combo[item] = [combo[item]]

        # 3) Run reconstructions for each algorithm on baobab
        for i, combo in enumerate(combos):
            wrapped_combo = {"name" : algorithm["name"], "params": combo}
            #print(name, i, wrapped_combo)

            config_file = "configs/config" + name + "-" + str(i) + ".yaml"

            print("Running", config_file)

            qsub_file = "/home/mitchw94/Desktop/pathway-reconstruction-pipeline/qsub/qsub" + name + "-" + str(i)

            cop = copy.deepcopy(config_map)
            cop["input_settings"]["algorithms"] = [wrapped_combo]
            #print(cop)

            with open(config_file, 'w') as f:
                yaml.dump(cop, f, default_flow_style=True)

            baobab.writeQsubFile(
                ["cd ~/Desktop/pathway-reconstruction-pipeline"
                ,"source venv/bin/activate"
                ,"venv/bin/python pipeline.py " 
                "--config=" + config_file + " "
                ],
                qsub_file,
                name="Run:" + config_file,
                ppn=8
                )
            baobab.submitQsubFile(qsub_file)


def parse_arguments():
    parser = get_parser()
    opts = parser.parse_args()

    return opts 


def get_parser():

    parser = argparse.ArgumentParser(                                                                                                                                        
    description='Run pathway reconstruction pipeline...on baobab!!!')                     

    parser.add_argument('--config', default='config.yaml',
        help='Configuration file')     

    return parser


if __name__ == '__main__':
    main()
