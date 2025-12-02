# beeline_config_utils.py
from ruamel.yaml import YAML
from ruamel.yaml.scalarstring import DoubleQuotedScalarString
from pathlib import Path

yaml = YAML()
Quoted = DoubleQuotedScalarString

def quoted_constructor(loader, node):
    return Quoted(loader.construct_scalar(node))

yaml.constructor.add_constructor('tag:yaml.org,2002:str', quoted_constructor)
yaml.default_flow_style = False


# -----------------------------
# Master dictionary of all algorithms with defaults
# -----------------------------
ALL_ALGOS = {
    "GENIE3": {"should_run": [True]},
    "GRNBOOST2": {"should_run": [True]},
    "PIDC": {"should_run": [True]},
    "GRNVBEM": {"should_run": [False]},
    "PPCOR": {"should_run": [False], "pVal": [0.05]},
    "SCODE": {"should_run": [False], "z": [10], "nIter": [1000], "nRep": [6]},
    "SCNS": {"should_run": [False]},
    "SINCERITIES": {"should_run": [False], "nBins": [10]},
    "LEAP": {"should_run": [False], "maxLag": [0.33]},
    "GRISLI": {"should_run": [False], "L": [10], "R": [3000], "alphaMin": [0.0]},
    "SINGE": {
        "should_run": [False],
        "lambda": [0.01],
        "dT": [15],
        "num_lags": [5],
        "kernel_width": [0.5],
        "prob_zero_removal": [0],
        "prob_remove_samples": [0.0],
        "family": ["gaussian"],
        "num_replicates": [6],
    },
    "SCRIBE": {
        "should_run": [False],
        "delay": ["5"],
        "method": ["ucRDI"],
        "lowerDetectionLimit": [0],
        "expressionFamily": ["uninormal"],
        "log": [False],
        "ignorePT": [True],
    },
    "SCSGL": {
        "should_run": [False],
        "pos_density": [0.45],
        "neg_density": [0.45],
        "assoc": ["correlation"],
    },
    "CELLORACLE": {"should_run": [False]},
    "CELLORACLE_ATAC": {"should_run": [False]},
    "CELLORACLE_alpha01": {"should_run": [False]},
    "CELLORACLE_alpha1": {"should_run": [False]},
}


# -----------------------------
# Helper functions
# -----------------------------

def load_config(path):
    """Load a BEELINE YAML config file."""
    with open(path) as f:
        return yaml.load(f)


def save_config(config, path):
    """Save a BEELINE config dictionary to a YAML file."""
    with open(path, "w") as f:
        yaml.dump(config, f)
    print(f"Config saved to {Path(path).resolve()}")


def add_datasets(config, new_datasets):
    """Add new dataset names with consistent quoted formatting."""
    Quoted = DoubleQuotedScalarString

    for ds in new_datasets:
        config["input_settings"]["datasets"].append({
            "name": Quoted(ds),
            "exprData": Quoted("ExpressionData.csv"),
            "cellData": Quoted("PseudoTime.csv"),
            "trueEdges": Quoted("refNetwork.csv"),
        })

def set_algorithms(config, algos_to_run, all_algos_dict):
    """Replace algorithms list, keeping all formatting consistent."""
    Quoted = DoubleQuotedScalarString

    new_algo_list = []
    for algo in algos_to_run:
        if algo in all_algos_dict:
            params = all_algos_dict[algo].copy()
            # Ensure should_run is explicitly True and has the correct style
            params["should_run"] = yaml.seq([True])
            params["should_run"].fa.set_flow_style()

            # Apply flow style to all other parameter lists
            for key, value in params.items():
                if isinstance(value, list):
                    seq = yaml.seq(value)
                    seq.fa.set_flow_style()
                    params[key] = seq

            new_algo_list.append({
                "name": Quoted(algo),
                "params": params
            })
    config["input_settings"]["algorithms"] = new_algo_list


def toggle_algorithm(config, algo_name, run=True):
    """Toggle should_run, keeping inline list formatting."""
    Quoted = DoubleQuotedScalarString

    algos = config["input_settings"]["algorithms"]
    for algo in algos:
        if algo["name"] == algo_name:
            # Create a new list and apply flow style to it
            new_val = yaml.seq([run])
            new_val.fa.set_flow_style()
            algo["params"]["should_run"] = new_val
            break


def init_config(datasets, algos_to_run, all_algos_dict=ALL_ALGOS,
                dataset_dir="default_dataset_dir",
                input_dir="inputs",
                output_dir="outputs",
                output_prefix="default_output_prefix"):
    """
    Create a BEELINE config dictionary with quotes on strings and inline lists.
    """
    
    # Use an alias for DoubleQuotedScalarString to keep the code clean.
    Quoted = DoubleQuotedScalarString

    # --- Build the Algorithms Section ---
    # This loop handles the special formatting for algorithm parameters.
    algo_entries = []
    for algo_name in algos_to_run:
        if algo_name in all_algos_dict:
            params = all_algos_dict[algo_name].copy()
            
            # This part puts parameter lists on a single line (e.g., [True]).
            for key, value in params.items():
                if isinstance(value, list):
                    seq = yaml.seq(value)
                    seq.fa.set_flow_style()
                    params[key] = seq
            
            algo_entries.append({"name": Quoted(algo_name), "params": params})

    # --- Build the Datasets Section ---
    dataset_entries = [
        {
            "name": Quoted(ds),
            "exprData": Quoted("ExpressionData.csv"),
            "cellData": Quoted("PseudoTime.csv"),
            "trueEdges": Quoted("refNetwork.csv"),
        }
        for ds in datasets
    ]
    
    # --- Assemble the Final Config Dictionary ---
    config = {
        "input_settings": {
            "input_dir": Quoted(input_dir),
            "dataset_dir": Quoted(dataset_dir),
            "datasets": dataset_entries,
            "algorithms": algo_entries,
        },
        "output_settings": {
            "output_dir": Quoted(output_dir),
            "output_prefix": Quoted(output_prefix),
        },
    }
    return config