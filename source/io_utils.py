from pathlib import Path
import yaml

def load_config_yaml(cfg_file: str):
    config_path = Path.cwd() / cfg_file
    with open(config_path, 'r') as f:
        cfg = yaml.safe_load(f)["calculation_parameters"]
    
    if cfg["grid"]["r_min"] < 0 or cfg["grid"]["r_max"]<0:
        raise ValueError(f"r_min and r_max must be greater or equal to 0")

    if cfg["grid"]["r_min"] >= cfg["grid"]["r_max"]:
        raise ValueError(f"r_min must be less than r_max")
    
    if cfg["grid"]["h"] <= 0:
        raise ValueError(f"h must be positive")
    
    if cfg["bisect"]["E_min"] >= cfg["bisect"]["E_max"]:
        raise ValueError(f"E_min must be less than E_max")
    
    if cfg["bisect"]["rough_step"] <= 0:
        raise ValueError(f"rough_step must be positive")
    
    if not (0 < cfg["scf"]["mix_alpha"] <= 1):
        raise ValueError(f"mix_alpha must be in (0, 1]")
    
    if cfg["scf"]["max_iter"] <= 0:
        raise ValueError(f"max_iter must be positive")
    
    if cfg["scf"]["TOTEN_threshold"] <= 0:
        raise ValueError(f"TOTEN_convergence_threshold must be positive")
    
    if cfg["grid"]["r_min"] == 0:
        cfg["grid"]["r_min"] = 1e-12

    return cfg    