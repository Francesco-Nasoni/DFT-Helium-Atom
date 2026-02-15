from pathlib import Path
import numpy as np
import yaml

def load_config_yaml(cfg_file: str):
    config_path = Path.cwd() / cfg_file
    with open(config_path, 'r') as f:
        cfg = yaml.safe_load(f)["config"]
    
    if cfg["grid"]["r_max"]<0:
        raise ValueError(f"r_max must be greater or equal to 0")
    
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

    return cfg    
    

def append_csv(path, it, E_1, Etot, dE):
    if not Path(path).exists():
        Path(path).write_text("iter,eps_1s,E_tot,dE\n", encoding="utf-8")
    with open(path, "a", encoding="utf-8") as f:
        f.write(f"{it},{E_1:.12e},{Etot:.12e},{dE:.12e}\n")


def save_profiles(path, r, **fields):
    path = Path(path)

    names = ["r"]
    arrays = [r]

    # Unpack columns passed as keyword arg
    for name, arr in fields.items():
        if arr is None:
            continue
        names.append(name)
        arrays.append(arr)

    # Check lengths
    N = r.shape[0]
    for name, arr in zip(names, arrays):
        if arr.ndim != 1:
            raise ValueError(f"'{name}' is not 1D (ndim={arr.ndim})")
        if arr.shape[0] != N:
            raise ValueError(f"Different length: '{name}' has {arr.shape[0]} but r has {N}")

    values = np.column_stack(arrays)

    header = "  ".join(names)
    np.savetxt(
        path,
        values,
        header=header,
        fmt="%.12e",
    )