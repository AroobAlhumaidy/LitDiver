import os
import yaml

# Load YAML config file or create a default one
def load_config(config_file="config.yml"):
    default_config = {
        "email": "your-email@example.com",
        "max_results": 10000,
        "date_range": None,
        "download_pdfs": True,
        "output_dir": "output_litdiver"
    }

    if not os.path.exists(config_file):
        with open(config_file, "w") as f:
            yaml.dump(default_config, f)
        print(f"Default config.yml created. Please update your email.")
        return default_config

    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    # Merge missing keys from default config
    for key, value in default_config.items():
        if key not in config:
            config[key] = value

    return config
