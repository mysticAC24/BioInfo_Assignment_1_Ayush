import json
import os

def load_restriction_sites():
    """
    Load restriction enzyme recognition sequences from JSON file
    """
    base_dir = os.path.dirname(__file__)
    json_file = os.path.join(base_dir, "restriction_sites.json")

    with open(json_file, "r") as f:
        return json.load(f)
