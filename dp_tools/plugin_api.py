import importlib
import os
import sys
import re
from pathlib import Path

def load_all_plugins(plugin_dir: Path):
    sys.path.append(str(plugin_dir))

    print ("DP TOOLS VALIDATION PLUGIN LOADER starting.  loading plugins:")
    plugins = {}
    for plugin_file in sorted(os.listdir(plugin_dir)):
        print(f"- {plugin_file}")
        if re.search('^dp_tools__', plugin_file):
            sys.path.append(plugin_file)
            module_name = plugin_file
            module = importlib.import_module(module_name)
            module_key = re.sub('^dp_tools__', '', module_name)
            plugins[module_key] = module
            print ("   %s" % ( module_name ))
    return plugins

def load_plugin(plugin_dir: Path):
    plugin_str = plugin_dir.name
    sys.path.append(str(plugin_dir))

    print ("DP TOOLS VALIDATION PLUGIN LOADER starting.  loading plugins:")
    print(f"- {plugin_str}")
    if re.search('^dp_tools__', plugin_str):
        module_name = plugin_str
        module = importlib.import_module(module_name)
        return module
    raise ValueError(f"Could not find validate plugin in {plugin_dir}")