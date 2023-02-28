from kstar import config
import os




config.NETWORK_DIR = config.update_network_directory("C:/Users/poss982/Documents/GitHub/Exp18_cell_lines/KSTAR/NETWORKS/NetworKIN")
print(config.NETWORK_DIR)

config.create_network_pickles()

config.install_resource_files()

config.check_configuration()