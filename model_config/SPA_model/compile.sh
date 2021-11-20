#!/bin/sh

ifort -fpp -vec-report=0 -i-dynamic -O3 global_variables.f90 log_tools.f90 math_tools.f90 crops.f90 spa_io_csv.f90 config_tools.f90 spa_cmd_line.f90 spa_config.f90 soil_air.f90 spa_initialise.f90 linked_lists.f90 spa_io.f90 light.f90 allocate_carbon.f90 soil_functions.f90 leaf.f90 canopy.f90 main.f90 -o spa

mv spa ../
