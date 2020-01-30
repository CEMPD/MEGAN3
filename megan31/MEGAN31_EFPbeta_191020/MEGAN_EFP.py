#!/usr/bin/env python
""" MEGAN EF Processor

Description:
    The MEGAN Emission Factor Processor (MEGAN EFP)
    calculate emission factors (EF) and light dependence factors (LDF),
    for specific vegetation types and calculated the grid average EF and LDF.
    This is accomplished by synthesizing available measurements and is weighted by
    landcover. This program is built on
    Python (version 2.7) and SQLite.

Required Python Packages:
    pandas
    numpy

Local Modules (located in ./src/):
    run_M3EFP
    db_functions
    M3GEFP
    M3VTEF
    M3LDF

Usage:
    $ python MEGAN_EFP.py

Development:
    This program was developed by Alex Guenther (alex.guenther@uci.edu) in ACCESS VBA  
    and then converted to python + SQLite by Andy Wentland (Ramboll) in 2017
    Revised by Alex Guenther in August 2019 (This is MEGAN_EFP Version 3.1_beta)    


"""
# Import libraries and modules
import sys
sys.path.append("./src")
import run_M3EFP as efp



# ~~~~~~~~~~~~~ BEGIN USER OPTIONS ~~~~~~~~~~~~~ #

# ~~~~ SCENARIO NAME ~~~~ #

scen_name = "TCEQ12.2019b.J4"

# ~~~~ DATABASE ~~~~ #
# Specify path and name of M3VTEF, M3GEFP, M3LDF databases to be created
M3VTEF_database = "./database/M3VTEF_database."+scen_name+".db"
M3GEFP_database = "./database/M3GEFP_database."+scen_name+".db"
M3LDF_database = "./database/M3LDF_database."+scen_name+".db"
M3GLDF_database = "./database/M3GLDF_database."+scen_name+".db"


# ~~~~ INPUT TABLES ~~~~ #
# Directories with premade CSV files
# to be converted into database tables
M3VTEF_inputs_path = "./inputs/M3VTEF/"
M3GEFP_inputs_path = "./inputs/M3GEFP/"
M3LDF_inputs_path = "./inputs/M3LDF/"
M3GLDF_inputs_path = "./inputs/M3GLDF/"

VegDB_tables = "./inputs/VegDB_Tables/"
GridDB_tables = "./inputs/GridDB_Tables/"

# input file names
# under VegDB_tables
Description_Vegetation = "Description_Vegetation.2019a.csv"
#Description_Vegetation = "Description_Vegetation.TEOW19.csv"
#Description_Vegetation = "Description_Vegetation.2x2.csv"
#DB_SLA = "DB_SLA.TEOW19.csv"
#DB_LDF = "DB_LDF.TEOW19.csv"
DB_SLA = "DB_SLA.2019a.csv"
DB_LDF = "DB_LDF.2019a.csv"
DB_emissions = "DB_emissions.2019b.csv"
#DB_emissions = "DB_emissions.TEOW19.csv"
#DB_emissions = "DB_emissions.2x2.csv"

# Under GridDB_tables
#Ecotype_Crop = "Speciation.Crop.TEOW19.csv"
#Ecotype_Herb = "Speciation.Herb.TEOW19.csv"
#Ecotype_Shrub = "Speciation.Shrub.TEOW19.csv"
#Ecotype_Tree = "Speciation.Tree.TEOW19.csv"
#grid_ecotype = "grid.ecotype.Global_1deg.TEOW19b.csv"
#grid_growth_form = "grid.growth_form.Global_1deg.G2012.csv"
Ecotype_Crop = "Speciation.Crop.2019a.csv"
Ecotype_Herb = "Speciation.Herb.2019a.csv"
Ecotype_Shrub = "Speciation.Shrub.2019a.csv"
Ecotype_Tree = "Speciation.Tree.2019a.csv"
#grid_ecotype = "grid.ecotype.tceq_36km.TCEQ17.csv"
#grid_growth_form = "grid.growth_form.tceq_36km.Y2017.csv"
grid_ecotype = "grid.ecotype.tceq_12km.TCEQ17.csv"
grid_growth_form = "grid.growth_form.tceq_12km.Y2017.csv"

# ~~~~ SETTINGS ~~~~ #
# Total number of classes to loop through
# Default: 20
# If non-default number used, updates must be made in all submodules' concat_*_tables function
TotalClasses = 20

# Specify lowest J-rating to use in calculations
# 0 (lowest confidence) to 4 (highest confidence)
# Default: 0
#Jrating = 0
Jrating = 4


# ~~~~ OUTPUTS ~~~~ #
# Output directory
outputs_path = "./output/"

# ~~~~~~~~~~~~~ END USER OPTIONS ~~~~~~~~~~~~~ #

# Call m3efp_driver function in the run_M3EFP module to run program
if __name__ == "__main__":
    efp.m3efp_driver(scen_name, VegDB_tables, GridDB_tables, Description_Vegetation,
                     DB_LDF, DB_SLA, DB_emissions, Ecotype_Crop, Ecotype_Herb, Ecotype_Tree,
		     Ecotype_Shrub, grid_ecotype, grid_growth_form,
                     M3VTEF_database, M3GEFP_database, M3LDF_database, M3GLDF_database,
                     TotalClasses, Jrating, outputs_path)

