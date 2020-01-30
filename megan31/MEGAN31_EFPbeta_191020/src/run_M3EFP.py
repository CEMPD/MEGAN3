#!/usr/bin/env python

# Import libraries, packages, and modules
from __future__ import division
import sqlite3
import os
import errno
import warnings
import numpy as np
import db_functions as dbf
import M3VTEF as mvtef
import M3GEFP as mgefp
import M3LDF as mldf
import M3GLDF as mgldf
warnings.simplefilter(action="ignore", category=FutureWarning)

# Metadata
__version__ = "1.5"

def make_dir(path):
    """

    Function to check if output path exists
    and create directories if needed
    :param path: directory path
    :return: new directory if one does not previously exists

    """
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def m3efp_driver(scen_name, VegDB_tables, GridDB_tables, Description_Vegetation,
                 DB_LDF, DB_SLA, DB_emissions, Ecotype_Crop, Ecotype_Herb, Ecotype_Tree,
                 Ecotype_Shrub, grid_ecotype, grid_growth_form,
                 M3VTEF_database_path, M3GEFP_database_path, M3LDF_database_path, M3GLDF_database_path,
                 TotalClasses, Jrating, outputs_path):
    """
    Main function to aggregate user options and call functions to create databases and output
    :param scen_name: Scenario Name
    :param VegDB_tables: Path to the CSVs to import into the M3VTEF & M3LDF database
    :param GridDB_tables: Path to the CSVs to import into the M3GEFP & M3GLDF database
    :param Description_Vegetation: vegetation description input file name
    :param DB_LDF: LDF input file name
    :param DB_SLA: SLA input file name
    :param DB_emissions: emissions input file name
    :Ecotype_Crop: input csv file name
    :Ecotype_Herb: input csv file name
    :Ecotype_Tree: input csv file name
    :Ecotype_Shrub: input csv file name
    :grid_ecotype: input csv file name
    :grid_growth_form: input csv file name
    :param M3VTEF_database_path: Path to SQLite database to be created for M3VTEF DB
    :param M3GEFP_database_path: Path to SQLite database to be created for M3GEFP DB
    :param M3LDF_database_path: Path to SQLite database to be created for M3LDF DB
    :param TotalClasses: Total number of classes to loop through when generating output
    :param Jrating: Jrating (user specified option)
    :param outputs_path: Path to output CSVs

    :return: CSVs by class and 1 CSV per database with concatenated tables containing all classes
    """

    print("\n MEGAN EF Processor")
    print(" Version %s \n" % __version__)

    # Check if CSV output directory exists and create if necessary
    make_dir(outputs_path)
    make_dir(outputs_path+"vegtype_EF_byClass")
    make_dir(outputs_path + "grid_EF_byClass")
    make_dir(outputs_path + "vegtype_LDF_byClass")

    # Print user specification
    print("\n\nUser Settings")
    print("Scenario Name: %s" % scen_name)
    print("M3VTEF Input Directory: %s" % VegDB_tables)
    print("M3GEFP Input Directory: %s" % GridDB_tables)
    print("M3LDF Input Directory: %s" % VegDB_tables)
    print("M3VTEF Database: %s" % M3VTEF_database_path)
    print("M3GEFP Database: %s" % M3GEFP_database_path)
    print("M3LDF Database: %s" % M3LDF_database_path)
    print("CSV Output Directory: %s" % outputs_path)
    print("J-rating: %s \n\n\n\n" % Jrating)

    # Define number of classes to loop through
    classids = np.arange(1, TotalClasses+1, 1)

    # Loop through each class
    print("Creating VegType EF and LDF databases")
    for cid in classids:
        vegtype_ef_out = outputs_path + "vegtype_EF_byClass/Vegtype_EF%s.%s.csv" % (cid, scen_name)
        vegtype_ldf_out = outputs_path + "vegtype_LDF_byClass/Vegtype_LDF%s.%s.csv" % (cid, scen_name)

        print("ClassID: %s" % cid)

        # Remove old base database (if exists)
        try:
            os.remove(M3VTEF_database_path)
            os.remove(M3LDF_database_path)
        except OSError:
            pass

        # Create M3VTEF database and initialize functions
        M3VTEF_connection = sqlite3.connect(M3VTEF_database_path)
        M3VTEF_connection.text_factory = str

        # Initalize SQL functions
        dbf.initDB(M3VTEF_connection)

        # Load CSV files as tables into the M3VTEF database
        mvtef.make_M3VTEF_tables(M3VTEF_connection, VegDB_tables, Description_Vegetation, DB_SLA, DB_emissions)

        # Run base database driver
        mvtef.run_M3VTEF_DB(M3VTEF_connection, vegtype_ef_out, cid, Jrating)

        print("M3VTEF Database created: %s\n" % M3VTEF_database_path)

        # Close database connection
        M3VTEF_connection.close()
        print("\n M3VTEF Database connection closed\n\n")

        # Create M3LDF database and initialize functions
        M3LDF_connection = sqlite3.connect(M3LDF_database_path)
        M3LDF_connection.text_factory = str

        # Initalize SQL functions
        dbf.initDB(M3LDF_connection)

        # Load CSV files as tables into the M3VTEF database
        mldf.make_M3LDF_tables(M3LDF_connection, VegDB_tables, Description_Vegetation, DB_LDF)

        # Run base database driver
        mldf.run_M3LDF_DB(M3LDF_connection, vegtype_ldf_out, cid, Jrating)

        print("M3LDF Database created: %s\n" % M3LDF_database_path)

        # Close database connection
        M3LDF_connection.close()
        print("\n M3LDF Database connection closed\n\n")

    print("Concatenating Vegtype EF and Vegtype LDF Tables")
    M3VTEF_connection = sqlite3.connect(M3VTEF_database_path)
    M3VTEF_connection.text_factory = str

    M3LDF_connection = sqlite3.connect(M3LDF_database_path)
    M3LDF_connection.text_factory = str

    for cid in classids:
        vegtype_ef_out = outputs_path + "vegtype_EF_byClass/Vegtype_EF%s.%s.csv" % (cid, scen_name)
        # Load each vegtype_ef csv output back into database
        mvtef.make_vegtyp_ef_tables(M3VTEF_connection, vegtype_ef_out, cid)
        vegtype_ldf_out = outputs_path + "vegtype_LDF_byClass/Vegtype_LDF%s.%s.csv" % (cid, scen_name)
        # Load each vegtype_ef csv output back into database
        mldf.make_vegtyp_ldf_tables(M3LDF_connection, vegtype_ldf_out, cid)

    Vegtype_EF_All_df = mvtef.concat_vegtype_ef_tables(M3VTEF_connection)
    Vegtype_EF_All_df.to_csv(outputs_path+"Vegtype_EF.%s.csv" % scen_name, index=False)
    M3VTEF_connection.close()

    Vegtype_LDF_All_df = mldf.concat_vegtype_ldf_tables(M3LDF_connection)
    Vegtype_LDF_All_df.to_csv(outputs_path+"Vegtype_LDF.%s.csv" % scen_name, index=False)
    M3LDF_connection.close()

    print("Vegtype EF and LDF Concatenated Output:")
    print(outputs_path+"Vegtype_EF.%s.csv" % scen_name)
    print(outputs_path+"Vegtype_LDF.%s.csv" % scen_name)
    print("\n\n\n\n\n\n\n\n\n")

    print("Creating grid EF and grid LDF databases")
    for cid in classids:
        vegtype_ef_input = outputs_path+"vegtype_EF_byClass/Vegtype_EF%s.%s.csv" % (cid, scen_name)
        grid_EF_output = outputs_path+"grid_EF_byClass/grid_EF%s.%s.csv" % (cid, scen_name)

        vegtype_ldf_input = outputs_path+"vegtype_LDF_byClass/Vegtype_LDF%s.%s.csv" % (cid, scen_name)
        grid_LDF_output = outputs_path+"grid_LDF_byClass/grid_LDF%s.%s.csv" % (cid, scen_name)

        print("ClassID: %s" % cid)

        # Remove old M3GEFP and M3GLDF databases (if exists)
        try:
            os.remove(M3GEFP_database_path)
            os.remove(M3GLDF_database_path)
        except OSError:
            pass

        # Create base database and initialize functions
        M3GEFP_connection = sqlite3.connect(M3GEFP_database_path)
        M3GEFP_connection.text_factory = str

        # Load CSV files as tables into the base database
        mgefp.make_M3GEFP_tables(M3GEFP_connection, GridDB_tables, Ecotype_Crop, Ecotype_Shrub,
                                 Ecotype_Herb, Ecotype_Tree, grid_ecotype, grid_growth_form,vegtype_ef_input)

        # Run M3GEF database driver
        mgefp.run_M3GEFP_DB(M3GEFP_connection, grid_EF_output)
        print("M3GEFP Database created: %s" % M3GEFP_database_path)

        # Close database connection
        M3GEFP_connection.close()
        print("\n M3GEFP Database connection closed\n\n")

        # Create base database and initialize functions
        M3GLDF_connection = sqlite3.connect(M3GLDF_database_path)
        M3GLDF_connection.text_factory = str

        # Load CSV files as tables into the base database
        mgldf.make_M3GLDF_tables(M3GLDF_connection, GridDB_tables, Ecotype_Crop, Ecotype_Shrub, 
			         Ecotype_Herb, Ecotype_Tree, grid_ecotype, grid_growth_form,vegtype_ldf_input)

        # Run M3GEF database driver
        mgldf.run_M3GLDF_DB(M3GLDF_connection, grid_LDF_output)
        print("M3GLDF Database created: %s" % M3GLDF_database_path)

        # Close database connection
        M3GLDF_connection.close()
        print("\n M3GLDF Database connection closed\n\n")

    M3GEFP_connection = sqlite3.connect(M3GEFP_database_path)
    M3GEFP_connection.text_factory = str
    M3GLDF_connection = sqlite3.connect(M3GLDF_database_path)
    M3GLDF_connection.text_factory = str
    print("Concatenating grid EF and grid LDF Tables")

    for cid in classids:
        grid_EF_out = outputs_path+"grid_EF_byClass/grid_EF%s.%s.csv" % (cid, scen_name)
        # Load each vegtype_ef csv output back into database
        mgefp.make_grid_ef_tables(M3GEFP_connection, grid_EF_out, cid)
        grid_LDF_out = outputs_path+"grid_LDF_byClass/grid_LDF%s.%s.csv" % (cid, scen_name)
        # Load each vegtype_ef csv output back into database
        mgldf.make_grid_ldf_tables(M3GLDF_connection, grid_LDF_out, cid)

    grid_EF_all_df = mgefp.concat_grid_ef_tables(M3GEFP_connection)
    grid_EF_all_df.to_csv(outputs_path+"grid_EF.%s.csv" % scen_name, index=False)
    grid_LDF_all_df = mgldf.concat_grid_ldf_tables(M3GLDF_connection)
    grid_LDF_all_df.to_csv(outputs_path+"grid_LDF.%s.csv" % scen_name, index=False)
    print("grid EF and grid LDF Concatenated Output:")
    print(outputs_path+"grid_EF.%s.csv" % scen_name)
    print(outputs_path+"grid_LDF.%s.csv" % scen_name)
    M3GEFP_connection.close()
