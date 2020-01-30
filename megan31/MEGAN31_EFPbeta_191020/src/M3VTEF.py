#!/usr/bin/env python
"""

Module that contains functions related to creating the Vegtype EF database

Functions include
    run_M3VTEF_DB: Executes the functions and queries to generate the database output for each class
    make_M3VTEF_tables: Create tables based on input CSVs
    SLA_* & EF_*: SQL queries to run the database for the given vegetative class
    make_grid_ef_tables: Load the individual classes back into the database
    concat_grid_ef_tables: Concatenate the individual classes into 1 table

"""

# Import libraries, packages, and modules
from __future__ import division
import warnings
import pandas as pd
warnings.simplefilter(action="ignore", category=FutureWarning)


def make_M3VTEF_tables(conn, csv_input_dir, Description_Vegetation, DB_SLA, DB_emissions):

    """

    Function to load CSV inputs to create primary tables of the M3VTEF database
    :param conn: SQLite database connection
    :param csv_input_dir: Directory of CSV inputs
    :param Description_Vegetation: input csv file name
    :param DB_SLA: input csv file name
    :param DB_emissions: input csv file name
    :return: A SQLite database with primary tables

    """

    print("Creating M3VTEF DB tables from: %s\n" % csv_input_dir)
    csv_Canopy_Position_adjust = pd.read_csv(csv_input_dir + 'Canopy_Position_adjust.csv')
    csv_Description_Class = pd.read_csv(csv_input_dir + 'Description_Class.csv')
    csv_Description_Vegetation = pd.read_csv(csv_input_dir + Description_Vegetation)
    csv_Description_Compounds = pd.read_csv(csv_input_dir + 'Description_Compounds.csv')
    csv_Description_References = pd.read_csv(csv_input_dir + 'Description_References.csv')
    csv_DB_SLA = pd.read_csv(csv_input_dir + DB_SLA)
    csv_DB_emissions = pd.read_csv(csv_input_dir + DB_emissions)

    csv_Canopy_Position_adjust.to_sql("Canopy Position adjust", conn, flavor='sqlite', if_exists='replace')
    print("'Canopy Position adjust' Table Loaded from: %s" % csv_input_dir + 'Canopy_Position_adjust.csv')

    csv_Description_Class.to_sql("Description Class", conn, flavor='sqlite', if_exists='replace')
    print("'Description Class' Table Loaded from: %s" % csv_input_dir + 'Description_Class.csv')

    csv_Description_Vegetation.to_sql("Description Vegetation", conn, flavor='sqlite', if_exists='replace')
    print("'Description Vegetation' Table Loaded from: %s" % csv_input_dir + Description_Vegetation)

    csv_Description_Compounds.to_sql("Description Compounds", conn, flavor='sqlite', if_exists='replace')
    print("'Description Compounds' Table Loaded from: %s" % csv_input_dir + 'Description_Compounds.csv')

    csv_Description_References.to_sql("Description References", conn, flavor='sqlite', if_exists='replace')
    print("'Description References' Table Loaded from: %s" % csv_input_dir + 'Description_References.csv')

    csv_DB_SLA.to_sql("DB SLA", conn, flavor='sqlite', if_exists='replace')
    print("'DB SLA' Table Loaded from: %s" % csv_input_dir + DB_SLA)

    csv_DB_emissions.to_sql("DB emissions", conn, flavor='sqlite', if_exists='replace')
    print("'DB emissions' Table Loaded from: %s" % csv_input_dir + DB_emissions)


def SLA_vegID(conn, Jrating):
    c = conn.cursor()
    c.execute("CREATE TABLE 'SLA vegID' AS SELECT 'Description Vegetation'.[VegID] AS [VegID], "
              "'Description Vegetation'.[genus] AS [genus], 'Description Vegetation'.[Family] as [Family], "
              "'Description Vegetation'.[GrowthForm] AS [GrowthForm], "
              "AVG([SLAcm2g]*[SLAadj]) AS [VEGIDSLA], "
              "'DB SLA'.[J-rating] FROM 'Canopy Position adjust' "
              "INNER JOIN ('DB SLA' INNER JOIN 'Description Vegetation' "
              "ON 'DB SLA'.[VegID] = 'Description Vegetation'.VegID COLLATE NOCASE) "
              "ON 'Canopy Position adjust'.[Position] = 'DB SLA'.[Position] COLLATE NOCASE "
              "WHERE ((('DB SLA'.[J-rating])>=?));", (Jrating,))

    print("'SLA vegID' Table Created")


def SLA_genus(conn):
    c = conn.cursor()
    c.execute("CREATE TABLE 'SLA genus' AS SELECT 'SLA vegID'.[genus] AS [genus], "
              "AVG('SLA vegID'.[VEGIDSLA]) AS genusSLA FROM 'SLA vegID' GROUP BY 'SLA vegID'.[genus];")

    print("'SLA genus' Table Created")


def SLA_family(conn):
    c = conn.cursor()
    c.execute(
        "CREATE TABLE 'SLA family' AS SELECT 'SLA vegID'.[Family] AS [Family], AVG('SLA vegID'.[VEGIDSLA]) AS familySLA "
        "FROM 'SLA vegID' GROUP BY 'SLA vegID'.[Family];")

    print("'SLA family' Table Created")


def SLA_growthform(conn):
    c = conn.cursor()
    c.execute("CREATE TABLE 'SLA growthform' AS SELECT 'SLA vegID'.[growthform] AS [growthform], "
              "AVG('SLA vegID'.[VEGIDSLA]) AS growthformSLA "
              "FROM 'SLA vegID' GROUP BY 'SLA vegID'.[growthform];")

    print("'SLA growthform' Table Created")


def SLA(conn):
    c = conn.cursor()
    c.execute("CREATE TABLE 'SLA' AS SELECT 'Description Vegetation'.[VegID], "
              "SelectSLA([VEGIDSLA],[genusSLA],[FamilySLA],[growthformSLA]) AS SLA "
              "FROM ((('Description Vegetation' "
              "LEFT JOIN  'SLA vegID' ON 'Description Vegetation'.[VegID] = 'SLA vegID'.[VegID] COLLATE NOCASE) "
              "LEFT JOIN 'SLA genus' ON 'Description Vegetation'.[genus] = 'SLA genus'.[genus] COLLATE NOCASE)"
              "LEFT JOIN 'SLA family' ON 'Description Vegetation'.[Family] = 'SLA family'.[Family] COLLATE NOCASE)"
              "LEFT JOIN 'SLA growthform' ON 'Description Vegetation'.[growthform] = 'SLA growthform'.[growthform] COLLATE NOCASE;")

    print("'SLA' Table Created")


def EF_VEGID_Jrating(conn):
    c = conn.cursor()
    c.execute("CREATE TABLE 'EF VEGID Jrating' AS SELECT 'DB emissions'.[VegID] AS [VegID], "
              "'DB emissions'.[CAS] AS [CAS], "
              "MAX('DB emissions'.[J-rating]) AS MaxJ, "
              "MIN('DB emissions'.[J-rating]) AS MinJ, "
              "COUNT('DB emissions'.[J-rating]) AS [CountOfJ-rating] "
              "FROM 'DB emissions' GROUP BY 'DB emissions'.[VegID], 'DB emissions'.[CAS];")

    print("'EF VEGID Jrating' Table Created")


def EF_VEGID(conn, ClassID, Jrating):
    c = conn.cursor()
    c.execute("CREATE TABLE 'EF VEGID' AS SELECT 'Description Vegetation'.[VegID] AS [VegID], "
              "'Description Vegetation'.[genus] AS [genus], "
              "'Description Vegetation'.[Family] AS [Family], "
              "'Description Vegetation'.[growthform] AS [growthform], "
              "AVG(EFarea([SLA],[MW],[ERuggh],[ERnmm2s])*[EFadj]) AS VEGIDEF "
              "FROM 'Canopy Position adjust' INNER JOIN ('EF VEGID Jrating' INNER JOIN (('Description Vegetation' "
              "INNER JOIN ('DB emissions' "
              "INNER JOIN 'Description Compounds' ON 'DB emissions'.[CAS] = 'Description Compounds'.[CAS]) "
              "ON 'Description Vegetation'.[VegID] = 'DB emissions'.[VegID]) "
              "INNER JOIN 'SLA' ON 'Description Vegetation'.[VegID] = 'SLA'.[VegID]) "
              "ON ('EF VEGID Jrating'.[MaxJ] = 'DB emissions'.[J-rating]) "
              "AND ('EF VEGID Jrating'.[CAS] = 'DB emissions'.[CAS]) "
              "AND ('EF VEGID Jrating'.[VegID] = 'DB emissions'.[VegID])) "
              "ON 'Canopy Position adjust'.[Position] = 'DB emissions'.[Position] "
              "WHERE ((('Description Compounds'.[ClassPlusID])=?) AND ('DB emissions'.[J-rating]>=?)) "
              "GROUP BY 'Description Vegetation'.[VegID], 'Description Vegetation'.[genus], "
              "'Description Vegetation'.[Family], 'Description Vegetation'.[growthform];", (ClassID, Jrating,))

    print("'EF VEGID' Table Created")


def EF_genus(conn):
    c = conn.cursor()
    c.execute("CREATE TABLE 'EF genus' AS SELECT 'EF VEGID'.[genus] AS [genus], AVG('EF VEGID'.[VEGIDEF]) "
              "AS genusEF FROM 'EF VEGID' GROUP BY 'EF VEGID'.[genus];")

    print("'EF genus' Table Created")


def EF_family(conn):
    c = conn.cursor()
    c.execute("CREATE TABLE 'EF family' AS SELECT 'EF VEGID'.[Family] AS [Family], AVG('EF VEGID'.[VEGIDEF]) "
              "AS [Family EF] FROM 'EF VEGID' GROUP BY 'EF VEGID'.[Family];")

    print("'EF family' Table Created")


def EF_growthform(conn):
    c = conn.cursor()
    c.execute("CREATE TABLE 'EF growthform' AS SELECT 'EF VEGID'.[growthform] AS [growthform], "
              "AVG('EF VEGID'.[VEGIDEF]) AS [growthformEF] FROM 'EF VEGID' GROUP BY 'EF VEGID'.[growthform];")

    print("'EF growthform' Table Created")


def Vegtype_EF(conn):
    c = conn.cursor()
    c.execute("CREATE TABLE 'Vegtype EF' AS SELECT 'Description Vegetation'.[VegID] AS [VegID], "
              "SelectEF([VEGIDEF],[genusEF],[Family EF],[growthformEF]) AS VegEF "
              "FROM ((('Description Vegetation' LEFT JOIN 'EF Family' "
              "ON 'Description Vegetation'.[Family] = 'EF Family'.[Family]) "
              "LEFT JOIN 'EF genus' ON 'Description Vegetation'.[genus] = 'EF genus'.[genus]) "
              "LEFT JOIN 'EF VEGID' ON 'Description Vegetation'.[VegID] = 'EF VEGID'.[VegID]) "
              "INNER JOIN 'EF growthform' ON 'Description Vegetation'.[growthform] = 'EF growthform'.[growthform];")

    print("'Vegtype EF' Table Created\n")

    Vegtype_EF_query_df = pd.read_sql_query("SELECT * FROM 'Vegtype EF';", conn)

    return Vegtype_EF_query_df


def make_vegtyp_ef_tables(conn, Vegtype_EF_csv_output, ClassID):
    classID_name = "Vegtype EF%s" % ClassID
    csv_Vegtype_EF = pd.read_csv(Vegtype_EF_csv_output)
    csv_Vegtype_EF.to_sql(classID_name, conn, flavor='sqlite', if_exists='replace')
    print("'Vegtype EF%s' Table Loaded from: %s" % (ClassID, Vegtype_EF_csv_output))


def concat_vegtype_ef_tables(conn):
    c = conn.cursor()
    c.execute("CREATE TABLE 'Vegtype EF All' AS SELECT DISTINCT 'SLA'.[VegID] AS [VegID], "
              "'Vegtype EF1'.[VegEF1] AS [VegEF1], "
              "'Vegtype EF2'.[VegEF2] AS [VegEF2], "
              "'Vegtype EF3'.[VegEF3] AS [VegEF3], "
              "'Vegtype EF4'.[VegEF4] AS [VegEF4], "
              "'Vegtype EF5'.[VegEF5] AS [VegEF5], "
              "'Vegtype EF6'.[VegEF6] AS [VegEF6], "
              "'Vegtype EF7'.[VegEF7] AS [VegEF7], "
              "'Vegtype EF8'.[VegEF8] AS [VegEF8], "
              "'Vegtype EF9'.[VegEF9] AS [VegEF9], "
              "'Vegtype EF10'.[VegEF10] AS [VegEF10], "
              "'Vegtype EF11'.[VegEF11] AS [VegEF11], "
              "'Vegtype EF12'.[VegEF12] AS [VegEF12], "
              "'Vegtype EF13'.[VegEF13] AS [VegEF13], "
              "'Vegtype EF14'.[VegEF14] AS [VegEF14], "
              "'Vegtype EF15'.[VegEF15] AS [VegEF15], "
              "'Vegtype EF16'.[VegEF16] AS [VegEF16], "
              "'Vegtype EF17'.[VegEF17] AS [VegEF17], "
              "'Vegtype EF18'.[VegEF18] AS [VegEF18], "
              "'Vegtype EF19'.[VegEF19] AS [VegEF19], "
              "'Vegtype EF20'.[VegEF20] AS [VegEF20] "
              "FROM 'SLA' "
              "LEFT JOIN 'Vegtype EF1' ON 'SLA'.[VegID] = 'Vegtype EF1'.[VegID] "
              "LEFT JOIN 'Vegtype EF2' ON 'SLA'.[VegID] = 'Vegtype EF2'.[VegID] "
              "LEFT JOIN 'Vegtype EF3' ON 'SLA'.[VegID] = 'Vegtype EF3'.[VegID] "
              "LEFT JOIN 'Vegtype EF4' ON 'SLA'.[VegID] = 'Vegtype EF4'.[VegID] "
              "LEFT JOIN 'Vegtype EF5' ON 'SLA'.[VegID] = 'Vegtype EF5'.[VegID] "
              "LEFT JOIN 'Vegtype EF6' ON 'SLA'.[VegID] = 'Vegtype EF6'.[VegID] "
              "LEFT JOIN 'Vegtype EF7' ON 'SLA'.[VegID] = 'Vegtype EF7'.[VegID] "
              "LEFT JOIN 'Vegtype EF8' ON 'SLA'.[VegID] = 'Vegtype EF8'.[VegID] "
              "LEFT JOIN 'Vegtype EF9' ON 'SLA'.[VegID] = 'Vegtype EF9'.[VegID] "
              "LEFT JOIN 'Vegtype EF10' ON 'SLA'.[VegID] = 'Vegtype EF10'.[VegID] "
              "LEFT JOIN 'Vegtype EF11' ON 'SLA'.[VegID] = 'Vegtype EF11'.[VegID] "
              "LEFT JOIN 'Vegtype EF12' ON 'SLA'.[VegID] = 'Vegtype EF12'.[VegID] "
              "LEFT JOIN 'Vegtype EF13' ON 'SLA'.[VegID] = 'Vegtype EF13'.[VegID] "
              "LEFT JOIN 'Vegtype EF14' ON 'SLA'.[VegID] = 'Vegtype EF14'.[VegID] "
              "LEFT JOIN 'Vegtype EF15' ON 'SLA'.[VegID] = 'Vegtype EF15'.[VegID] "
              "LEFT JOIN 'Vegtype EF16' ON 'SLA'.[VegID] = 'Vegtype EF16'.[VegID] "
              "LEFT JOIN 'Vegtype EF17' ON 'SLA'.[VegID] = 'Vegtype EF17'.[VegID] "
              "LEFT JOIN 'Vegtype EF18' ON 'SLA'.[VegID] = 'Vegtype EF18'.[VegID] "
              "LEFT JOIN 'Vegtype EF19' ON 'SLA'.[VegID] = 'Vegtype EF19'.[VegID] "
              "LEFT JOIN 'Vegtype EF20' ON 'SLA'.[VegID] = 'Vegtype EF20'.[VegID];")

    Vegtype_EF_All_query_df = pd.read_sql_query("SELECT * FROM 'Vegtype EF All';", conn)

    return Vegtype_EF_All_query_df


def run_M3VTEF_DB(db_connection, vegtype_ef_out, ClassID, Jrating):
    """
    M3VTEF creates the M3VTEF database and associated tables
    :param db_connection: Database connection
    :param output_csv_path: Path to CSV output
    :param ClassID: ClassID (user option)
    :param Jrating: Jrating (user option)
    :return: A SQLite database and copy of the Vegtype EF table as a CSV file
    """
    # Queries need to be executed in specific order
    # due to dependency of tables generated from previous queries
    print("\n BEGINNING M3VTEF DB GENERATION")
    SLA_vegID(db_connection, Jrating)
    SLA_genus(db_connection)
    SLA_family(db_connection)
    SLA_growthform(db_connection)
    SLA(db_connection)
    EF_VEGID_Jrating(db_connection)
    EF_VEGID(db_connection, ClassID, Jrating)
    EF_genus(db_connection)
    EF_family(db_connection)
    EF_growthform(db_connection)
    Vegtype_EF_output = Vegtype_EF(db_connection)

    Vegtype_EF_output.to_csv(vegtype_ef_out, index=False, header=['VegID', 'VegEF%s' % ClassID])  # Save Vegtype EF table to CSV
    print("Vegtype EF Table CSV Generated: %s" % vegtype_ef_out)
