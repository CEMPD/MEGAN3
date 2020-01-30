#!/usr/bin/env python
"""

Module that contains functions related to createing the Vegtype EF database

Functions include
    run_M3LDF_DB: Executes the functions and queries to generate the database output for each class
    make_M3LDF_tables: Create tables based on input CSVs
    LDF_*: SQL queries to run the database for the given vegetative class
    make_grid_ldf_tables: Load the individual classes back into the database
    concat_grid_ldf_tables: Concatenate the individual classes into 1 table

"""
import pandas as pd


def make_M3LDF_tables(conn, csv_input_dir, Description_Vegetation, DB_LDF):

    """

    Function to load CSV inputs to create primary tables of the M3LDF database
    :param conn: SQLite database connection
    :param csv_input_dir: Directory of CSV inputs
    :param Description_Vegetation: input csv file name
    :param DB_LDF: input csv file name
    :return: A SQLite database with primary tables

    """

    print("Creating M3LDF DB tables from: %s\n" % csv_input_dir)
    csv_Canopy_Position_adjust = pd.read_csv(csv_input_dir + 'Canopy_Position_adjust.csv')
    csv_Description_Class = pd.read_csv(csv_input_dir + 'Description_Class.csv')
    csv_Description_Vegetation = pd.read_csv(csv_input_dir + Description_Vegetation)
    csv_Description_Compounds = pd.read_csv(csv_input_dir + 'Description_Compounds.csv')
    csv_Description_References = pd.read_csv(csv_input_dir + 'Description_References.csv')
    csv_DB_LDF = pd.read_csv(csv_input_dir + DB_LDF)

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

    csv_DB_LDF.to_sql("DB LDF", conn, flavor='sqlite', if_exists='replace')
    print("'DB LDF' Table Loaded from: %s" % csv_input_dir + DB_LDF)


def LDF_VEGID_Jrating(conn):
    c = conn.cursor()
    c.execute("CREATE TABLE 'LDF VEGID Jrating' AS SELECT 'DB LDF'.[VegID] AS [VegID], "
              "'DB LDF'.[CAS] AS [CAS], "
              "MAX('DB LDF'.[J-rating]) AS MaxJ, "
              "MIN('DB LDF'.[J-rating]) AS MinJ, "
              "COUNT('DB LDF'.[J-rating]) AS [CountOfJ-rating] "
              "FROM 'DB LDF' GROUP BY 'DB LDF'.[VegID], 'DB LDF'.[CAS];")

    print("'LDF VEGID Jrating' Table Created")


def LDF_VEGID(conn, ClassID, Jrating):
    c = conn.cursor()
    c.execute("CREATE TABLE 'LDF VEGID' AS SELECT 'Description Vegetation'.[VegID] AS [VegID], "
              "'Description Vegetation'.[genus] AS [genus], "
              "'Description Vegetation'.[Family] AS [Family], "
              "'Description Vegetation'.[growthform] AS [growthform], "
              "'DB LDF'.[J-rating] AS [J-rating], "
              "AVG([LDfrac]) AS VEGIDLDF "
              "FROM 'Canopy Position adjust' "
              "INNER JOIN ('LDF VEGID Jrating' "
              "INNER JOIN (('Description Vegetation' "
              "INNER JOIN ('DB LDF' "
              "INNER JOIN 'Description Compounds' ON 'DB LDF'.[CAS] = 'Description Compounds'.[CAS]) "
              "ON 'Description Vegetation'.[VegID] = 'DB LDF'.[VegID])) "
              "ON ('LDF VEGID Jrating'.[MaxJ] = 'DB LDF'.[J-rating]) "
              "AND ('LDF VEGID Jrating'.[CAS] = 'DB LDF'.[CAS]) "
              "AND ('LDF VEGID Jrating'.[VegID] = 'DB LDF'.[VegID])) "
              "WHERE ((('Description Compounds'.[ClassPlusID])=?) AND (('DB LDF'.[J-rating])>=?)) "
              "GROUP BY 'Description Vegetation'.[VegID], "
              "'Description Vegetation'.[genus], "
              "'Description Vegetation'.[Family], "
              "'Description Vegetation'.[growthform], "
              "'DB LDF'.[J-rating];", (ClassID, Jrating,))


def LDF_genus(conn):
    c = conn.cursor()
    c.execute("CREATE TABLE 'LDF genus' AS SELECT 'LDF VEGID'.[genus] AS [genus], AVG('LDF VEGID'.[VEGIDLDF]) "
              "AS genusLDF FROM 'LDF VEGID' GROUP BY 'LDF VEGID'.[genus];")

    print("'LDF genus' Table Created")


def LDF_family(conn):
    c = conn.cursor()
    c.execute("CREATE TABLE 'LDF family' AS SELECT 'LDF VEGID'.[Family] AS [Family], AVG('LDF VEGID'.[VEGIDLDF]) "
              "AS [Family LDF] FROM 'LDF VEGID' GROUP BY 'LDF VEGID'.[Family];")

    print("'LDF family' Table Created")


def LDF_growthform(conn):
    c = conn.cursor()
    c.execute("CREATE TABLE 'LDF growthform' AS SELECT 'LDF VEGID'.[growthform] AS [growthform], "
              "AVG('LDF VEGID'.[VEGIDLDF]) AS [growthformLDF] FROM 'LDF VEGID' GROUP BY 'LDF VEGID'.[growthform];")

    print("'LDF growthform' Table Created")


def Vegtype_LDF(conn):
    c = conn.cursor()
    c.execute("CREATE TABLE 'Vegtype LDF' AS SELECT 'Description Vegetation'.[VegID] AS [VegID], "
              "SelectLDF([VEGIDLDF],[genusLDF],[Family LDF],[growthformLDF]) AS VegLDF "
              "FROM ((('Description Vegetation' LEFT JOIN 'LDF Family' "
              "ON 'Description Vegetation'.[Family] = 'LDF Family'.[Family]) "
              "LEFT JOIN 'LDF genus' ON 'Description Vegetation'.[genus] = 'LDF genus'.[genus]) "
              "LEFT JOIN 'LDF VEGID' ON 'Description Vegetation'.[VegID] = 'LDF VEGID'.[VegID]) "
              "INNER JOIN 'LDF growthform' ON 'Description Vegetation'.[growthform] = 'LDF growthform'.[growthform];")

    print("'Vegtype LDF' Table Created\n")

    Vegtype_EF_query_df = pd.read_sql_query("SELECT * FROM 'Vegtype LDF';", conn)

    return Vegtype_EF_query_df


def make_vegtyp_ldf_tables(conn, Vegtype_LDF_csv_output, ClassID):
    classID_name = "Vegtype LDF%s" % ClassID
    csv_Vegtype_LDF = pd.read_csv(Vegtype_LDF_csv_output)
    csv_Vegtype_LDF.to_sql(classID_name, conn, flavor='sqlite', if_exists='replace')
    print("'Vegtype LDF%s' Table Loaded from: %s" % (ClassID, Vegtype_LDF_csv_output))


def concat_vegtype_ldf_tables(conn):
    c = conn.cursor()
    c.execute("CREATE TABLE 'Vegtype LDF All' AS SELECT DISTINCT 'Description Vegetation'.[VegID] AS [VegID], "
              "'Vegtype LDF1'.[VegLDF] AS [VegLDF1], "
              "'Vegtype LDF2'.[VegLDF] AS [VegLDF2], "
              "'Vegtype LDF3'.[VegLDF] AS [VegLDF3], "
              "'Vegtype LDF4'.[VegLDF] AS [VegLDF4], "
              "'Vegtype LDF5'.[VegLDF] AS [VegLDF5], "
              "'Vegtype LDF6'.[VegLDF] AS [VegLDF6], "
              "'Vegtype LDF7'.[VegLDF] AS [VegLDF7], "
              "'Vegtype LDF8'.[VegLDF] AS [VegLDF8], "
              "'Vegtype LDF9'.[VegLDF] AS [VegLDF9], "
              "'Vegtype LDF10'.[VegLDF] AS [VegLDF10], "
              "'Vegtype LDF11'.[VegLDF] AS [VegLDF11], "
              "'Vegtype LDF12'.[VegLDF] AS [VegLDF12], "
              "'Vegtype LDF13'.[VegLDF] AS [VegLDF13], "
              "'Vegtype LDF14'.[VegLDF] AS [VegLDF14], "
              "'Vegtype LDF15'.[VegLDF] AS [VegLDF15], "
              "'Vegtype LDF16'.[VegLDF] AS [VegLDF16], "
              "'Vegtype LDF17'.[VegLDF] AS [VegLDF17], "
              "'Vegtype LDF18'.[VegLDF] AS [VegLDF18], "
              "'Vegtype LDF19'.[VegLDF] AS [VegLDF19], "
              "'Vegtype LDF20'.[VegLDF] AS [VegLDF20] "
              "FROM 'Description Vegetation' "
              "LEFT JOIN 'Vegtype LDF1' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF1'.[VegID] "
              "LEFT JOIN 'Vegtype LDF2' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF2'.[VegID] "
              "LEFT JOIN 'Vegtype LDF3' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF3'.[VegID] "
              "LEFT JOIN 'Vegtype LDF4' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF4'.[VegID] "
              "LEFT JOIN 'Vegtype LDF5' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF5'.[VegID] "
              "LEFT JOIN 'Vegtype LDF6' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF6'.[VegID] "
              "LEFT JOIN 'Vegtype LDF7' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF7'.[VegID] "
              "LEFT JOIN 'Vegtype LDF8' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF8'.[VegID] "
              "LEFT JOIN 'Vegtype LDF9' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF9'.[VegID] "
              "LEFT JOIN 'Vegtype LDF10' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF10'.[VegID] "
              "LEFT JOIN 'Vegtype LDF11' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF11'.[VegID] "
              "LEFT JOIN 'Vegtype LDF12' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF12'.[VegID] "
              "LEFT JOIN 'Vegtype LDF13' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF13'.[VegID] "
              "LEFT JOIN 'Vegtype LDF14' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF14'.[VegID] "
              "LEFT JOIN 'Vegtype LDF15' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF15'.[VegID] "
              "LEFT JOIN 'Vegtype LDF16' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF16'.[VegID] "
              "LEFT JOIN 'Vegtype LDF17' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF17'.[VegID] "
              "LEFT JOIN 'Vegtype LDF18' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF18'.[VegID] "
              "LEFT JOIN 'Vegtype LDF19' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF19'.[VegID] "
              "LEFT JOIN 'Vegtype LDF20' ON 'Description Vegetation'.[VegID] = 'Vegtype LDF20'.[VegID];")

    Vegtype_LDF_All_query_df = pd.read_sql_query("SELECT * FROM 'Vegtype LDF All';", conn)

    return Vegtype_LDF_All_query_df


def run_M3LDF_DB(db_connection, output_csv_path, ClassID, Jrating):
    LDF_VEGID_Jrating(db_connection)
    LDF_VEGID(db_connection, ClassID, Jrating)
    LDF_genus(db_connection)
    LDF_family(db_connection)
    LDF_growthform(db_connection)
    query_df = Vegtype_LDF(db_connection)
    query_df.to_csv(output_csv_path, index=False)  # Save LDF table to CSV
    print("Vegtype LDF Table CSV Generated: %s" % output_csv_path)
