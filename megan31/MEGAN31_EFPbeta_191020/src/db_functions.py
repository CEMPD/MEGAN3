#!/usr/bin/env python

"""

This module contains functions that are used in SQL queries.
These functions are intializaed when the initDB function is called with the database connection

"""

# ~~~~~~~~~~~~~ BEGIN SQL FUNCTIONS ~~~~~~~~~~~~~ #
def EFarea(SLA, MW, EFuggh, EFnmm2s):
    """Calculate emission factor area
    :param SLA: Specific leaf area
    :param MW: Molecular weight of compound (grams per mole)
    :param EFuggh: Emission rate in micrograms compound per gram dry weight foliage per hour
    :param EFnmm2s: Emission rate in nanomoles compound per square meter foliage area (one sided) per second
    :return: Calculated EF area or None if values are not available from inputs
    """
    try:
        EFnmm2s = float(EFnmm2s)
        if EFnmm2s > -9:
            EFarea_out = EFnmm2s  # already have per mass EF
        else:
            SLA = float(SLA *.0001)  #convert cm2/g to m2/g
            MW = float(MW)
            EFuggh = float(EFuggh)
            uh2ns = 3.6  # convert micrograms to nanograms and hours to seconds
            EFarea_out = (EFuggh) / (SLA * uh2ns * MW)  # convert using specific leaf area, molecular weight

    except:  # Except error is broad to catch anything that cant be converted to float dtype
        EFarea_out = None

    return EFarea_out


def SelectEF(vegEF, genEF, FamEF, GFEF):
    """Select most representative Emission Factor
    :param vegEF: Vegatative EF
    :param genEF: Genus EF
    :param FamEF: Family EF
    :param GFEF: Growth Form EF
    :return: Most representative EF
    """
    if vegEF > -9:
        SelectEF_out = vegEF
    elif genEF > -9:
        SelectEF_out = genEF
    elif FamEF > -9:
        SelectEF_out = FamEF
    else:
        SelectEF_out = GFEF

    return SelectEF_out


def SelectSLA(vegSLA, genSLA, FamSLA, GFSLA):
    """Select most representative plant specific leaf area (SLA)
    :param vegSLA: Vegetative SLA
    :param genSLA: Genus SLA
    :param FamSLA: Family SLA
    :param GFSLA: growth form
    :return: Most representative SLA
    """
    if vegSLA > 0:
        SelectSLA_out = vegSLA
    elif genSLA > 0:
        SelectSLA_out = genSLA
    elif FamSLA > 0:
        SelectSLA_out = FamSLA
    else:
        SelectSLA_out = GFSLA

    return SelectSLA_out

def SelectLDF(vegLDF, genLDF, FamLDF, GFLDF):
    """Select most representative LDF
    :param vegEF: Vegatative LDF
    :param genEF: Genus LDF
    :param FamEF: Family LDF
    :param GFEF: Growth Form LDF
    :return: Most representative LDF
    """
    if vegLDF > -9:
        SelectLDF_out = vegLDF
    elif genLDF > -9:
        SelectLDF_out = genLDF
    elif FamLDF > -9:
        SelectLDF_out = FamLDF
    else:
        SelectLDF_out = GFLDF

    return SelectLDF_out



# ~~~~~~~~~~~~~ END SQL FUNCTIONS ~~~~~~~~~~~~~ #


# ~~~~~~~~~~~~~ DATABASE INITIALIZATION FUNCTION ~~~~~~~~~~~~~ #

def initDB(conn):
    """
    Function used to initialize the 3 custom functions
     used in subsequent SQL queries
    :param conn: database connection
    :return: 3 function that can be used in SQL queries
    """
    print("\n Database Initialization")
    # Connect to the SQLite database

    # Initialize functions used in queries
    conn.create_function("EFarea", 4, EFarea)
    conn.create_function("SelectEF", 4, SelectEF)
    conn.create_function("SelectSLA", 4, SelectSLA)
    conn.create_function("SelectLDF", 4, SelectLDF)
    print(" Initialization Complete\n")
