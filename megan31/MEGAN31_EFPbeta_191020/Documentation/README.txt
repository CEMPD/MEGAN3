README

Alex Guenther, alex.guenther@uci.edu
January 18, 2019

MEGAN Emission Factor Processor (MEGAN EFP)

1. INTRODUCTION
2. INSTALL
3. USAGE
4. NOTES ON MODIFYING PROGRAM AND PRIMARY DATABASE TABLES


~~~ INTRODUCTION ~~~

The MEGAN Emission Factor Processor (MEGAN EFP) is a program to calculate vegetative specific emissions factors, vegetative light dependent fractions, and grid average emission factors based on measurements and known biogenic emission factors. This program is built on Python (version 2.7) and SQLite. 

The following scripts and modules should be included:

MEGAN_EFP.py: 
Main script user specifies all paths and options in.

Included in ./src directory

run_M3EFP.py: Contains function that loops through classes and drives each database
db_functions.py: Module containing SQL functions and initialization function
M3VTEF.py: Module containing SQL queries to create M3VTEF DB
M3GEFP.py: Module containing SQL queries to create M3GEFP DB
M3LDF.py: Module containing SQL queries to create M3LDF DB
M3GLDF.py: Module containing SQL queries to create M3GLDF DB


Included in ./database directory

M3VTEF_DB_Description_README.txt
M3GEFP_DB_Description_README.txt

*Note the M3LDF and M3GLDF databases are structured in a similar way

Included in ./inputs folder

VegDB_Tables and GridDB_Tables directories containing CSV files that will be loaded into the respective database as individual tables









~~~ INSTALL ~~~

The MEGAN EFP program is built on the Python programming language. If your computer does not already have Python installed, we recommend installing the Anaconda distribution of Python via

https://www.anaconda.com/download/

Alternatively, Python can also be installed from source via

https://www.python.org/downloads

This program has been tested and successfully run with Python 2.7 on Windows 10, Centos 6 and Mac OS 10.10.

Several Python libraries and packages are required to run the MEGAN EFP program. Most are standard Python libraries that do not require any additional installation. All other need packages are available for download and installation via pip or conda. Included in this folder is a requirements.txt file that can be run as instructed below to install all additional packages 

Required Python Libraries and Packages

       pandas (Not included in standard Python build)
       numpy
  
To install any packages not included in the standard Python build, users can install them via pip or conda.

To install via pip run:
$ pip install -r requirements.txt

To install via conda:
$ conda install --file requirements.txt



~~~ USAGE ~~~

To run the MEGAN EFP program, verify user options and directory paths have been correctly configured in the MEGAN_EFP.py script. After verification, simply call from terminal:

$ python MEGAN_EFP.py

Running this program will generate 3 SQLite databases (M3VTEF, M3 LDF, and M3GEFP) based on the primary tables loaded from CSVs in the ./input directory and queried tables based on SQL queries run in the program. The program loops through each class overwriting each database each time while generating CSV output found in ./output/*byClass/. This output is then reloaded into the database when each class has been run and concatenated into a single table in the database and then that is written out as a CSV in ./outout/


~~~ NOTES ON MODIFYING PROGRAM ~~~

If you would like to make modifications to the program, please note the following:


a. General Design of Program

See included PowerPoint flowchart: M3EFP_DesignDocument.pptx 


b. Adding/Modifying SQL queries

All SQL queries are stored individually as functions within its respective sub module and executed in consecutive order by the run_*_DB function located at the bottom of the sub module. 

c. Adding/Modifying functions used in SQL queries

A number of functions are used in SQL queries to manipulate data. These are stored in the db_functions.py sub module. Before the functions can be added to the SQL queries, they must first be initialized in the initDB function in the db_functions.py module. Format for initializing a function is as follows:

conn.create_function("functionNameInSQL", #NumberOfInputs, functionNameInPython)


d. Modification and Visualization of Database

Users may find it helpful to examine the original database or the user created database directly. Beyond using SQLite command line tools, several third party programs exist to examine SQLite databases including the DB Browser for SQLite:

http://sqlitebrowser.org/


e. Using Custom Inputs/Measurements

Users may want to use custom measurements or modify the primary tables the program uses to make calculations. To make adjustments to the input tables, simply modify the CSVs that are loaded into each respective database in the ./inputs directory.

*NOTE: If the user changes the names of the CSV files, they will need to accordingly update the expected name in its respective sub module under the make_*_tables function

