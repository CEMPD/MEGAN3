MEGAN3 Grid Emission Factor Processor (M3GEFP) 

This processor integrates vegetation-type emission factors with landcover 
(growth form and ecotype fractions) and vegetation-type speciation data to estimate
grid average emission factors
   
INPUT FILES
------------
Emission Factor files (input from M3VEFP database)
E1) Vegtype EF01: vegetation type emission factors (nmol compound per m2 leaf area per second) for emission class 1 (Isoprene)
	VegID: unique number identifying each vegetation type
	VegEF01: Isoprene emission factor (nmol isoprene per m2 leaf area per second) for specified vegetation type
E2) Additional files can be created for each emission class

Ecotype Speciation files
S1) Ecotype Tree Speciation: vegetation type composition of trees in each ecotype
	EcotypeID: unique number identifying each ecotype
	VegID: unique number identifying each vegetation type
	TreeSpecFrac: Fraction of the trees in specified ecotype comprised of the specified vegetation type
S2) Ecotype Shrub Speciation: vegetation type composition of shrubs in each ecotype
	EcotypeID: unique number identifying each ecotype
	VegID: unique number identifying each vegetation type
	ShrubSpecFrac: Fraction of the shrubs in specified ecotype comprised of the specified vegetation type
S3) Ecotype Herb Speciation: vegetation type composition of herb in each ecotype
	EcotypeID: unique number identifying each ecotype
	VegID: unique number identifying each vegetation type
	HerbSpecFrac: Fraction of the herb in specified ecotype comprised of the specified vegetation type
S4) Ecotype Crop Speciation: vegetation type composition of crops in each ecotype
	EcotypeID: unique number identifying each ecotype
	VegID: unique number identifying each vegetation type
	CropSpecFrac: Fraction of the crops in specified ecotype comprised of the specified vegetation type

Grid Files
G1) grid growth form: Fraction of grid surface area covered by each major growth form
	gridID: unique number identifying each grid location in a modeling domain
	TreeFrac: Fraction of grid covered by trees
	CropFrac: Fraction of grid covered by crops
	ShrubFrac: Fraction of grid covered by shrubs
	HerbFrac: Fraction of grid covered by herbs

G2) grid ecotype: Fraction of grid surface area covered by each ecotype
	gridID: unique number identifying each grid location in a modeling domain
	EcotypeID: unique number identifying each ecotype
	EcoTypeFrac: Fraction of grid covered by specified ecotype

OUTPUT FILES
------------
O1) Grid EF01: Grid average emission factors (nmol compound per m2 leaf area per second) for emission class 1 (isoprene)
	gridID: unique number identifying each grid location in a modeling domain
	EF01: grid average emission factor (nmol isoprene per m2 leaf area per second) for emission class 1 (isoprene)
O2) Additional files can be created for each emission class
