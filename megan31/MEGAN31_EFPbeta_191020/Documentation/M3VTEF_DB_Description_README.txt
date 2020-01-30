MEGAN3 Vegetation Type EF DataBase (VTEFDB) input files
Description files
D1) Description Compounds: Description of emitted compounds included in MEGAN3
	CAS (Chemical Abstracts Service) registry number: unique identifier for each compound
	CompoundName: Name of chemical compound (or group of compounds)
	MW: Molecular weight of compound (grams per mole)
	ClassPlusID: Unique identifier for MEGAN3 classes plus some total EF categories (links to Description Class)
D2) Description Class: Description of MEGAN3 classes plus total EF categories 
	ClassPlusID: Unique identifier for classes plus total EF categories (links to "Description Compounds" and "Speciate EF categories")
	ClassPlusDescription: Description of MEG3 compound classes
D3) Speciate EF categories: Scheme to map summed categories in emissions database (e.g., total monoterpenes) to MEGAN3 classes
	ClassPlusID: links to "Description Class"
	ClassID: Unique identifier for MEGAN3 emission classes 
	SpeciateFrac:Fraction of total (e.g. total monoterpenes) assocaited with each MEGAN3 class
D4) Description References: contains citation and other information on publications containing emissions and SLM (specific leaf mass) observations 
	REFID: Unique identifier for each publication
	REFauthor: last name of first author of publication
	REFYear: year of publication
        Citation: Citation information 
D4) Description Vegetation: Information on all of the vegetation types
	VegID: Unique identifier for each vegetation type 
	Veg_Description: Description of vegetation type
	Common_Name: common name of vegetation type
	Genus: Genus of vegetation type ("NA" if vegetation type is not a specific species or genus)
	Family: Family of vegetation type ("NA" if vegetation type is not a specific species or genus)
	Division: Division of vegetation type (e.g., gymnospern, dicot, fern, unknown, other, etc.)
	GrowthForm: Growth form of vegetation type (e.g., tree, shrub, vine, graminoid, etc.)
	Lifespan: Short (annual plant), Moderate (unknown, Perennial grass, etc), Long (woody plant)  

Measurement files
M1) Measurement DB emissions: Observations of vegetation biogenic emissions
	REFID: Unique identifier for publication (link to "Description References")
	VegID: Unique identifier for each vegetation type (link to "Description Vegetation")
	CAS (Chemical Abstracts Service) registry number: unique identifier for each compound (link to "Description Compounds")
	LDfrac: Fraction of emission that is light dependent
	ERuggh: Emission rate in micrograms compound per gram dry weight foliage per hour
	ERnmm2s: Emission rate in nanomoles compound per square meter foliage area (one sided) per second
	D-rating: "grade" based on number of replicates 0 (1 sample), 1 (2-4 replicates), 3: (5-11), 4: (12-23), 5 (24-99), (6 >99)
	V-rating: "grade" based on vintage of measurements 1: year 2000 or earlier, 2: year 2001 to 2010, 3: year 2011 or later
	J-rating: 0 (lowest confidence) to 4 (highest confidence) 
	position: position in canopy 1 for shade leaves, 3 for sun leaves, 2 for mixture or unknown
	Comments: additional information

M2) Measurement DB SLA: Observations of plant specific leaf area (leaf area per dry mass)
	REFID: Unique identifier for publication (link to "Description References")
	VegID: Unique identifier for each vegetation type (link to "Description Vegetation")
	SLAcm2g: Specific Leaf Area (square cm per gram dry weight)
	Position: position in canopy 1 for shade leaves, 3 for sun leaves, 2 for mixture or unknown
	J-rating: 0 (lowest confidence) to 4 (highest confidence) 

