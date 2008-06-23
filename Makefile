#
#  Makefile for PCR resequencing data analysis applications 
#           Processing Objects perl modules
#  $Id$
#

CURDIR      = .
INSTALLDIR  = /usr/local/devel/DAS/software/reseq/data_analysis/ProcessingObjects

install:
	@mkdir -p $(INSTALLDIR)
	@cp SafeIO.pm $(INSTALLDIR)/
	@cp Metrics.pm $(INSTALLDIR)/
	@cp WriteXML.pm $(INSTALLDIR)/
	@cp ReadStudyXML.pm $(INSTALLDIR)/
