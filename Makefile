#
#  Makefile for PCR resequencing data analysis applications 
#           Processing Objects perl modules
#  $Id$
#

CURDIR      = .
INSTALLDIR  = /usr/local/devel/DAS/software/reseq/data_analysis/ProcessingObjects
INSBINDIR   = $(INSTALLDIR)/
INSLIBDIR   = $(INSTALLDIR)/lib
INSETCDIR   = $(INSTALLDIR)/etc
INSDOCDIR   = $(INSTALLDIR)/doc

install:
	@mkdir -p $(INSTALLDIR)
	@mkdir -p $(INSBINDIR)
	@cp SafeIO.pm $(INSBINDIR)/
	@cp Metrics.pm $(INSBINDIR)/
	@cp WriteXML.pm $(INSBINDIR)/
	@cp ReadStudyXML.pm $(INSBINDIR)/
