#
#  Makefile for PCR resequencing data analysis applications 
#           Processing Objects perl modules
#  $Id$
#

CURDIR      = .
INSTALLDIR  = /usr/local/reseq/data_analysis/ProcessingObjects
INSBINDIR   = $(INSTALLDIR)/
INSLIBDIR   = $(INSTALLDIR)/lib
INSETCDIR   = $(INSTALLDIR)/etc
INSDOCDIR   = $(INSTALLDIR)/doc

install:
	@mkdir -p $(INSBINDIR)
	@cp SafeIO.pm $(INSBINDIR)/
	@cp Metrics.pm $(INSBINDIR)/
	@cp ../../primer_design/PrimerDesigner/LIMSTools/WriteXML.pm $(INSBINDIR)/
